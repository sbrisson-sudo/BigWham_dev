module BigWham

function flatten_c_order(matrix::AbstractMatrix)
    # Get the dimensions of the matrix
    rows, cols = size(matrix)
    # Create an empty vector to hold the flattened matrix
    flat_vector = Vector{eltype(matrix)}(undef, rows * cols)
    # Fill the vector in row-major order
    idx = 1
    for i in 1:rows
        for j in 1:cols
            flat_vector[idx] = matrix[i, j]
            idx += 1
        end
    end
    return flat_vector
end

using CxxWrap
@wrapmodule(() -> joinpath(@__DIR__, "libjl_bigwham"))

function __init__()
    @initcxx
end

using LinearMaps, LinearAlgebra

abstract type AbstractBEMatrix <: LinearMaps.LinearMap{Float64} end

struct BEMatrix <: AbstractBEMatrix
    bigwham_obj::BigWhamIOAllocated
    col_pts::Vector{Float64}
    normals::Vector{Float64}
    rotation_matrix::Vector{Float64}
    size::Dims{2}

    function BEMatrix(coor::Vector{Float64}, conn::Vector{Int32},
        kernel::String, prop::Vector{Float64}, num_threads::Int64)
        h = BigWhamIO(coor, conn, kernel, prop, num_threads)
        colpts = get_collocation_points(h)
        normals = get_element_normals(h)
        rotation_matrix = get_rotation_matrix(h)
        dim = size(colpts)[1]
        return new(h, colpts, normals, rotation_matrix, (dim, dim))
    end

end


function BEMatrix(coor::Matrix{Float64}, conn::Matrix{Int64},
    kernel::String, prop::Vector{Float64}, num_threads::Int64=8)
    coor = flatten_c_order(coor)
    conn = convert(Vector{Int32}, flatten_c_order(conn))
    return BEMatrix(coor, conn, kernel, prop, num_threads)
end

Base.size(hmat::BEMatrix) = hmat.size

function LinearMaps._unsafe_mul!(y::AbstractVector, A::BEMatrix, x::AbstractVector)
    tmp = matvec(A.bigwham_obj, x[:])
    copyto!(y, tmp)
    return tmp
end

function LinearMaps._unsafe_mul!(y::Vector{Float64}, A::BEMatrix, x::Vector{Float64})
    # tmp = matvec(A.bigwham_obj, x)
    # copyto!(y, tmp)
    matvec!(A.bigwham_obj, x, y)
    return y
end

# import Base: *
# function *(H::BEMatrix, dd::Vector{Float64})::Vector{Float64}
#     tmp = similar(dd)
#     mul!(H.bigwham_obj, dd, tmp)
#     return tmp
# end

# import Base: *
# function *(hmat::BEMatrix, dd::Vector{Float64})::Vector{Float64}
#     return matvec(hmat.bigwham_obj, dd)
# end

# function *(hmat::BEMatrix, dd::AbstractVector)::Vector{Float64}
#     return matvec(hmat.bigwham_obj, dd[:])
# end

function build_hierarchical_matrix(hmat::BEMatrix, max_leaf_size::Int64, eta::Float64, eps_aca::Float64)
    build_hierarchical_matrix(hmat.bigwham_obj, max_leaf_size, eta, eps_aca)
end


end # module BigWham
