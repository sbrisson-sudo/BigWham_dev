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
struct BEMatrix <: LinearMaps.LinearMap{Float64}
    bigwham_obj::BigWhamIO
    col_pts::Vector{Float64}
    size::Dims{2}

    function BEMatrix(coor::Vector{Float64}, conn::Vector{Int32}, kernel::String, prop::Vector{Float64})
        h = BigWhamIO(coor, conn, kernel, prop)
        colpts = get_collocation_points(h)
        dim = size(colpts)[1]
        return new(h, colpts, (dim, dim))
    end

    function BEMatrix(coor::Matrix{Float64}, conn::Matrix{Int64}, kernel::String, prop::Vector{Float64})
        coor = flatten_c_order(coor)
        conn = convert(Vector{Int32}, flatten_c_order(conn))
        return BEMatrix(coor, conn, kernel, prop)
    end
end

Base.size(hmat::BEMatrix) = hmat.size

function LinearMaps._unsafe_mul!(y, A::BEMatrix, x::AbstractVector)
    tmp = matvec(A.bigwham_obj, x[:])
    copyto!(y, tmp)
end

import Base: *
function *(hmat::BEMatrix, dd::Vector{Float64})::Vector{Float64}
    return matvec(hmat.bigwham_obj, dd)
end

function *(hmat::BEMatrix, dd::AbstractVector)::Vector{Float64}
    return matvec(hmat.bigwham_obj, dd[:])
end

function build_hierarchical_matrix(hmat::BEMatrix, max_leaf_size::Int64, eta::Float64, eps_aca::Float64)
    build_hierarchical_matrix(hmat.bigwham_obj, max_leaf_size, eta, eps_aca)
end

function get_collocation_points(hmat::BEMatrix)::Vector{Float64}
    return hmat.col_pts
end




end # module BigWham