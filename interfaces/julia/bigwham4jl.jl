module BigWham

using CxxWrap
@wrapmodule(() -> joinpath(@__DIR__,"libjl_bigwham"))

function __init__()
	@initcxx
end

end # module