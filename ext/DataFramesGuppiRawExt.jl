module DataFramesGuppiRawExt

import Blio.GuppiRaw.Header

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, Tables, push!
else
    import ..DataFrames: DataFrame, Tables, push!
end

function DataFrame(v::Vector{Header})
    DataFrame(Tables.dictcolumntable(v))
end

function push!(df::DataFrame, grh::Header)
    push!(df, getfield(grh, :dict), cols=:union)
end

end # module DataFramesGuppiRawExt
