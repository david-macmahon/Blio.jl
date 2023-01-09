module DataFramesGuppiRawExt

import Blio.GuppiRaw.Header

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, Tables, push!, select!
else
    import ..DataFrames: DataFrame, Tables, push!, select!
end

function DataFrame(v::Vector{Header})
    df = DataFrame(Tables.dictcolumntable(v))
    select!(df, sort(names(df)))
end

function push!(df::DataFrame, grh::Header)
    push!(df, getfield(grh, :dict), cols=:union)
end

end # module DataFramesGuppiRawExt
