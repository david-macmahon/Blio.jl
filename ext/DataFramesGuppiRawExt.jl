module DataFramesGuppiRawExt

import Blio.GuppiRaw.Header
import Blio.GuppiRaw.load

if isdefined(Base, :get_extension)
    import DataFrames: DataFrame, DataFrameRow, Tables, push!, select!
else
    import ..DataFrames: DataFrame, DataFrameRow, Tables, push!, select!
end

# Outer-constructor to create Header from DataFrameRow
function Header(h::DataFrameRow)
    grh = Header()
    # If backend field exists, put it first
    haskey(h, :backend)  && (grh[:backend] = h[:backend])
    for (k,v) in zip(keys(h), values(h))
        k === :backend && continue
        grh[k] = v
    end
    grh
end

function DataFrame(v::Vector{Header})
    df = DataFrame(Tables.dictcolumntable(v))
    select!(df, sort(names(df)))
end

function push!(df::DataFrame, grh::Header)
    push!(df, getfield(grh, :dict), cols=:union)
end

"""
  load(io::IO ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])
  load(fn::AbstractString ::Type{DataFrame}; datablocks=Array{<:Complex{<:Integer}}[])

Equivalent to `load(io; headers=DataFrame(), datablocks)` except that the
returned `headers` DataFrame will have column names in sorted order.  Same for
the methods taking `fn::AbstractString` and
`fns::AbstractVector{<:AbstractString}`.
"""
function load(io::IO, ::Type{DataFrame};
              headers::DataFrame = DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
    load(io; headers, datablocks)
    select!(headers, sort(names(headers))), datablocks
end

function load(rawname::AbstractString, ::Type{DataFrame};
              headers::DataFrame = DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
  open(rawname) do io
    load(io, DataFrame; headers, datablocks)
  end
end

function load(rawnames::AbstractVector{<:AbstractString}, ::Type{DataFrame};
              headers::DataFrame = DataFrame(),
              datablocks=Array{<:Complex{<:Integer}}[])
  foreach(rawnames) do rawname
    load(rawname, DataFrame; headers, datablocks)
  end
  (headers, datablocks)
end

end # module DataFramesGuppiRawExt
