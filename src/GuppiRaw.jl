"""
Module for interacting with GuppiRaw files.

See also:
[`GuppiRaw.Header`](@ref),
[`read!(io::IO, gh::GuppiRaw.Header)`](@ref),
[`data_array(gh::GuppiRaw.Header[, nchan::Int=1])`](@ref)
"""
module GuppiRaw

export Header
export read!
export data_array

using OrderedCollections
import Base.read!

# Header record size
const HEADER_REC_SIZE = 80

# Maximum number of records in a header
const HEADER_MAX_RECS = 256

# Used to detect END record
const END = map(Int, collect("END "))

"""
Type used to hold a GuppiRaw header
"""
struct Header
  # OrderedDict that holds the header key=value pairs
  dict::OrderedDict{Symbol, Any}

  # Buffer used to read in header bytes (and then some)
  buf::Array{UInt8,2}

  # Inner constructor to ensure that dict isn't set directly
  function Header(maxrecs::Integer=HEADER_MAX_RECS)::Header
    # Create an OrderedDict
    dict = OrderedDict()
    # Create an Array to buffer header data
    buf = Array{UInt8}(undef, HEADER_REC_SIZE, maxrecs);
    # Initialize fields and return
    new(dict, buf)
  end
end

function Base.setindex!(gh::Header, val::Any, key::Symbol)
  setindex!(getfield(gh, :dict), val, key)
end

function Base.setindex!(gh::Header, val::Any, key::String)
  setindex!(gh, val, Symbol(lowercase(key)))
end

function Base.getindex(gh::Header, key::Symbol)
  getindex(getfield(gh, :dict), key)
end

function Base.getindex(gh::Header, key::String)
  getindex(gh, Symbol(lowercase(key)))
end

function Base.get(gh::Header, key::Symbol, default=nothing)
  get(getfield(gh, :dict), key, default)
end

function Base.get(gh::Header, key::String, default=nothing)
  get(gh, Symbol(lowercase(key)), default)
end

function Base.getproperty(gh::Header, sym::Symbol)
  get(gh, sym, nothing)
end

function Base.length(gh::Header)
  length(getfield(gh, :dict))
end

function Base.show(gh::Header)
  show(getfield(gh, :dict))
end

function Base.display(gh::Header)
  display(getfield(gh, :dict))
end

function Base.propertynames(gh::Header)
  Tuple(keys(getfield(gh, :dict)))
end

"""
    read!(io::IO, gh::Header)::Header

Read a GUPPI header from `io` and populate `gh`.

If not enough bytes remain in the file, empty gh to indicate EOF.  Otherwise
parse header, seek `io` to the start of the data block following the header.
Always return gh (or throw error if GUPPI header is not found).
"""
function read!(io::IO, gh::Header)
  # Get fields
  dict = getfield(gh, :dict)
  buf = getfield(gh, :buf)
  @assert size(buf, 1) == HEADER_REC_SIZE

  # Make dict empty
  empty!(dict)

  # If not enough bytes remaining, bail out
  if filesize(io) - position(io) < sizeof(buf)
    return gh
  end

  # Read bytes into buf
  read!(io, buf)

  # Find END record
  endidx = findfirst(c->buf[1:4,c] == END, 1:size(buf,2))

  if isnothing(endidx)
    error("GUPPI RAW header not found in $(length(buf)) bytes")
  end

  # Parse first endidx-1 records
  for i in 1:endidx-1
    rec = String(buf[:,i])
    k, v = split(rec, '=', limit=2)
    k = Symbol(lowercase(strip(k)))
    v = strip(v)
    if v[1] == '\''
      v = strip(v, [' ', '\''])
    elseif !isnothing(tryparse(Int, v))
      v = parse(Int, v)
    elseif !isnothing(tryparse(Float64, v))
      v = parse(Float64, v)
    end
    dict[k] = v
  end

  # Seek io to just after END rec
  skip(io, HEADER_REC_SIZE*endidx - sizeof(buf))

  # If DIRECTIO exists and is non-zero, seek past padding
  if get(gh, :DIRECTIO, 0) != 0
    skip(io, mod(-position(io), 512))
  end

  gh
end

"""
    data_array(gh::GuppiRaw.Header, nchan::Int=0)::Array{Complex{Integer},3}

Return an uninitialized 3 dimensional Array sized for `nchan` channels of RAW
data as specified by metadata in `header`, specifically the `BLOCSIZE`,
`OBSNCHAN`, `NPOL`, and `NBITS` fields.  `nchan <= 0` implies `OBSNCHAN`.  The
data type of the Array elements will be `Complex{Int8}` when `NBITS == 8` or
`Complex{Int16}` when `NBITS == 16`.

The Array will be dimensioned as [pol, time, chan] to match the RAW data block
layout.
"""
function data_array(gh::Header, nchan::Int=0)::Array{Complex{Integer},3}
  blocsize = gh["BLOCSIZE"]
  obsnchan = gh["OBSNCHAN"]
  if nchan <= 0
    nchan = obsnchan
  end
  npol = get(gh, "NPOL", 1) < 2 ? 1 : 2
  nbits = get(gh, "NBITS", 8)
  @assert nbits == 8 || nbits == 16
  eltype = nbits == 8 ? Int8 : Int16
  ntime, rem = divrem(8 * blocsize, 2 * obsnchan * npol * nbits)
  @assert rem == 0
  @assert typeof(npol) <: Int
  @assert typeof(ntime) <: Int
  Array{Complex{eltype}}(undef, npol, ntime, nchan)
end

end # module GuppiRaw
