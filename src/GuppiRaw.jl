"""
Module for interacting with GuppiRaw files.

See also:
[`GuppiRaw.Header`](@ref),
[`read!(io::IO, grh::GuppiRaw.Header)`](@ref),
[`read(io::IO, ::Type{GuppiRaw.Header})`](@ref),
[`write(io::IO, grh::GuppiRaw.Header)`](@ref),
[`Array(grh::GuppiRaw.Header, nchan::Int=0)`](@ref)
"""
module GuppiRaw

export Header

using Printf
using OrderedCollections

# Header record size
const HEADER_REC_SIZE = 80

# Maximum number of records in a header
const HEADER_MAX_RECS = 256

# Number of columns for numeric header values
const HEADER_NUMERIC_COLS = 20

# Used to detect END record
const END = map(Int, collect("END "))

"""
Type used to hold a GuppiRaw header
"""
struct Header <: AbstractDict{Symbol, Any}
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

function Base.setindex!(h::Header, val::Any, key::Symbol)
  setindex!(getfield(h, :dict), val, key)
end

function Base.setindex!(h::Header, val::Any, key::AbstractString)
  setindex!(h, val, Symbol(lowercase(key)))
end

function Base.getindex(h::Header, key::Symbol)
  getindex(getfield(h, :dict), key)
end

function Base.getindex(h::Header, key::AbstractString)
  getindex(h, Symbol(lowercase(key)))
end

function Base.get(h::Header, key::Symbol, default=nothing)
  get(getfield(h, :dict), key, default)
end

function Base.get(h::Header, key::AbstractString, default=nothing)
  get(h, Symbol(lowercase(key)), default)
end

function Base.getproperty(h::Header, sym::Symbol)
  get(h, sym, nothing)
end

function Base.length(h::Header)
  length(getfield(h, :dict))
end

function Base.propertynames(h::Header)
  Tuple(keys(getfield(h, :dict)))
end

function Base.iterate(h::Header, state...)
  iterate(getfield(h, :dict), state...)
end

function Base.empty!(h::Header)
  empty!(getfield(h, :dict))
end

"""
    read!(io::IO, grh::GuppiRaw.Header; skip_padding::Bool=true)::GuppiRaw.Header

Read a GUPPI header from `io` and populate `grh`.

If not enough bytes remain in the file, empty `grh` to indicate EOF.  Otherwise
parse header, seek `io` to the start of the data block following the header.
Always return `grh` (or throw error if GUPPI header is not found).

If `skip_padding` is true and `grh[:directio]` is non-zero, padding bytes will
be skipped after the header such that the file position will be a multiple of
512 bytes.
"""
function Base.read!(io::IO, grh::GuppiRaw.Header;
                    skip_padding::Bool=true
                   )::GuppiRaw.Header
  # Get buf from Header
  buf = getfield(grh, :buf)
  @assert size(buf, 1) == HEADER_REC_SIZE

  # Make grh empty
  empty!(grh)

  # If not enough bytes remaining (EOF), bail out
  if filesize(io) - position(io) < sizeof(buf)
    return grh
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
    elseif !isnothing(match(r"^[+-]?[0-9]+$", v))
      v = parse(Int, v)
    elseif !isnothing(tryparse(Float64, v))
      v = parse(Float64, v)
    end
    grh[k] = v
  end

  # Seek io to just after END rec
  skip(io, HEADER_REC_SIZE*endidx - sizeof(buf))

  # If directio exists and is non-zero, seek past padding
  if skip_padding && get(grh, :directio, 0) != 0
    skip(io, mod(-position(io), 512))
  end

  grh
end

"""
    read(io::IO, ::Type{GuppiRaw.Header};
         skip_padding::Bool=true)::GuppiRaw.Header

Create a `GuppiRaw::Header` object, then call `read!()` to populate it.
"""
function Base.read(io::IO, ::Type{Header}; skip_padding::Bool=true)
  read!(io, Header(), skip_padding=skip_padding)
end

# write_header_item is an internal function whose methods are used to write
# different types of values as a GUPPI raw record.
function write_header_item(io::IO, kw::Symbol, val::Integer)
  s = uppercase(rpad(kw,8)[1:8]) * "= " * lpad(string(val), HEADER_NUMERIC_COLS)
  write(io, rpad(s, 80)[1:80])
end

# Not sure the best/approriate formatting to use here,
# but this seems like a reasonable approach to start with.
function write_header_item(io::IO, kw::Symbol, val::AbstractFloat)
  valstr = @sprintf("%.16G", val)
  # If no decimal point was output, redo to force one
  if isnothing(findfirst('.', valstr))
    if isnothing(findfirst('E', valstr))
      valstr = @sprintf("%.1f", val)
    else
      valstr = @sprintf("%.1E", val)
    end
  end
  s = uppercase(rpad(kw,8)[1:8]) *  "= " * lpad(valstr, HEADER_NUMERIC_COLS)
  write(io, rpad(s, 80)[1:80])
end

function write_header_item(io::IO, kw::Symbol, val::AbstractIrrational)
  write_header_item(io, kw, Float64(val))
end

function write_header_item(io::IO, kw::Symbol, val::Rational)
  write_header_item(io, kw, Float64(val))
end

function write_header_item(io::IO, kw::Symbol, val::AbstractString)
  s = uppercase(rpad(kw,8)[1:8]) * "= '" * rpad(val, 8) * "'"
  write(io, rpad(s, 80)[1:80])
end

"""
    write(io::IO, grh::GuppiRaw.Header; skip_padding::Bool=true)

Write `grh` as a GUPPI header to `io`.  If `skip_padding` is true and
`grh[:directio]` is non-zero, padding bytes will be written after the header
such that the file position will be a multiple of 512 bytes.
"""
function Base.write(io::IO, grh::GuppiRaw.Header;
                    skip_padding::Bool=true
                   )
  bytes_written = 0

  for (k,v) in grh
    bytes_written += write_header_item(io, k, v)
  end

  bytes_written += write(io, rpad("END", 80));

  # Handle DIRECTIO padding
  if skip_padding && get(grh, :directio, 0) != 0
    padding = mod(-bytes_written, 512)
    if padding != 0
      bytes_written += write(io, zeros(UInt8, padding))
    end
  end

  return bytes_written
end

# This is a type alias for possible GuppiRaw data Arrays
RawArray = Union{Array{Complex{Int8},3},Array{Complex{Int16},3}}

"""
    Array(grh::GuppiRaw.Header, nchan::Int=0)::Array{Complex{Integer},3}

Return an uninitialized 3 dimensional Array sized for `nchan` channels of RAW
data as specified by metadata in `header`, specifically the `BLOCSIZE`,
`OBSNCHAN`, `NPOL`, and `NBITS` fields.  `nchan <= 0` implies `OBSNCHAN`.  The
data type of the Array elements will be `Complex{Int8}` when `NBITS == 8` or
`Complex{Int16}` when `NBITS == 16`.

The Array will be dimensioned as [pol, time, chan] to match the RAW data block
layout.
"""
function Base.Array(grh::Header, nchan::Int=0)::RawArray
  blocsize = grh.blocsize
  obsnchan = grh.obsnchan
  if nchan <= 0
    nchan = obsnchan
  end
  npol = get(grh, :npol, 1) < 2 ? 1 : 2
  nbits = get(grh, :nbits, 8)
  @assert nbits == 8 || nbits == 16 "unsupported nbits ($nbits)"
  eltype = nbits == 8 ? Int8 : Int16
  ntime, rem = divrem(8 * blocsize, 2 * obsnchan * npol * nbits)
  @assert rem == 0
  @assert typeof(npol) <: Int
  @assert typeof(ntime) <: Int
  Array{Complex{eltype}}(undef, npol, ntime, nchan)
end

end # module GuppiRaw
