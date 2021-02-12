#=
Module for copying a Filterbank file to another Filterbank file or to an HDF5 file.
The input file can be subsetted by f_start to f_stop and/or t_start to t_stop.
=#

module FbfUtil

export PgmParameters

import Base.Filesystem
import Base.get
import Printf: @printf, @sprintf
using Dates
using ArgParse
using OrderedCollections
using AstroTime

const TFMT = dateformat"HH:MM:SS.sss"
const VERSION = "0.20210206"

const Telescopes = Dict(
        0 => "Fake data",
        1 => "Arecibo",
        2 => "Ooty",
        3 => "Nancay",
        4 => "Parkes",
        5 => "Jodrell",
        6 => "GBT",
        8 => "Effelsberg",
        10 => "SRT",
        64 => "MeerKAT",
        65 => "KAT7"
        )

mutable struct PgmParameters
    infile::String
    outfile::String
    flag_debug::Bool
    freq_last::Float64
    tsamp_last::Float64
    f_start::Float64
    f_stop::Float64
    t_start::Float64
    t_stop::Float64
###    maxbufsz::Int64
    # Instantiation function:
    PgmParameters() = new()
end

include("Filterbank.jl")
include("Fb2FbCopy.jl")


validate_verb(str::AbstractString) = str == "info" || str == "copy"


function parse_commandline()
    settings = ArgParseSettings()
    settings.description = "Copy a Filterbank file to another Filterbank file or to an HDF5 file."
    @add_arg_table! settings begin
        "verb"
            help = "Command verb: info, copy."
            arg_type = String
            range_tester = validate_verb
            required = true
        "infile"
            help = "Input Filterbank file for all commands."
            arg_type = String
            required = true
        "outfile"
            help = "Output Filterbank or HDF5 file if performing a full or subset copy."
            arg_type = String
            required = false
        "--debug", "-d"
            help = "DEBUG output? (true/false)"
            action = :store_true
        "--fstart"
            help = "Start frequency in MHz.  Default: Frequency of first input file channel."
            arg_type = Float64
        "--fstop"
            help = "Stop frequency in MHz.  Default: Frequency of last input file channel."
            arg_type = Float64
        "--tstart"
            help = "Start integration elapsed time in seconds, relative to 0.  Default: 0."
            arg_type = Float64
        "--tstop"
            help = "Stop integration elapsed time in seconds, relative to 0.  Default: Last input file time integration."
            arg_type = Float64
###        "--maxbufsz", "-M"
###            help = "Maximum I/O buffer size in bytes."
###            arg_type = Int64
###            default = 100000000
    end
    parsed_args = parse_args(settings)
    return parsed_args
end


function format_angle(ddorhh::Float64)
    nanos = 0.0
    try
        nanos = Dates.Nanosecond(ddorhh * 3600.0 * 1e9)
    catch
        println("*** FbfUtil WARNING: Invalid angle detected ($ddorhh).  Using 0.0 instead.")
    end
    mytime = Dates.Time(nanos)
    return Dates.format(mytime, TFMT)
end


function file_overview(label::String, header::Filterbank.Header, parms::PgmParameters)
    println("$label header:")
    for (kw, value) in header
        if kw in Set([:src_dej, :src_raj])
            angfmt = format_angle(value)
            println("  $kw : $angfmt")
        elseif kw == :tstart
            isot = TTEpoch(value * AstroTime.days, origin=:modified_julian)
            println("  tstart (ISOT) : $isot")
            println("  tstart (MJD) : $value")
        elseif kw == :telescope_id
            telescope_name = get(Telescopes, value, "UNDEFINED")
            println("  telescope : $telescope_name ($value)")
        else
            println("  $kw : $value")
        end
    end
    if header.foff < 0  # descending values
        minfreq = parms.freq_last
        maxfreq = header.fch1
    else  # ascending values
        minfreq = header.fch1
        maxfreq = parms.freq_last
    end
    @printf("Data Shape: (%d, %d, %d)\n", header.nsamps, header.nifs, header.nchans)
    @printf("%s : %f\n", "Minimum file frequency (MHz)", minfreq)
    @printf("%s : %f\n", "Maximum file frequency (MHz)", maxfreq)

    return nothing
end


function oops(msg::AbstractString)
    println("\n*** Oops, $msg\n")
    #ccall(:jl_exit, Cvoid, (Int32,), 86) # Exit to O/S with status code 86
    exit(86)
end


function main(args)
    # Parse the command line.
    println("FbfUtil version $VERSION")
    parms = PgmParameters()
    parsed_args = parse_commandline()
    parms.flag_debug = parsed_args["debug"]
    parms.infile = abspath(parsed_args["infile"])

    # Loop through parameters.  If a value was specified, print kw & value.
    println("FbfUtil: Command line parameters:")
    for (key, value) in sort(parsed_args)
        ! isnothing(value) && println("  $key : $value")           
    end

    # Get the header of the input file.
    if ! isfile(parms.infile)
        msg = @sprintf("no such file exists: %s", parms.infile)
        oops(msg)
    end
    infh = open(parms.infile, "r")
    header::Filterbank.Header = read(infh, Filterbank.Header)
    close(infh)
    parms.freq_last = header.fch1 + (header.nchans - 1) * header.foff
    parms.tsamp_last = header.tsamp * (Float64(header.nsamps) - 1) + 0.999

    # Set up f_start.
    if isnothing(parsed_args["fstart"])
        parms.f_start = header.fch1
    else
        parms.f_start = parsed_args["fstart"]
    end

    # Set up f_stop.
    if isnothing(parsed_args["fstop"])
        parms.f_stop = header.fch1 + (header.nchans - 1) * header.foff
    else
        parms.f_stop = parsed_args["fstop"]
    end

    # Validate f_start and f_stop.
    if header.foff > 0
        if ! ( header.fch1 <= parms.f_start <= parms.f_stop <= parms.freq_last )
            msg = @sprintf("required: %f <= fstart (%f) <= fstop (%f) <= %f.", 
                           header.fch1, parms.f_start, parms.f_stop, parms.freq_last)
            oops(msg)
        end
    else #= header.foff < 0 =#
        if ! ( header.fch1 >= parms.f_start >= parms.f_stop >= parms.freq_last )
            msg = @sprintf("required: %f >= fstart (%f) >= fstop (%f) >= %f.", 
                           header.fch1, parms.f_start, parms.f_stop, parms.freq_last)
            oops(msg)
        end
    end

    # Set up t_start.
    if isnothing(parsed_args["tstart"])
        parms.t_start = 0.0
    else
        parms.t_start = parsed_args["tstart"]
    end

    # Set up t_stop.  Make sure it isn't too large.'
    if isnothing(parsed_args["tstop"])
        parms.t_stop = parms.tsamp_last
    else
        parms.t_stop = parsed_args["tstop"]
        if parms.tsamp_last <= parms.t_stop 
            @printf("*** FbfUtil WARNING: The provided --tstop parameter (%f) indicates\n", parms.t_stop)
            println("\t\t\ta time that is beyond the last time integration of the input file.")
            println("FbfUtil: Proceeding as if that parameter was not specified.")
            parms.t_stop = parms.tsamp_last
        end
    end
    
    # Validate t_start and t_stop.
    if ! ( parms.tsamp_last >= parms.t_stop >= parms.t_start >= 0.0 )
        msg = @sprintf("required: %f >= tstop (%f) >= tstart (%f) >= 0.", 
                       parms.tsamp_last, parms.t_stop, parms.t_start, )
        oops(msg)
    end
    
    # Based on command verb, dispatch to the appropriate function.
    if parsed_args["verb"] == "info"
        # Get an overview of the input file.
        file_overview("Input file", header, parms)
    else # "copy"
        # Validate nifs.
        if header.nifs != 1
            msg = @sprintf("header nifs values other than 1 are not yet supported: %d", header.nifs)
            oops(msg)
        end
        # Make sure that an output file was specified.
        isnothing(parsed_args["outfile"]) && oops("required: output file specification")
        parms.outfile = abspath(parsed_args["outfile"])
        # Validate output file extension.
        rubbish, outfile_ext = splitext(basename(parsed_args["outfile"]))
        if outfile_ext == ".h5"
            msg = @sprintf("HDF5 output files are not yet supported: %s", parsed_args["outfile"])
            oops(msg)
        end
        fb2fb_copy(header, parms)
    end

    return nothing
end


main(ARGS)


end # module fbfcopy
