import Base.Filesystem: filesize, isfile


function eof_watch(io::IO, context::String)
    flush(io)
    println("eof_watch: $context, output file size = ", position(io))
end


function edited_size(num_size::Int64)::String
    if num_size > 999999999
        str = @sprintf("%.9f GB", Float64(num_size) / 1.0e9)
    elseif num_size > 999999
        str = @sprintf("%.6f MB", Float64(num_size) / 1.0e6)
    elseif num_size > 999
        str = @sprintf("%.3f KB", Float64(num_size) / 1.0e3)
    else
        str = @sprintf("%d bytes", num_size)
    end
    return str
end


function fb2fb_copy(header_in::Filterbank.Header, parms::PgmParameters)
    t1 = time()
    SEC_PER_DAY = 3600.0 * 24.0
    println("fb2fb_copy: Input file path = ", parms.infile)
    println("fb2fb_copy: Output file path = ", parms.outfile)
    if parms.flag_debug
        println("fb2fb_copy DEBUG: Input data size = ", edited_size(header_in.data_size))
        @printf("fb2fb_copy DEBUG: Input data Shape = (%d, %d, %d)\n", 
                header_in.nsamps, header_in.nifs, header_in.nchans)
        @printf("fb2fb_copy DEBUG: Header size = %d bytes\n", header_in.header_size)
        @printf("%s : %f\n", "Selected Start frequency (MHz)", parms.f_start)
        @printf("%s : %f\n", "Selected Stop frequency (MHz)", parms.f_stop)
    end
    @printf("Selected Start time (in s, relative to 0) : %.3f\n", parms.t_start)
    @printf("Selected Stop time (in s, relative to 0) : %.3f\n", parms.t_stop)
   
    # Initialize input and output files.
    in_stream = open(parms.infile, "r")
    seek(in_stream, header_in.header_size)
    out_stream = open(parms.outfile, "w")

    # Convert time start/stop in e.t. seconds to time sample index values.
    ix_t_first = trunc(Int64, parms.t_start / header_in.tsamp) + 1
    ix_t_last = trunc(Int64, parms.t_stop / header_in.tsamp) + 1
    if ix_t_last < ix_t_first 
        msg = @sprintf("fb2fb_copy FATAL ix_t_last(%d) < ix_t_first (%d)", 
                        ix_t_last, ix_t_first)
        oops(msg)
    end
    if parms.flag_debug
        @printf("fb2fb_copy DEBUG: f_start: %f, f_stop: %f, f_ch1: %f, foff: %e\n", 
                parms.f_start, parms.f_stop, header_in.fch1, header_in.foff)
        println("fb2fb_copy DEBUG: ix_t_first = ", ix_t_first)
        println("fb2fb_copy DEBUG: ix_t_last = ", ix_t_last)
    end
    
    # Convert frequency start/stop to offsets within the entire spectrum (1 time sample).
    bytesize, remainder = divrem(header_in.nbits, 8) # nbits must be a multiple of 8
    if remainder != 0
        msg = @sprintf("Header.nbits = %d is invalid; myst be a multiple of 8", header_in.nbits)
        oops(msg)
    end
    spoff_f_start = bytesize * trunc(Int64, abs((parms.f_start - header_in.fch1) / header_in.foff))
    spoff_f_uprbnd = bytesize * (trunc(Int64, abs((parms.f_stop + header_in.foff - header_in.fch1) / header_in.foff)))
    sample_size_out = spoff_f_uprbnd - spoff_f_start
    if parms.flag_debug
        println("fb2fb_copy DEBUG: spoff_f_start = ", spoff_f_start)
        println("fb2fb_copy DEBUG: spoff_f_uprbnd = ", spoff_f_uprbnd)
        println("fb2fb_copy DEBUG: sample_size_out = ", sample_size_out)
    end

    # Write updated header to output stream.
    header_out::Filterbank.Header = deepcopy(header_in)
    header_out["tstart"] = header_in.tstart + Float64(ix_t_first - 1) * header_in.tsamp / SEC_PER_DAY
    header_out["nsamps"] = ix_t_last - ix_t_first + 1
    header_out["fch1"] = parms.f_start
    header_out["nchans"] = Int64(sample_size_out / bytesize)
    header_out["sample_size"] = sample_size_out
    header_out["data_size"] = header_out.sample_size * header_out.nsamps
    if header_in.tstart != header_out.tstart
        isot_in = TTEpoch(header_in.tstart * AstroTime.days, origin=:modified_julian)
        isot_out = TTEpoch(header_out.tstart * AstroTime.days, origin=:modified_julian)
        println("fb2fb_copy: Header tstart $isot_in ==> $isot_out")
    end
    if header_in.nsamps != header_out.nsamps
        println("fb2fb_copy: Header nsamps ", header_in.nsamps, " ==> ", header_out.nsamps)
    end
    if header_in.fch1 != header_out.fch1
        println("fb2fb_copy: Header fch1 ", header_in.fch1, " ==> ", header_out.fch1)
    end
    if header_in.nchans != header_out.nchans
        println("fb2fb_copy: Header nchans ", header_in.nchans, " ==> ", header_out.nchans)
    end
    if header_in.sample_size != header_out.sample_size
        println("fb2fb_copy: Header sample_size ", header_in.sample_size, " ==> ", header_out.sample_size)
    end
    if header_in.data_size != header_out.data_size
        pct::Float64 = 100.0 * Float64(header_out.data_size) / Float64(header_in.data_size)
        @printf("fb2fb_copy: Header data_size %d ==> %d (%.1f%% of original)\n", 
                header_in.data_size, header_out.data_size, pct)
    end
    parms.flag_debug && (file_overview("Output file (DEBUG)", header_out, parms))
    write(out_stream, header_out)
    parms.flag_debug && (eof_watch(out_stream, "DEBUG wrote header to output stream"))

    # For each selected time sample, copy databytes from the selected frequency range.
    buffer = Array{UInt8, 1}(undef, header_in.sample_size) # parms.maxbufsz ???
    for ndx in (1 : header_in.nsamps)
        if ndx < ix_t_first # Not yet ... skip this spectrum sample.
            skip(in_stream, header_in.sample_size)
            continue
        end
        ndx > ix_t_last && break # Done!
        # Process the current spectrum.
        spoff_f_start > 0 && skip(in_stream, spoff_f_start)
        readbytes!(in_stream, buffer, sample_size_out)
        write(out_stream, buffer)
        parms.flag_debug && (eof_watch(out_stream, "fb2fb_copy DEBUG ndx=$ndx, wrote $sample_size_out bytes"))
        next_skipsz = header_in.sample_size - spoff_f_start - sample_size_out
        next_skipsz > 0 && (skip(in_stream, next_skipsz))
    end
    
    # Close files.
    close(out_stream)
    close(in_stream)
    println("fb2fb_copy: Output file size = ", edited_size(filesize(parms.outfile)))
    @printf("fb2fb_copy: Output data Shape = (%d, %d, %d)\n", 
            header_out.nsamps, header_out.nifs, header_out.nchans)
    
    # Bye bye.
    t2 = time()
    @printf("fb2fb_copy: Elapsed time = %.2f seconds\n", t2 - t1)
    println("fb2fb_copy: End")

end
