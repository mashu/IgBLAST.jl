"""
    run_igblast(
        igblast_type::Type{T},
        query_file::String, 
        v_database::String, 
        d_database::String, 
        j_database::String,
        aux_file::String,
        output_file::String,
        num_threads::Int = Base.Threads.nthreads();
        outfmt::Int = 19,
        additional_params::Dict{String, Any} = Dict()
    ) where T <: AbstractIgBLAST

Run IgBLAST with the specified parameters.

Required parameters:
- `igblast_type`: Type of IgBLAST to run (IgBLASTn or IgBLASTp)
- `query_file`: Path to the input query file
- `v_database`: Path to the V gene database
- `d_database`: Path to the D gene database
- `j_database`: Path to the J gene database
- `aux_file`: Path to the auxiliary file
- `output_file`: Path for the output file

Optional parameters:
- `num_threads`: Number of threads to use (default: all available threads)
- `outfmt`: Output format (default: 19 for AIRR format)
- `additional_params`: Dictionary of additional IgBLAST parameters

Example:
```julia
run_igblast(IgBLASTn, "query.fasta", "V.fasta", "D.fasta", "J.fasta", "aux.txt", "output.txt",
            additional_params = Dict(
                "organism" => "human",
                "domain_system" => "imgt"
            ))
```
"""
function run_igblast(
    igblast_type::Type{T},
    query_file::String, 
    v_database::String, 
    d_database::String, 
    j_database::String,
    aux_file::String,
    output_file::String,
    num_threads::Int = Base.Threads.nthreads();
    outfmt::Int = 19,
    additional_params::Dict{String, Any} = Dict{String, Any}()
) where T <: AbstractIgBLAST

    if !is_igblast_installed()
        error("IgBLAST is not installed. Please run install_igblast() first.")
    end

    @assert isfile(query_file) "Query file does not exist: $query_file"
    @assert isfile(v_database) "V database file does not exist: $v_database"
    @assert isfile(d_database) "D database file does not exist: $d_database"
    @assert isfile(j_database) "J database file does not exist: $j_database"
    @assert isfile(aux_file) "Auxiliary file does not exist: $aux_file"
    @assert num_threads > 0 "Number of threads must be positive"

    output_dir = dirname(output_file)
    if !isdir(output_dir)
        @info "Creating output directory: $output_dir"
        mkpath(output_dir)
    end

    igblast_path = artifact"IgBLAST"
    bin_path = joinpath(igblast_path, "ncbi-igblast-$IGBLAST_VERSION", "bin")
    igblast_exe = joinpath(bin_path, executable(igblast_type))
    makeblastdb = joinpath(bin_path, "makeblastdb")

    @assert isfile(igblast_exe) "IgBLAST executable does not exist: $igblast_exe"
    @assert isfile(makeblastdb) "makeblastdb executable does not exist: $makeblastdb"

    ENV["IGDATA"] = joinpath(igblast_path, "ncbi-igblast-$IGBLAST_VERSION")

    temp_dir = mktempdir()
    process = nothing
    monitor_task = nothing

    try
        temp_query = joinpath(temp_dir, "query.fasta")
        copy_and_decompress(query_file, temp_query)

        temp_v_db = prepare_db(makeblastdb, v_database, "V", temp_dir)
        temp_d_db = prepare_db(makeblastdb, d_database, "D", temp_dir)
        temp_j_db = prepare_db(makeblastdb, j_database, "J", temp_dir)

        temp_aux = joinpath(temp_dir, "aux_file")
        cp(aux_file, temp_aux)

        cmd = `$igblast_exe -germline_db_V $temp_v_db -germline_db_D $temp_d_db -germline_db_J $temp_j_db
            -auxiliary_data $temp_aux -query $temp_query -outfmt $outfmt
            -num_threads $num_threads -out $output_file`

        for (key, value) in additional_params
            cmd = `$cmd -$key $value`
        end

        total_sequences = count_fasta_sequences(temp_query)
        p = Progress(total_sequences; desc="Running IgBLAST... ", color=:green)
        done_channel = Channel{Bool}(1)

        monitor_task = @async begin
            try
                processed_sequences = 0
                last_position = 0
                while true
                    if isfile(output_file)
                        open(output_file, "r") do f
                            seekend(f)
                            if position(f) > last_position
                                seek(f, last_position)
                                for line in eachline(f)
                                    if !startswith(line, '#')
                                        processed_sequences += 1
                                        update!(p, processed_sequences)
                                    end
                                end
                                last_position = position(f)
                            end
                        end
                    end
                    if processed_sequences >= total_sequences
                        finish!(p)
                        put!(done_channel, true)
                        break
                    end
                    sleep(0.1)
                end
            catch e
                if !(e isa InterruptException)
                    @error "Error in progress monitoring task" exception=(e, catch_backtrace())
                end
            end
        end

        process = run(cmd, wait=false)

        while process_running(process)
            if isready(done_channel)
                break
            end
            sleep(0.1)
        end

        wait(process)

        if !istaskdone(monitor_task)
            schedule(monitor_task, InterruptException(), error=true)
        end

        if !success(process)
            error_output = read(stderr, String)
            @error "IgBLAST failed to run successfully" error_output
            error("IgBLAST execution failed. Check the error message above.")
        end

        @info "IgBLAST analysis completed. Output saved to $output_file"

    catch e
        if e isa InterruptException
            @warn "IgBLAST analysis interrupted by user."
        else
            @error "Error running IgBLAST" exception=(e, catch_backtrace())
        end
        rethrow(e)
    finally
        if process !== nothing && process_running(process)
            kill(process)
        end
        if monitor_task !== nothing && !istaskdone(monitor_task)
            schedule(monitor_task, InterruptException(), error=true)
        end
        rm(temp_dir, recursive=true)
    end
end