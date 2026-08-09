"""
    IgBLASTRunner{T,A}

Callable configuration for repeated IgBLAST runs of variant `T` with auxiliary data `A`.

# Example
```julia
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner(query, v, d, j, "out.tsv")

runner_aux = IgBLASTRunner(IgBLASTn; aux="human_gl.aux")
runner_aux(query, v, d, j, "out.tsv")
```
"""
struct IgBLASTRunner{T<:AbstractIgBLAST,A<:AbstractAuxiliary}
    aux::A
    num_threads::Int
    outfmt::Int
    additional_params::Dict{String,String}
end

function IgBLASTRunner(
    ::Type{T};
    aux=NoAuxiliary(),
    num_threads::Integer=Base.Threads.nthreads(),
    outfmt::Integer=19,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST
    normalized = normalize_auxiliary(aux)
    return IgBLASTRunner{T,typeof(normalized)}(
        normalized,
        Int(num_threads),
        Int(outfmt),
        additional_params,
    )
end

function (runner::IgBLASTRunner{T})(
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString,
) where T <: AbstractIgBLAST
    return run_igblast(
        T,
        query,
        v,
        d,
        j,
        runner.aux,
        output,
        runner.num_threads;
        outfmt=runner.outfmt,
        additional_params=runner.additional_params,
    )
end

"""
    run_igblast(::Type{T}, query, v, d, j, output; kwargs...)

Run IgBLAST **without** an auxiliary file.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString;
    num_threads::Integer=Base.Threads.nthreads(),
    outfmt::Integer=19,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST
    return run_igblast(
        T, query, v, d, j, NoAuxiliary(), output, num_threads;
        outfmt=outfmt, additional_params=additional_params,
    )
end

"""
    run_igblast(::Type{T}, query, v, d, j, output, num_threads; kwargs...)

Run IgBLAST **without** an auxiliary file, with an explicit thread count.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString,
    num_threads::Integer;
    outfmt::Integer=19,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST
    return run_igblast(
        T, query, v, d, j, NoAuxiliary(), output, num_threads;
        outfmt=outfmt, additional_params=additional_params,
    )
end

"""
    run_igblast(::Type{T}, query, v, d, j, aux, output, num_threads=...; kwargs...)

Run IgBLAST with optional auxiliary data.

`aux` may be:
- omitted (use the 5-path methods above)
- `nothing` or `NoAuxiliary()` / `noauxiliary`
- an `AuxiliaryFile`
- a non-empty path `AbstractString`
- `""` (treated as no auxiliary, for backward compatibility)
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    aux::Union{AbstractAuxiliary,Nothing,AbstractString},
    output::AbstractString,
    num_threads::Integer=Base.Threads.nthreads();
    outfmt::Integer=19,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST
    return run_igblast(
        T,
        query,
        GermlineDatabases(v, d, j),
        normalize_auxiliary(aux),
        output,
        num_threads;
        outfmt=outfmt,
        additional_params=additional_params,
    )
end

"""
    run_igblast(::Type{T}, query, germlines, aux, output, num_threads=...; kwargs...)

Core method: typed germline databases and normalized auxiliary data.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    germlines::GermlineDatabases,
    aux::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer=Base.Threads.nthreads();
    outfmt::Integer=19,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST

    is_igblast_installed() || error("IgBLAST is not installed. Please run install_igblast() first.")
    isfile(query) || throw(ArgumentError("Query file does not exist: $query"))
    validate_inputs(germlines)
    validate_inputs(aux)
    num_threads > 0 || throw(ArgumentError("Number of threads must be positive"))

    output_dir = dirname(output)
    if !isempty(output_dir) && !isdir(output_dir)
        @info "Creating output directory: $output_dir"
        mkpath(output_dir)
    end

    igblast_exe = executable_path(T)
    makeblastdb = makeblastdb_path()
    isfile(igblast_exe) || error("IgBLAST executable does not exist: $igblast_exe")
    isfile(makeblastdb) || error("makeblastdb executable does not exist: $makeblastdb")

    set_igdata!()

    mktempdir() do temp_dir
        temp_query = joinpath(temp_dir, "query.fasta")
        copy_and_decompress(query, temp_query)

        v_db, d_db, j_db = prepare_databases(T, makeblastdb, germlines, temp_dir)
        staged_aux = stage_auxiliary(aux, temp_dir)

        cmd = build_command(
            T,
            igblast_exe,
            temp_query,
            v_db,
            d_db,
            j_db,
            staged_aux,
            output,
            num_threads,
            outfmt,
            additional_params,
        )

        total_sequences = count_fasta_sequences(temp_query)
        run_process(cmd, output, total_sequences)
    end

    @info "IgBLAST analysis completed. Output saved to $output"
    return output
end
