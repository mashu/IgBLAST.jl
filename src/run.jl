"""
    IgBLASTRunner{T,A}

Callable configuration for repeated IgBLAST runs of variant `T`.

Auxiliary data defaults to [`NoAuxiliary`](@ref); pass a custom path /
[`AuxiliaryFile`](@ref) only when needed.

# Example
```julia
runner = IgBLASTRunner(IgBLASTn; additional_params=Dict("organism"=>"human"))
runner(query, v, d, j, "out.tsv")

runner_aux = IgBLASTRunner(IgBLASTn; aux="human_gl.aux")
runner_aux(query, VDJGermlines(v, d, j), "out.tsv")
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
    outfmt::Union{Integer,Nothing}=nothing,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST
    normalized = normalize_auxiliary(aux)
    fmt = outfmt === nothing ? default_outfmt(T) : Int(outfmt)
    return IgBLASTRunner{T,typeof(normalized)}(
        normalized,
        Int(num_threads),
        fmt,
        additional_params,
    )
end

function (runner::IgBLASTRunner{T})(
    query::AbstractString,
    germlines::AbstractGermlines,
    output::AbstractString,
) where T <: AbstractIgBLAST
    return run_igblast(
        T,
        query,
        germlines,
        runner.aux,
        output;
        num_threads=runner.num_threads,
        outfmt=runner.outfmt,
        additional_params=runner.additional_params,
    )
end

function (runner::IgBLASTRunner{T})(
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString,
) where T <: AbstractIgBLAST
    return runner(query, germlines_for(T, v, d, j), output)
end

function (runner::IgBLASTRunner{IgBLASTp})(
    query::AbstractString,
    v::AbstractString,
    output::AbstractString,
)
    return runner(query, VGermlines(v), output)
end

"""
    run_igblast(::Type{T}, query, germlines, aux, output; kwargs...)

Core entry point. `aux` defaults to no auxiliary data; supply a custom
[`AuxiliaryFile`](@ref) (or path via convenience methods) when needed.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    germlines::AbstractGermlines,
    aux::AbstractAuxiliary,
    output::AbstractString;
    num_threads::Integer=Base.Threads.nthreads(),
    outfmt::Union{Integer,Nothing}=nothing,
    additional_params::Dict{String,String}=Dict{String,String}(),
) where T <: AbstractIgBLAST

    is_igblast_installed() || error("IgBLAST is not installed. Call install_igblast() or reload the package.")
    isfile(query) || throw(ArgumentError("Query file does not exist: $query"))
    validate_inputs(germlines)
    validate_inputs(aux)
    num_threads > 0 || throw(ArgumentError("Number of threads must be positive"))

    fmt = outfmt === nothing ? default_outfmt(T) : Int(outfmt)

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

    out_encoding = file_encoding(output)

    mktempdir() do temp_dir
        temp_query = joinpath(temp_dir, "query.fasta")
        stage_query(query, temp_query)

        prepared = prepare_databases(T, makeblastdb, germlines, temp_dir)
        staged_aux = stage_auxiliary(aux, temp_dir)

        igblast_out = igblast_output_path(output, temp_dir, out_encoding)

        cmd = build_command(
            T,
            igblast_exe,
            temp_query,
            prepared,
            staged_aux,
            igblast_out,
            num_threads,
            fmt,
            additional_params,
        )

        total_sequences = count_fasta_sequences(temp_query)
        run_process(cmd, igblast_out, total_sequences, output_format(T, fmt))
        finalize_output(igblast_out, output, out_encoding)
    end

    @info "IgBLAST analysis completed. Output saved to $output"
    return output
end

"""
    run_igblast(::Type{T}, query, germlines, output; kwargs...)

Run without an auxiliary file.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    germlines::AbstractGermlines,
    output::AbstractString;
    kwargs...,
) where T <: AbstractIgBLAST
    return run_igblast(T, query, germlines, NoAuxiliary(), output; kwargs...)
end

# --- Convenience string APIs (optional aux) ---

"""
    run_igblast(::Type{T}, query, v, d, j, output; kwargs...)

Run with V/D/J paths and **no** auxiliary file.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString;
    num_threads::Integer=Base.Threads.nthreads(),
    kwargs...,
) where T <: AbstractIgBLAST
    return run_igblast(
        T, query, germlines_for(T, v, d, j), NoAuxiliary(), output;
        num_threads=num_threads, kwargs...,
    )
end

"""
    run_igblast(::Type{T}, query, v, d, j, output, num_threads; kwargs...)

No auxiliary file, positional thread count.
"""
function run_igblast(
    ::Type{T},
    query::AbstractString,
    v::AbstractString,
    d::AbstractString,
    j::AbstractString,
    output::AbstractString,
    num_threads::Integer;
    kwargs...,
) where T <: AbstractIgBLAST
    return run_igblast(
        T, query, germlines_for(T, v, d, j), NoAuxiliary(), output;
        num_threads=num_threads, kwargs...,
    )
end

"""
    run_igblast(::Type{T}, query, v, d, j, aux, output, num_threads=...; kwargs...)

Optional auxiliary: pass a path, [`AuxiliaryFile`](@ref), [`noauxiliary`](@ref),
`nothing`, or `""` (compat).
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
    kwargs...,
) where T <: AbstractIgBLAST
    return run_igblast(
        T, query, germlines_for(T, v, d, j), normalize_auxiliary(aux), output;
        num_threads=num_threads, kwargs...,
    )
end

"""
    run_igblast(::Type{IgBLASTp}, query, v, output; kwargs...)

Protein IgBLAST with only a V germline database.
"""
function run_igblast(
    ::Type{IgBLASTp},
    query::AbstractString,
    v::AbstractString,
    output::AbstractString;
    kwargs...,
)
    return run_igblast(IgBLASTp, query, VGermlines(v), NoAuxiliary(), output; kwargs...)
end
