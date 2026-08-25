"""
    IgBLASTSession{T,DBs,A,E,N,O,D,W}

Prepared IgBLAST databases and run settings for variant `T`.

Created only inside [`prepare`](@ref); BLAST databases live for the `do`
block and are deleted when it returns. Call the session on a query/output
path, or [`command`](@ref) to get the raw `igblastn`/`igblastp` `Cmd`.
"""
struct IgBLASTSession{
    T<:AbstractIgBLAST,
    DBs<:AbstractPreparedDB,
    A<:AbstractAuxiliary,
    E<:AbstractString,
    N<:Integer,
    O<:Integer,
    D<:AbstractDict,
    W<:AbstractString,
}
    exe::E
    dbs::DBs
    aux::A
    num_threads::N
    outfmt::O
    additional_params::D
    workdir::W
end

function IgBLASTSession{T}(
    exe::E,
    dbs::DBs,
    aux::A,
    num_threads::N,
    outfmt::O,
    additional_params::D,
    workdir::W,
) where {
    T<:AbstractIgBLAST,
    DBs<:AbstractPreparedDB,
    A<:AbstractAuxiliary,
    E<:AbstractString,
    N<:Integer,
    O<:Integer,
    D<:AbstractDict,
    W<:AbstractString,
}
    return IgBLASTSession{T,DBs,A,E,N,O,D,W}(
        exe, dbs, aux, num_threads, outfmt, additional_params, workdir,
    )
end

"""
    command(session, query, output) -> Cmd

Build the IgBLAST command for **uncompressed** `query` and `output` paths.

This is the timing surface: `run(command(ig, query, output))` is NCBI IgBLAST
writing a plain file, with no gzip, progress monitor, or extra copies.
"""
function command(
    session::IgBLASTSession,
    query::AbstractString,
    output::AbstractString,
)
    isfile(query) || throw(ArgumentError("Query file does not exist: $query"))
    return build_command(
        session.exe,
        query,
        session.dbs,
        session.aux,
        output,
        session.num_threads,
        session.outfmt,
        session.additional_params,
    )
end

"""
    (session::IgBLASTSession)(query, output)

Run IgBLAST. Gzip query/output (`.gz`) is staged around the command; plain
paths are passed through to IgBLAST unchanged.
"""
function (session::IgBLASTSession)(query::AbstractString, output::AbstractString)
    isfile(query) || throw(ArgumentError("Query file does not exist: $query"))
    ensure_parent_directory(output)
    execute(session, query, output, file_encoding(query), file_encoding(output))
    return output
end

function execute(
    session::IgBLASTSession,
    query::AbstractString,
    output::AbstractString,
    ::PlainEncoding,
    ::PlainEncoding,
)
    run(command(session, query, output))
    return nothing
end

function execute(
    session::IgBLASTSession,
    query::AbstractString,
    output::AbstractString,
    ::GzipEncoding,
    out_enc::AbstractFileEncoding,
)
    plain_query = joinpath(session.workdir, "query.fasta")
    stage_query(query, plain_query, GzipEncoding())
    execute(session, plain_query, output, PlainEncoding(), out_enc)
    return nothing
end

function execute(
    session::IgBLASTSession,
    query::AbstractString,
    output::AbstractString,
    ::PlainEncoding,
    ::GzipEncoding,
)
    plain_output = joinpath(session.workdir, "igblast_out")
    execute(session, query, plain_output, PlainEncoding(), PlainEncoding())
    finalize_output(plain_output, output, GzipEncoding())
    return nothing
end

"""
    with_session(f, ::Type{T}, germlines; kwargs...)

Normalize keywords, then `open_session`.
"""
function with_session(
    f::F,
    ::Type{T},
    germlines::G;
    aux=noauxiliary,
    num_threads::Integer=Base.Threads.nthreads(),
    outfmt::Integer=default_outfmt(T),
    additional_params::AbstractDict=Dict{String,String}(),
) where {F, T<:AbstractIgBLAST, G<:AbstractGermlines}
    return open_session(
        f,
        T,
        germlines,
        normalize_auxiliary(aux),
        num_threads,
        outfmt,
        additional_params,
    )
end

"""
    open_session(f, ::Type{T}, germlines, aux, num_threads, outfmt, params)

Prepare databases and call `f(session)`. `aux` is a concrete
[`AbstractAuxiliary`](@ref) so the session type is known at compile time.
"""
function open_session(
    f::F,
    ::Type{T},
    germlines::G,
    aux::A,
    num_threads::Integer,
    outfmt::Integer,
    additional_params::D,
) where {
    F,
    T<:AbstractIgBLAST,
    G<:AbstractGermlines,
    A<:AbstractAuxiliary,
    D<:AbstractDict,
}
    is_igblast_installed() ||
        error("IgBLAST is not installed. Call install_igblast() or reload the package.")
    validate_inputs(germlines)
    validate_auxiliary(T, aux)
    num_threads > 0 || throw(ArgumentError("Number of threads must be positive"))

    exe = executable_path(T)
    makeblastdb = makeblastdb_path()
    isfile(exe) || error("IgBLAST executable does not exist: $exe")
    isfile(makeblastdb) || error("makeblastdb executable does not exist: $makeblastdb")

    set_igdata!()

    mktempdir() do workdir
        prepared = prepare_databases(makeblastdb, germlines, workdir, molecule(T))
        session = IgBLASTSession{T}(
            exe,
            prepared,
            aux,
            num_threads,
            outfmt,
            copy(additional_params),
            workdir,
        )
        return f(session)
    end
end

"""
    prepare(f, ::Type{T}, germlines; kwargs...)

Prepare BLAST databases, call `f(session)`, then delete the temporary
databases. Use do-block syntax:

```julia
prepare(IgBLASTn, VDJGermlines(v, d, j); additional_params=Dict("organism"=>"human")) do ig
    ig("query.fasta", "out.tsv")
end
```

Returns the value of `f`. Setup (`makeblastdb`) is outside the session body,
so `@elapsed ig(query, output)` or `@elapsed run(command(ig, query, output))`
times only IgBLAST.

`IgBLASTn` requires [`VDJGermlines`](@ref); `IgBLASTp` requires [`VGermlines`](@ref).

# Keywords
- `aux`: [`NoAuxiliary`](@ref) (default), [`AuxiliaryFile`](@ref), or a path string
- `num_threads`: IgBLAST `-num_threads` (default `Threads.nthreads()`)
- `outfmt`: IgBLAST `-outfmt` (default 19 for IgBLASTn, 7 for IgBLASTp)
- `additional_params`: extra CLI flags (`Dict`; empty value = flag)
"""
prepare(
    f::F,
    ::Type{IgBLASTn},
    germlines::VDJGermlines;
    kwargs...,
) where F = with_session(f, IgBLASTn, germlines; kwargs...)

prepare(
    f::F,
    ::Type{IgBLASTp},
    germlines::VGermlines;
    kwargs...,
) where F = with_session(f, IgBLASTp, germlines; kwargs...)

prepare(
    f::F,
    ::Type{IgBLASTn},
    v::AbstractString,
    d::AbstractString,
    j::AbstractString;
    kwargs...,
) where F = prepare(f, IgBLASTn, VDJGermlines(v, d, j); kwargs...)

prepare(
    f::F,
    ::Type{IgBLASTp},
    v::AbstractString;
    kwargs...,
) where F = prepare(f, IgBLASTp, VGermlines(v); kwargs...)
