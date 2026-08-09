"""
    append_params(cmd, params)

Append additional IgBLAST CLI parameters. Empty values become flags.
"""
function append_params(cmd::Cmd, params::AbstractDict{<:AbstractString,<:AbstractString})
    for (key, value) in params
        cmd = isempty(value) ? `$cmd -$key` : `$cmd -$key $value`
    end
    return cmd
end

"""
    apply_auxiliary(cmd, aux)

Attach `-auxiliary_data` when a custom auxiliary file is present.
"""
apply_auxiliary(cmd::Cmd, ::NoAuxiliary) = cmd

function apply_auxiliary(cmd::Cmd, aux::AuxiliaryFile)
    return `$cmd -auxiliary_data $(aux.path)`
end

"""
    build_command(::Type{IgBLASTn}, exe, query, dbs, aux, output, num_threads, outfmt, params)

Build an `igblastn` command. Auxiliary data is optional.
"""
function build_command(
    ::Type{IgBLASTn},
    exe::AbstractString,
    query::AbstractString,
    dbs::PreparedVDJ,
    aux::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer,
    outfmt::Integer,
    additional_params::AbstractDict{<:AbstractString,<:AbstractString},
)
    cmd = `$exe -germline_db_V $(dbs.v) -germline_db_D $(dbs.d) -germline_db_J $(dbs.j)
           -query $query -outfmt $outfmt -num_threads $num_threads -out $output`
    cmd = apply_auxiliary(cmd, aux)
    return append_params(cmd, additional_params)
end

"""
    build_command(::Type{IgBLASTp}, exe, query, dbs, output, num_threads, params)

Build an `igblastp` command. Only the V germline database is used.
"""
function build_command(
    ::Type{IgBLASTp},
    exe::AbstractString,
    query::AbstractString,
    dbs::PreparedV,
    ::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer,
    ::Integer,
    additional_params::AbstractDict{<:AbstractString,<:AbstractString},
)
    outfmt = default_outfmt(IgBLASTp)
    cmd = `$exe -germline_db_V $(dbs.v) -query $query -outfmt $outfmt
           -num_threads $num_threads -out $output`
    return append_params(cmd, additional_params)
end
