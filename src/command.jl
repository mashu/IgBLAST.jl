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
    apply_auxiliary(cmd, supports, aux)

Attach `-auxiliary_data` when supported and an auxiliary file is present.
"""
apply_auxiliary(cmd::Cmd, ::Val{false}, ::AbstractAuxiliary) = cmd
apply_auxiliary(cmd::Cmd, ::Val{true}, ::NoAuxiliary) = cmd

function apply_auxiliary(cmd::Cmd, ::Val{true}, aux::AuxiliaryFile)
    return `$cmd -auxiliary_data $(aux.path)`
end

"""
    build_command(::Type{IgBLASTn}, ...)

Build an `igblastn` command including V/D/J germline databases.
"""
function build_command(
    ::Type{IgBLASTn},
    exe::AbstractString,
    query::AbstractString,
    v_db::AbstractString,
    d_db::AbstractString,
    j_db::AbstractString,
    aux::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer,
    outfmt::Integer,
    additional_params::AbstractDict{<:AbstractString,<:AbstractString},
)
    cmd = `$exe -germline_db_V $v_db -germline_db_D $d_db -germline_db_J $j_db
           -query $query -outfmt $outfmt -num_threads $num_threads -out $output`
    cmd = apply_auxiliary(cmd, supports_auxiliary(IgBLASTn), aux)
    return append_params(cmd, additional_params)
end

"""
    build_command(::Type{IgBLASTp}, ...)

Build an `igblastp` command. Only the V germline database is used; auxiliary data is ignored.
"""
function build_command(
    ::Type{IgBLASTp},
    exe::AbstractString,
    query::AbstractString,
    v_db::AbstractString,
    d_db::AbstractString,
    j_db::AbstractString,
    ::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer,
    outfmt::Integer,
    additional_params::AbstractDict{<:AbstractString,<:AbstractString},
)
    # IgBLASTp does not use D/J germline DBs or auxiliary_data
    cmd = `$exe -germline_db_V $v_db -query $query -outfmt 7
           -num_threads $num_threads -out $output`
    return append_params(cmd, additional_params)
end
