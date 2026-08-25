"""
    append_param(cmd, key, value)

Append one IgBLAST CLI argument. An empty `value` is a flag (`-key`).
Emptiness is string content, not a type, so it cannot be dispatched.
"""
append_param(cmd::Cmd, key::AbstractString, value::AbstractString) =
    append_param(cmd, key, value, Val(isempty(value)))

append_param(cmd::Cmd, key::AbstractString, ::AbstractString, ::Val{true}) =
    `$cmd -$key`

append_param(cmd::Cmd, key::AbstractString, value::AbstractString, ::Val{false}) =
    `$cmd -$key $value`

"""
    append_params(cmd, params)

Append additional IgBLAST CLI parameters.
"""
function append_params(cmd::Cmd, params::AbstractDict{<:AbstractString,<:AbstractString})
    for (key, value) in params
        cmd = append_param(cmd, key, value)
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
    build_command(exe, query, dbs, aux, output, num_threads, outfmt, params)

Build an IgBLAST `Cmd`, dispatching on the prepared database kind.
"""
function build_command(
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

function build_command(
    exe::AbstractString,
    query::AbstractString,
    dbs::PreparedV,
    aux::AbstractAuxiliary,
    output::AbstractString,
    num_threads::Integer,
    outfmt::Integer,
    additional_params::AbstractDict{<:AbstractString,<:AbstractString},
)
    cmd = `$exe -germline_db_V $(dbs.v) -query $query -outfmt $outfmt
           -num_threads $num_threads -out $output`
    cmd = apply_auxiliary(cmd, aux)
    return append_params(cmd, additional_params)
end
