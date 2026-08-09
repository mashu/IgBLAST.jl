"""
    igblastn_artifact_path(sha; artifact_path_fn=artifact_path)

Absolute path to `igblastn` inside an installed artifact identified by `sha`.
"""
function igblastn_artifact_path(sha; artifact_path_fn=artifact_path)
    return joinpath(
        artifact_path_fn(sha),
        "ncbi-igblast-$IGBLAST_VERSION",
        "bin",
        native_executable("igblastn"),
    )
end

"""
    verify_igblast_installation(sha; artifact_path_fn=artifact_path)

Ensure `igblastn` exists under artifact `sha`, otherwise throw.
"""
function verify_igblast_installation(sha; artifact_path_fn=artifact_path)
    path = igblastn_artifact_path(sha; artifact_path_fn)
    isfile(path) || error("IgBLAST installation failed. Could not find igblastn executable.")
    return path
end

"""
    is_igblast_installed()

Return `true` if the IgBLAST artifact is present and contains `igblastn`.
"""
function is_igblast_installed()
    artifact_toml = artifact_toml_path()
    haskey(Artifacts.load_artifacts_toml(artifact_toml), "IgBLAST") || return false

    igblast_hash = artifact_hash("IgBLAST", artifact_toml)
    igblast_hash === nothing && return false

    return isfile(igblastn_artifact_path(igblast_hash))
end

"""
    install_igblast(; kwargs...)

Install the IgBLAST binary artifact declared in `Artifacts.toml` via
`ensure_artifact_installed`, then verify `igblastn` is present.

Keyword hooks (for testing):
- `already_installed`
- `artifact_toml`
- `ensure_fn`
- `artifact_hash_fn`
- `artifact_path_fn`
"""
function install_igblast(;
    already_installed=is_igblast_installed,
    artifact_toml::AbstractString=artifact_toml_path(),
    ensure_fn=ensure_artifact_installed,
    artifact_hash_fn=artifact_hash,
    artifact_path_fn=artifact_path,
)
    if already_installed()
        @info "IgBLAST version $IGBLAST_VERSION is already installed."
        return nothing
    end

    haskey(Artifacts.load_artifacts_toml(artifact_toml), "IgBLAST") ||
        error("Artifacts.toml has no IgBLAST entry; cannot install.")

    @info "Installing IgBLAST version $IGBLAST_VERSION from Artifacts.toml"
    ensure_fn("IgBLAST", artifact_toml)

    sha = artifact_hash_fn("IgBLAST", artifact_toml)
    sha === nothing && error("IgBLAST artifact hash missing after install.")

    verify_igblast_installation(sha; artifact_path_fn)
    @info "IgBLAST successfully installed and verified."
    return sha
end
