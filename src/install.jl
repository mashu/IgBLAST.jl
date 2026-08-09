"""
    is_igblast_installed()

Return `true` if the IgBLAST artifact is present and contains `igblastn`.
"""
function is_igblast_installed()
    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    haskey(Artifacts.load_artifacts_toml(artifact_toml), "IgBLAST") || return false

    igblast_hash = artifact_hash("IgBLAST", artifact_toml)
    igblast_hash === nothing && return false

    igblastn_path = joinpath(
        artifact_path(igblast_hash),
        "ncbi-igblast-$IGBLAST_VERSION",
        "bin",
        "igblastn",
    )
    Sys.iswindows() && (igblastn_path *= ".exe")
    return isfile(igblastn_path)
end

"""
    install_igblast()

Download and install the IgBLAST binary artifact for the current platform.
"""
function install_igblast()
    if is_igblast_installed()
        @info "IgBLAST version $IGBLAST_VERSION is already installed."
        return nothing
    end

    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    igblast_url = get_igblast_url()

    @info "Downloading IgBLAST from $igblast_url"
    sha = add_artifact!(artifact_toml, "IgBLAST", igblast_url; force=true)
    @info "IgBLAST version $IGBLAST_VERSION installed with SHA: $sha"

    igblastn_path = joinpath(artifact_path(sha), "ncbi-igblast-$IGBLAST_VERSION", "bin", "igblastn")
    Sys.iswindows() && (igblastn_path *= ".exe")
    isfile(igblastn_path) || error("IgBLAST installation failed. Could not find igblastn executable.")

    @info "IgBLAST successfully installed and verified."
    return sha
end
