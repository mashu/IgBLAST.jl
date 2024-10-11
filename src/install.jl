function is_igblast_installed()
    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")

    if !haskey(Artifacts.load_artifacts_toml(artifact_toml), "IgBLAST")
        return false
    end

    igblast_hash = artifact_hash("IgBLAST", artifact_toml)
    if igblast_hash === nothing
        return false
    end

    igblast_path = artifact_path(igblast_hash)
    igblastn_path = joinpath(igblast_path, "ncbi-igblast-$IGBLAST_VERSION", "bin", "igblastn")
    if Sys.iswindows()
        igblastn_path *= ".exe"
    end

    return isfile(igblastn_path)
end

function install_igblast()
    if is_igblast_installed()
        @info "IgBLAST version $IGBLAST_VERSION is already installed."
        return nothing
    end

    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    igblast_url = get_igblast_url()

    @info "Downloading IgBLAST from $igblast_url"
    try
        sha = add_artifact!(
            artifact_toml,
            "IgBLAST",
            igblast_url,
            force=true
        )

        @info "IgBLAST version $IGBLAST_VERSION installed with SHA: $sha"

        igblast_path = artifact_path(sha)
        igblastn_path = joinpath(igblast_path, "ncbi-igblast-$IGBLAST_VERSION", "bin", "igblastn")
        if Sys.iswindows()
            igblastn_path *= ".exe"
        end

        if !isfile(igblastn_path)
            error("IgBLAST installation failed. Could not find igblastn executable.")
        end

        @info "IgBLAST successfully installed and verified."
        return sha
    catch e
        @error "Failed to install IgBLAST" exception=(e, catch_backtrace())
        rethrow(e)
    end
end