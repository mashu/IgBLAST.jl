"""
    platform_archive_suffix()

IgBLAST binary archive suffix for the current host platform.
"""
function platform_archive_suffix()
    abi = platform_key_abi()
    arch = abi.tags["arch"]

    Sys.islinux() && contains(arch, "x86_64") &&
        return "x64-linux.tar.gz"
    Sys.isapple() && contains(arch, "x86_64") &&
        return "x64-macosx.tar.gz"
    Sys.isapple() && contains(arch, "aarch64") &&
        error("IgBLAST does not provide a native ARM binary for macOS. Use Rosetta 2 or compile from source.")
    Sys.iswindows() && contains(arch, "x86_64") &&
        return "x64-win64.tar.gz"

    error("Unsupported platform: $arch. IgBLAST binaries are only available for x86_64 Linux, macOS, and Windows.")
end

"""
    get_igblast_url()

Download URL for the IgBLAST binary matching the current platform.
"""
function get_igblast_url()
    base = "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/$IGBLAST_VERSION/"
    return string(base, "ncbi-igblast-$IGBLAST_VERSION-", platform_archive_suffix())
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

    igblastn_path = joinpath(
        artifact_path(igblast_hash),
        "ncbi-igblast-$IGBLAST_VERSION",
        "bin",
        native_executable("igblastn"),
    )
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

    artifact_toml = artifact_toml_path()
    igblast_url = get_igblast_url()

    @info "Downloading IgBLAST from $igblast_url"
    sha = add_artifact!(artifact_toml, "IgBLAST", igblast_url; force=true)
    @info "IgBLAST version $IGBLAST_VERSION installed with SHA: $sha"

    igblastn_path = joinpath(
        artifact_path(sha),
        "ncbi-igblast-$IGBLAST_VERSION",
        "bin",
        native_executable("igblastn"),
    )
    isfile(igblastn_path) || error("IgBLAST installation failed. Could not find igblastn executable.")

    @info "IgBLAST successfully installed and verified."
    return sha
end
