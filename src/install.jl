"""
    detect_os() -> Symbol

Return `:linux`, `:apple`, `:windows`, or `:unknown` for the host.
"""
function detect_os()
    Sys.islinux() && return :linux
    Sys.isapple() && return :apple
    Sys.iswindows() && return :windows
    return :unknown
end

"""
    detect_arch() -> String

Return the host architecture tag from BinaryPlatforms.
"""
detect_arch() = platform_key_abi().tags["arch"]

"""
    platform_archive_suffix(os::Symbol, arch::AbstractString)

IgBLAST binary archive suffix for platforms supported by this package.
"""
function platform_archive_suffix(os::Symbol, arch::AbstractString)
    if os === :linux && contains(arch, "x86_64")
        return "x64-linux.tar.gz"
    elseif os === :apple && contains(arch, "x86_64")
        return "x64-macosx.tar.gz"
    elseif os === :apple && contains(arch, "aarch64")
        error("IgBLAST does not provide a native ARM binary for macOS. Use Rosetta 2 (x86_64) or compile from source.")
    elseif os === :windows
        error("Windows is not supported by IgBLAST.jl yet (Linux and macOS x86_64 only).")
    end
    error("Unsupported platform: $os/$arch. IgBLAST.jl currently supports x86_64 Linux and macOS only.")
end

"""
    platform_archive_suffix()

Archive suffix for the current host platform.
"""
platform_archive_suffix() = platform_archive_suffix(detect_os(), detect_arch())

"""
    get_igblast_url(os::Symbol, arch::AbstractString)

Download URL for the IgBLAST binary for `os`/`arch`.
"""
function get_igblast_url(os::Symbol, arch::AbstractString)
    base = "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/$IGBLAST_VERSION/"
    return string(base, "ncbi-igblast-$IGBLAST_VERSION-", platform_archive_suffix(os, arch))
end

"""
    get_igblast_url()

Download URL for the IgBLAST binary matching the current platform.
"""
get_igblast_url() = get_igblast_url(detect_os(), detect_arch())

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
