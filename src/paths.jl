"""
    artifact_root()

Return the installed IgBLAST artifact root directory.
"""
function artifact_root()
    return artifact"IgBLAST"
end

"""
    igblast_prefix()

Return the `ncbi-igblast-VERSION` directory inside the artifact.
"""
function igblast_prefix()
    return joinpath(artifact_root(), "ncbi-igblast-$IGBLAST_VERSION")
end

"""
    bin_dir()

Return the directory containing IgBLAST binaries.
"""
bin_dir() = joinpath(igblast_prefix(), "bin")

"""
    executable_path(::Type{T}) where T <: AbstractIgBLAST

Absolute path to the IgBLAST executable for variant `T`.
"""
function executable_path(::Type{T}) where T <: AbstractIgBLAST
    path = joinpath(bin_dir(), executable(T))
    Sys.iswindows() && (path *= ".exe")
    return path
end

"""
    makeblastdb_path()

Absolute path to `makeblastdb`.
"""
function makeblastdb_path()
    path = joinpath(bin_dir(), "makeblastdb")
    Sys.iswindows() && (path *= ".exe")
    return path
end

"""
    set_igdata!()

Set `ENV["IGDATA"]` to the IgBLAST data directory.
"""
function set_igdata!()
    ENV["IGDATA"] = igblast_prefix()
    return nothing
end

"""
    get_igblast_url()

Download URL for the IgBLAST binary matching the current platform.
"""
function get_igblast_url()
    abi = platform_key_abi()
    platform = abi.tags["arch"]
    base_url = "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/$IGBLAST_VERSION/"

    if Sys.islinux() && contains(platform, "x86_64")
        return base_url * "ncbi-igblast-$IGBLAST_VERSION-x64-linux.tar.gz"
    elseif Sys.isapple()
        if contains(platform, "x86_64")
            return base_url * "ncbi-igblast-$IGBLAST_VERSION-x64-macosx.tar.gz"
        elseif contains(platform, "aarch64")
            error("IgBLAST does not provide a native ARM binary for macOS. Use Rosetta 2 or compile from source.")
        end
    elseif Sys.iswindows() && contains(platform, "x86_64")
        return base_url * "ncbi-igblast-$IGBLAST_VERSION-x64-win64.tar.gz"
    end
    error("Unsupported platform: $platform. IgBLAST binaries are only available for x86_64 Linux, macOS, and Windows.")
end
