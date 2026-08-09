"""
    native_executable(name::AbstractString)

Append `.exe` on Windows; otherwise return `name` unchanged.
"""
native_executable(name::AbstractString) = Sys.iswindows() ? string(name, ".exe") : string(name)

"""
    artifact_root()

Return the installed IgBLAST artifact root directory.
"""
artifact_root() = artifact"IgBLAST"

"""
    igblast_prefix()

Return the `ncbi-igblast-VERSION` directory inside the artifact.
"""
igblast_prefix() = joinpath(artifact_root(), "ncbi-igblast-$IGBLAST_VERSION")

"""
    bin_dir()

Return the directory containing IgBLAST binaries.
"""
bin_dir() = joinpath(igblast_prefix(), "bin")

"""
    executable_path(::Type{T}) where T <: AbstractIgBLAST

Absolute path to the IgBLAST executable for variant `T`.
"""
executable_path(::Type{T}) where T <: AbstractIgBLAST =
    joinpath(bin_dir(), native_executable(executable(T)))

"""
    makeblastdb_path()

Absolute path to `makeblastdb`.
"""
makeblastdb_path() = joinpath(bin_dir(), native_executable("makeblastdb"))

"""
    set_igdata!()

Set `ENV["IGDATA"]` to the IgBLAST data directory.
"""
function set_igdata!()
    ENV["IGDATA"] = igblast_prefix()
    return nothing
end

"""
    artifact_toml_path()

Path to the package `Artifacts.toml`.
"""
artifact_toml_path() = joinpath(@__DIR__, "..", "Artifacts.toml")
