"""
    AbstractIgBLAST

Abstract supertype for IgBLAST program variants.
"""
abstract type AbstractIgBLAST end

"""
    IgBLASTn

Nucleotide IgBLAST (`igblastn`).
"""
struct IgBLASTn <: AbstractIgBLAST end

"""
    IgBLASTp

Protein IgBLAST (`igblastp`).
"""
struct IgBLASTp <: AbstractIgBLAST end

"""
    executable(::Type{<:AbstractIgBLAST})

Return the IgBLAST executable basename for the given variant.
"""
executable(::Type{IgBLASTn}) = "igblastn"
executable(::Type{IgBLASTp}) = "igblastp"

"""
    default_outfmt(::Type{<:AbstractIgBLAST})

Default `-outfmt` for the variant.
"""
default_outfmt(::Type{IgBLASTn}) = 19
default_outfmt(::Type{IgBLASTp}) = 7

# --- Molecule kind ---

"""
    AbstractMolecule

Molecule kind used when preparing BLAST databases.
"""
abstract type AbstractMolecule end

struct DNAMolecule <: AbstractMolecule end
struct ProteinMolecule <: AbstractMolecule end

"""
    molecule(::Type{<:AbstractIgBLAST})

Molecule kind required by germline databases for this variant.
"""
molecule(::Type{IgBLASTn}) = DNAMolecule()
molecule(::Type{IgBLASTp}) = ProteinMolecule()

"""
    blastdb_type(molecule) -> String

`makeblastdb -dbtype` for this molecule kind.
"""
blastdb_type(::DNAMolecule) = "nucl"
blastdb_type(::ProteinMolecule) = "prot"

# --- Auxiliary capability (trait) ---

"""
    AbstractAuxiliaryCapability

Whether an IgBLAST variant accepts `-auxiliary_data`.
"""
abstract type AbstractAuxiliaryCapability end

struct AuxiliarySupported <: AbstractAuxiliaryCapability end
struct AuxiliaryUnsupported <: AbstractAuxiliaryCapability end

"""
    auxiliary_capability(::Type{<:AbstractIgBLAST})

Auxiliary-data capability trait for the variant.
"""
auxiliary_capability(::Type{IgBLASTn}) = AuxiliarySupported()
auxiliary_capability(::Type{IgBLASTp}) = AuxiliaryUnsupported()

"""
    supports_auxiliary(::Type{<:AbstractIgBLAST}) -> Bool

Whether this variant accepts `-auxiliary_data`.
"""
supports_auxiliary(::Type{T}) where T <: AbstractIgBLAST =
    supports_auxiliary(auxiliary_capability(T))
supports_auxiliary(::AuxiliarySupported) = true
supports_auxiliary(::AuxiliaryUnsupported) = false

# --- Auxiliary data (optional) ---

"""
    AbstractAuxiliary

Optional IgBLAST auxiliary (J-gene) annotation data.
Defaults to [`NoAuxiliary`](@ref); pass [`AuxiliaryFile`](@ref) only when needed.
"""
abstract type AbstractAuxiliary end

"""
    NoAuxiliary

Sentinel meaning no auxiliary file is supplied.
"""
struct NoAuxiliary <: AbstractAuxiliary end

"""
    noauxiliary

Singleton [`NoAuxiliary`](@ref) instance.
"""
const noauxiliary = NoAuxiliary()

"""
    AuxiliaryFile{P}

Path to a custom IgBLAST auxiliary data file.
"""
struct AuxiliaryFile{P<:AbstractString} <: AbstractAuxiliary
    path::P
end

"""
    normalize_auxiliary(aux) -> AbstractAuxiliary

Normalize a path or auxiliary value to [`AbstractAuxiliary`](@ref).
"""
normalize_auxiliary(::NoAuxiliary) = NoAuxiliary()
normalize_auxiliary(aux::AuxiliaryFile) = aux
normalize_auxiliary(path::AbstractString) = AuxiliaryFile(path)

validate_inputs(::NoAuxiliary) = nothing

function validate_inputs(aux::AuxiliaryFile)
    isfile(aux.path) || throw(ArgumentError("Auxiliary file does not exist: $(aux.path)"))
    return nothing
end

validate_auxiliary(::Type{T}, aux::AbstractAuxiliary) where T <: AbstractIgBLAST =
    validate_auxiliary(auxiliary_capability(T), aux, T)

validate_auxiliary(::AuxiliarySupported, aux::AbstractAuxiliary, ::Type) = validate_inputs(aux)
validate_auxiliary(::AuxiliaryUnsupported, ::NoAuxiliary, ::Type) = nothing

function validate_auxiliary(::AuxiliaryUnsupported, ::AuxiliaryFile, ::Type{T}) where T
    throw(ArgumentError("$(executable(T)) does not accept auxiliary data"))
end

# --- Germline databases ---

"""
    AbstractGermlines

Abstract germline FASTA collection used by an IgBLAST variant.
"""
abstract type AbstractGermlines end

"""
    VGermlines{V}

V-only germline FASTA (used by [`IgBLASTp`](@ref)).
"""
struct VGermlines{V<:AbstractString} <: AbstractGermlines
    v::V
end

"""
    VDJGermlines{V,D,J}

V, D, and J germline FASTA paths (used by [`IgBLASTn`](@ref)).
"""
struct VDJGermlines{V<:AbstractString,D<:AbstractString,J<:AbstractString} <: AbstractGermlines
    v::V
    d::D
    j::J
end

function validate_inputs(db::VGermlines)
    isfile(db.v) || throw(ArgumentError("V database file does not exist: $(db.v)"))
    return nothing
end

function validate_inputs(db::VDJGermlines)
    isfile(db.v) || throw(ArgumentError("V database file does not exist: $(db.v)"))
    isfile(db.d) || throw(ArgumentError("D database file does not exist: $(db.d)"))
    isfile(db.j) || throw(ArgumentError("J database file does not exist: $(db.j)"))
    return nothing
end

"""
    germlines_for(::Type{IgBLASTn}, v, d, j)
    germlines_for(::Type{IgBLASTp}, v)

Build the germline collection appropriate for the IgBLAST variant.
"""
germlines_for(::Type{IgBLASTn}, v::AbstractString, d::AbstractString, j::AbstractString) =
    VDJGermlines(v, d, j)

germlines_for(::Type{IgBLASTp}, v::AbstractString) = VGermlines(v)

# --- Prepared BLAST databases (after makeblastdb) ---

"""
    AbstractPreparedDB

Prepared BLAST database prefixes produced for a specific IgBLAST variant.
"""
abstract type AbstractPreparedDB end

struct PreparedVDJ{V,D,J} <: AbstractPreparedDB
    v::V
    d::D
    j::J
end

struct PreparedV{V} <: AbstractPreparedDB
    v::V
end

# --- File encoding (query input / result output) ---

"""
    AbstractFileEncoding

Compression encoding inferred from a path (e.g. `.gz`).
"""
abstract type AbstractFileEncoding end

"""
    PlainEncoding

Uncompressed file encoding (default when the path does not end in `.gz`).
"""
struct PlainEncoding <: AbstractFileEncoding end

"""
    GzipEncoding

Gzip-compressed file encoding, used when a path ends in `.gz`.
"""
struct GzipEncoding <: AbstractFileEncoding end

"""
    file_encoding(path) -> AbstractFileEncoding

Reify the filename suffix as an encoding type (`*.gz` → [`GzipEncoding`](@ref)).
Callers then dispatch on that type.
"""
file_encoding(path::AbstractString) =
    gzip_suffix_encoding(Val(endswith(lowercase(path), ".gz")))

gzip_suffix_encoding(::Val{true}) = GzipEncoding()
gzip_suffix_encoding(::Val{false}) = PlainEncoding()
