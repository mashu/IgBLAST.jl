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

# --- Molecule kind (dispatch for database preparation) ---

"""
    AbstractMolecule

Abstract molecule kind used when preparing BLAST databases.
"""
abstract type AbstractMolecule end

struct DNAMolecule <: AbstractMolecule end
struct ProteinMolecule <: AbstractMolecule end

"""
    molecule(::Type{<:AbstractIgBLAST})

Return the molecule kind required by germline databases for this IgBLAST variant.
"""
molecule(::Type{IgBLASTn}) = DNAMolecule()
molecule(::Type{IgBLASTp}) = ProteinMolecule()

# --- Auxiliary data ---

"""
    AbstractAuxiliary

Abstract type for optional IgBLAST auxiliary (J-gene) annotation data.
"""
abstract type AbstractAuxiliary end

"""
    NoAuxiliary

Sentinel type meaning no auxiliary file is supplied.
"""
struct NoAuxiliary <: AbstractAuxiliary end

"""
    noauxiliary

Singleton instance of [`NoAuxiliary`](@ref). Pass this (or omit the aux argument)
when CDR3/auxiliary annotation is not needed.
"""
const noauxiliary = NoAuxiliary()

"""
    AuxiliaryFile{P<:AbstractString}

Path to an IgBLAST auxiliary data file.
"""
struct AuxiliaryFile{P<:AbstractString} <: AbstractAuxiliary
    path::P
end

"""
    supports_auxiliary(::Type{<:AbstractIgBLAST})

Trait: whether this IgBLAST variant accepts `-auxiliary_data`.
"""
supports_auxiliary(::Type{<:AbstractIgBLAST}) = Val{false}()
supports_auxiliary(::Type{IgBLASTn}) = Val{true}()

"""
    normalize_auxiliary(aux)

Convert `nothing`, empty string, path, or auxiliary types into `AbstractAuxiliary`.
"""
normalize_auxiliary(::Nothing) = NoAuxiliary()
normalize_auxiliary(::NoAuxiliary) = NoAuxiliary()
normalize_auxiliary(aux::AuxiliaryFile) = aux

function normalize_auxiliary(path::AbstractString)
    isempty(path) && return NoAuxiliary()
    return AuxiliaryFile(path)
end

"""
    validate_inputs(aux::AbstractAuxiliary)

Validate that an auxiliary file exists when one was requested.
"""
validate_inputs(::NoAuxiliary) = nothing

function validate_inputs(aux::AuxiliaryFile)
    isfile(aux.path) || throw(ArgumentError("Auxiliary file does not exist: $(aux.path)"))
    return nothing
end

# --- Germline databases ---

"""
    GermlineDatabases{P<:AbstractString}

V, D, and J germline FASTA paths.
"""
struct GermlineDatabases{P<:AbstractString}
    v::P
    d::P
    j::P
end

function validate_inputs(db::GermlineDatabases)
    isfile(db.v) || throw(ArgumentError("V database file does not exist: $(db.v)"))
    isfile(db.d) || throw(ArgumentError("D database file does not exist: $(db.d)"))
    isfile(db.j) || throw(ArgumentError("J database file does not exist: $(db.j)"))
    return nothing
end
