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
            error("IgBLAST does not provide a native ARM binary for macOS. You may need to use Rosetta 2 or compile from source.")
        end
    elseif Sys.iswindows() && contains(platform, "x86_64")
        return base_url * "ncbi-igblast-$IGBLAST_VERSION-x64-win64.tar.gz"
    else
        error("Unsupported platform: $platform. IgBLAST binaries are only available for x86_64 Linux, macOS, and Windows.")
    end
end

function count_fasta_sequences(filename::String)
    count = 0
    open(filename, "r") do file
        for line in eachline(file)
            if startswith(line, '>')
                count += 1
            end
        end
    end
    return count
end

"""
    prepare_db(makeblastdb::String, db_file::String, db_type::String, temp_dir::String;
              is_protein::Bool=false)

Prepare a BLAST database from a FASTA file. Can handle both nucleotide and protein sequences.

Arguments:
- makeblastdb: Path to makeblastdb executable
- db_file: Path to input FASTA file
- db_type: Type identifier for the database
- temp_dir: Directory for temporary files
- is_protein: Whether to create a protein database (will translate if input is nucleotide)

Returns:
- Path to the prepared database
"""
function prepare_db(makeblastdb::String, db_file::String, db_type::String, temp_dir::String;
                   is_protein::Bool=false)
    temp_db = joinpath(temp_dir, "$(db_type)_db.fasta")

    if !is_protein
        # For nucleotide sequences, just copy the file
        cp(db_file, temp_db, force=true)
        run(`$makeblastdb -in $temp_db -dbtype nucl -parse_seqids`)
    else
        # For protein database, translate if needed
        open(FASTA.Writer, temp_db) do writer
            reader = open(FASTA.Reader, db_file)
            for record in reader
                # Get DNA sequence and trim to be divisible by 3
                dna_seq = LongDNA{4}(sequence(record))
                trim_length = length(dna_seq) - (length(dna_seq) % 3)
                dna_seq = dna_seq[1:trim_length]
                protein_seq = translate(dna_seq)

                # Create new record with translated sequence
                new_record = FASTA.Record(FASTA.identifier(record), string(protein_seq))
                write(writer, new_record)
            end
            close(reader)
        end

        # Create protein BLAST database
        run(`$makeblastdb -in $temp_db -dbtype prot -parse_seqids`)
    end

    return temp_db
end

function copy_and_decompress(input_file::String, output_file::String)
    if endswith(input_file, ".gz")
        open(GzipDecompressorStream, input_file) do io
            write(output_file, read(io))
        end
    else
        cp(input_file, output_file)
    end
end