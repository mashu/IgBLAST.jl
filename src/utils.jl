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

function prepare_db(makeblastdb::String, db_file::String, db_type::String, temp_dir::String)
    temp_db = joinpath(temp_dir, "$(db_type)_db.fasta")
    cp(db_file, temp_db)
    run(`$makeblastdb -in $temp_db -dbtype nucl -parse_seqids`)
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