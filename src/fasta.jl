"""
    count_fasta_sequences(filename::AbstractString)

Count FASTA records in `filename` by counting header lines.
"""
function count_fasta_sequences(filename::AbstractString)
    count = 0
    open(filename, "r") do file
        for line in eachline(file)
            startswith(line, '>') && (count += 1)
        end
    end
    return count
end

"""
    copy_and_decompress(input_file::AbstractString, output_file::AbstractString)

Copy `input_file` to `output_file`, decompressing gzip when the name ends with `.gz`.
"""
function copy_and_decompress(input_file::AbstractString, output_file::AbstractString)
    if endswith(input_file, ".gz")
        open(GzipDecompressorStream, input_file) do io
            write(output_file, read(io))
        end
    else
        cp(input_file, output_file)
    end
    return output_file
end
