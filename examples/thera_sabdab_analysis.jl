"""
Example: Analyzing Therapeutic Antibody Sequences from TheraSAbDab

This example demonstrates how to:
1. Download the TheraSAbDab database
2. Process all entries (filtering for Whole mAb is commented out)
3. Extract heavy chain sequences
4. Run IgBLASTp to find V gene assignments
5. Save results in TSV format

Required packages (install with `using Pkg; Pkg.add(["CSV", "DataFrames", "Downloads"])`):
- CSV.jl
- DataFrames.jl
- Downloads.jl (standard library in Julia 1.6+)
"""

using IgBLAST
using CSV
using DataFrames
using Downloads
using FASTX

# Configuration
const THERASABDAB_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv"
const OUTPUT_DIR = "thera_sabdab_results"

# Database paths (from data folder) — IgBLASTp uses V only
const V_DB_HUMAN = "data/KI+1KGP-IGHV-SHORT.fasta"
const V_DB_MOUSE = "data/MouseIGHV.fasta"

function main(;
    v_db::String = V_DB_HUMAN,
    output_suffix::String = "",
    organism::String = "human"
)
    # Set up output file names with suffix
    suffix_str = output_suffix != "" ? "_$output_suffix" : ""
    query_fasta = joinpath(OUTPUT_DIR, "heavy_sequences$suffix_str.fasta")
    results_tsv = joinpath(OUTPUT_DIR, "igblast_results$suffix_str.tsv")
    
    println("=" ^ 60)
    println("TheraSAbDab IgBLASTp Analysis Example")
    if output_suffix != ""
        println("Database: $output_suffix")
    end
    println("=" ^ 60)
    
    # Create output directory
    mkpath(OUTPUT_DIR)
    
    # Step 1: Download TheraSAbDab database
    println("\n[1/5] Downloading TheraSAbDab database...")
    csv_file = joinpath(OUTPUT_DIR, "TheraSAbDab.csv")
    if !isfile(csv_file)
        Downloads.download(THERASABDAB_URL, csv_file)
        println("  ✓ Downloaded to: $csv_file")
    else
        println("  ✓ Using existing file: $csv_file")
    end
    
    # Step 2: Read and filter CSV
    println("\n[2/5] Reading CSV...")
    df = CSV.read(csv_file, DataFrame)
    println("  ✓ Total rows: $(nrow(df))")
    
    # Filter for Whole mAb (commented out - processing all entries)
    # filtered_df = filter(row -> row.Format == "Whole mAb", df)
    # println("  ✓ Whole mAb rows: $(nrow(filtered_df))")
    filtered_df = df  # Use all rows
    println("  ✓ Processing all rows: $(nrow(filtered_df))")
    
    # Step 3: Extract heavy sequences and save to FASTA
    println("\n[3/5] Extracting heavy chain sequences...")
    sequences_written = 0
    open(FASTA.Writer, query_fasta) do writer
        for (idx, row) in enumerate(eachrow(filtered_df))
            therapeutic = get(row, :Therapeutic, "Unknown_$idx")
            heavy_seq = get(row, :HeavySequence, "")
            
            # Skip if sequence is missing or empty
            if ismissing(heavy_seq) || heavy_seq == "" || heavy_seq == "na"
                continue
            end
            
            # Clean sequence (remove whitespace)
            heavy_seq = replace(heavy_seq, r"\s" => "")
            
            # Create FASTA record
            record = FASTA.Record("$therapeutic", heavy_seq)
            write(writer, record)
            sequences_written += 1
        end
    end
    println("  ✓ Written $sequences_written sequences to: $query_fasta")
    
    if sequences_written == 0
        error("No valid heavy chain sequences found!")
    end
    
    # Step 4: Run IgBLASTp (V germline only; aux not used)
    println("\n[4/5] Running IgBLASTp analysis...")
    println("  This may take a while depending on the number of sequences...")
    println("  Using V database: $v_db")
    
    # Install IgBLAST if needed
    if !is_igblast_installed()
        println("  Installing IgBLAST...")
        install_igblast()
    end
    
    # Run IgBLASTp (default outfmt 7; only V is prepared/used)
    igblast_output = joinpath(OUTPUT_DIR, "igblast_output$suffix_str.txt")
    run_igblast(
        IgBLASTp,
        query_fasta,
        v_db,
        igblast_output;
        additional_params = Dict("organism" => organism),
    )
    println("  ✓ IgBLASTp analysis completed")
    
    # Step 5: Convert IgBLAST output to TSV
    println("\n[5/5] Converting results to TSV format...")
    convert_igblast_to_tsv(igblast_output, results_tsv, query_fasta, v_db)
    println("  ✓ Results saved to: $results_tsv")
    
    println("\n" * "=" ^ 60)
    println("Analysis complete!")
    println("=" ^ 60)
    println("\nOutput files:")
    println("  - Query sequences: $query_fasta")
    println("  - IgBLAST output: $igblast_output")
    println("  - TSV results: $results_tsv")
end

"""
Convert IgBLAST output format 7 (tabular with comments) to TSV format.

IgBLASTp output format 7 contains:
- Comment lines starting with #
- Tab-separated data lines with columns:
  Query id, Subject id, % identity, alignment length, mismatches, 
  gap opens, q. start, q. end, s. start, s. end, e-value, bit score

Uses regex patterns to identify columns by content rather than position.
Calculates query_coverage based on query sequence lengths from the FASTA file.
Calculates database_coverage (V coverage) based on V database sequence lengths.
"""
function convert_igblast_to_tsv(igblast_output::String, tsv_output::String, query_fasta::String, v_db::String)
    results = []
    
    # Read query sequence lengths for query_coverage calculation
    query_lengths = Dict{String, Int}()
    if isfile(query_fasta)
        reader = open(FASTA.Reader, query_fasta)
        for record in reader
            query_id_fasta = FASTA.identifier(record)
            seq = FASTA.sequence(record)
            query_lengths[query_id_fasta] = length(seq)
        end
        close(reader)
        println("  ✓ Loaded query lengths for $(length(query_lengths)) sequences")
        if length(query_lengths) > 0
            sample_keys = collect(keys(query_lengths))[1:min(3, length(query_lengths))]
            println("  ✓ Sample query IDs: $sample_keys")
        end
    else
        println("  ⚠ Warning: Query FASTA file not found: $query_fasta")
    end
    
    # Read V database sequence lengths for database_coverage calculation
    # For IgBLASTp, nucleotide sequences are translated to protein, so we divide by 3
    v_lengths = Dict{String, Int}()
    if isfile(v_db)
        reader = open(FASTA.Reader, v_db)
        for record in reader
            v_id = FASTA.identifier(record)
            v_seq_nuc = FASTA.sequence(record)
            # For IgBLASTp, nucleotide sequences are translated to protein
            # Calculate protein length: trim to be divisible by 3, then divide by 3
            nuc_length = length(v_seq_nuc)
            trimmed_length = nuc_length - (nuc_length % 3)
            protein_length = trimmed_length ÷ 3
            v_lengths[v_id] = protein_length
        end
        close(reader)
        println("  ✓ Loaded V database lengths for $(length(v_lengths)) sequences")
        if length(v_lengths) > 0
            sample_keys = collect(keys(v_lengths))[1:min(3, length(v_lengths))]
            println("  ✓ Sample V gene IDs: $sample_keys")
        end
    else
        println("  ⚠ Warning: V database file not found: $v_db")
    end
    
    # Regex patterns to identify column types
    gene_pattern = r"^IGH[VDJ]|^IGK[VDJ]|^IGL[VDJ]|^TRA[VDJ]|^TRB[VDJ]|^TRG[VDJ]|^TRD[VDJ]"
    
    open(igblast_output, "r") do f
        current_query = ""
        
        for line in eachline(f)
            line = strip(line)
            isempty(line) && continue
            
            # Skip comment lines except for query info
            if startswith(line, "# Query:")
                # Extract query name from comment line like "# Query: Abagovomab"
                parts = split(line, ":", limit=2)
                if length(parts) >= 2
                    current_query = strip(parts[2])
                end
                continue
            elseif startswith(line, "#")
                continue
            end
            
            # Parse data line - split by tabs
            parts = split(line, '\t')
            query_id = current_query != "" ? current_query : "unknown"
            
            if length(parts) < 8
                continue  # Skip lines that don't have enough columns
            end
            
            # Use regex patterns to identify columns by content
            # This makes the parser robust to column order variations
            
            subject_id = ""
            identity = ""
            alignment_length = ""
            mismatches = ""
            gap_opens = ""
            gaps = ""
            query_start = ""
            query_end = ""
            subject_start = ""
            subject_end = ""
            e_value = ""
            bit_score = ""
            
            # Pattern definitions
            gene_pattern = r"^IGH[VDJ]|^IGK[VDJ]|^IGL[VDJ]|^TRA[VDJ]|^TRB[VDJ]|^TRG[VDJ]|^TRD[VDJ]"
            percentage_pattern = r"^\d+\.?\d*$"  # Number that could be a percentage
            integer_pattern = r"^\d+$"
            scientific_notation_pattern = r"^[\d\.]+[eE][\+\-]?\d+$"
            decimal_pattern = r"^[\d\.]+$"
            
            # Track which columns we've already assigned
            used_indices = Set{Int}()
            
            # 1. Find subject_id (gene name pattern)
            for (idx, part) in enumerate(parts)
                if occursin(gene_pattern, part)
                    subject_id = part
                    push!(used_indices, idx)
                    break
                end
            end
            
            if subject_id == ""
                continue  # Skip if we can't find a gene name
            end
            
            # 2. Find identity (percentage 0-100, decimal)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                val = tryparse(Float64, part)
                if val !== nothing && 0 <= val <= 100 && occursin(percentage_pattern, part)
                    identity = part
                    push!(used_indices, idx)
                    break
                end
            end
            
            # 3. Find alignment_length (integer, typically > 10)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(integer_pattern, part)
                    val = tryparse(Int, part)
                    if val !== nothing && val > 10
                        alignment_length = part
                        push!(used_indices, idx)
                        break
                    end
                end
            end
            
            # 4. Find mismatches (integer, typically < alignment_length)
            if alignment_length != ""
                aln_len = tryparse(Int, alignment_length)
                for (idx, part) in enumerate(parts)
                    if idx in used_indices
                        continue
                    end
                    if occursin(integer_pattern, part)
                        val = tryparse(Int, part)
                        if val !== nothing && val >= 0 && (aln_len === nothing || val <= aln_len)
                            mismatches = part
                            push!(used_indices, idx)
                            break
                        end
                    end
                end
            end
            
            # 5. Find gap_opens (integer, typically 0 or small)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(integer_pattern, part)
                    val = tryparse(Int, part)
                    if val !== nothing && val >= 0
                        gap_opens = part
                        push!(used_indices, idx)
                        break
                    end
                end
            end
            
            # 6. Find gaps (integer, typically 0, optional column)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(integer_pattern, part)
                    val = tryparse(Int, part)
                    if val !== nothing && val >= 0
                        gaps = part
                        push!(used_indices, idx)
                        break
                    end
                end
            end
            
            # 7-10. Find position integers (query_start, query_end, subject_start, subject_end)
            # These are typically 1-based positions
            position_fields = [query_start, query_end, subject_start, subject_end]
            position_names = [:query_start, :query_end, :subject_start, :subject_end]
            position_idx = 1
            
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(integer_pattern, part)
                    val = tryparse(Int, part)
                    if val !== nothing && val > 0 && position_idx <= length(position_fields)
                        if position_idx == 1
                            query_start = part
                        elseif position_idx == 2
                            query_end = part
                        elseif position_idx == 3
                            subject_start = part
                        elseif position_idx == 4
                            subject_end = part
                        end
                        push!(used_indices, idx)
                        position_idx += 1
                        if position_idx > 4
                            break
                        end
                    end
                end
            end
            
            # 11. Find e_value (scientific notation or very small decimal)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(scientific_notation_pattern, part) || occursin(decimal_pattern, part)
                    val = tryparse(Float64, part)
                    if val !== nothing && val >= 0
                        e_value = part
                        push!(used_indices, idx)
                        break
                    end
                end
            end
            
            # 12. Find bit_score (positive decimal/integer, typically > 50, last numeric field)
            for (idx, part) in enumerate(parts)
                if idx in used_indices
                    continue
                end
                if occursin(decimal_pattern, part)
                    val = tryparse(Float64, part)
                    if val !== nothing && val > 0
                        bit_score = part
                        push!(used_indices, idx)
                        break
                    end
                end
            end
            
            # Only add result if we have essential fields
            if subject_id != "" && identity != "" && alignment_length != ""
                # Calculate query_coverage (how much of the query sequence is covered)
                query_coverage = missing
                if query_start != "" && query_end != ""
                    if haskey(query_lengths, query_id)
                        query_len = query_lengths[query_id]
                        q_start = tryparse(Int, query_start)
                        q_end = tryparse(Int, query_end)
                        if q_start !== nothing && q_end !== nothing && query_len > 0
                            aligned_length = q_end - q_start + 1
                            query_coverage = round(aligned_length / query_len * 100, digits=2)
                            # Debug: print first few calculations
                            if length(results) < 3
                                println("    Debug query_coverage: query=$query_id, start=$q_start, end=$q_end, query_len=$query_len, query_coverage=$query_coverage%")
                            end
                        elseif length(results) < 3
                            println("    Debug: Failed to parse query positions for $query_id: start='$query_start', end='$query_end'")
                        end
                    elseif length(results) < 3
                        println("    Debug: Query ID '$query_id' not found in query_lengths. Available keys: $(collect(keys(query_lengths))[1:min(5, length(query_lengths))])")
                    end
                elseif length(results) < 3
                    println("    Debug: Missing query_start or query_end for $query_id")
                end
                
                # Calculate database_coverage (how much of the V database sequence is covered)
                # and get the translated database length
                database_coverage = missing
                database_length = missing
                
                # Get database length if we have the subject_id
                if haskey(v_lengths, subject_id)
                    database_length = v_lengths[subject_id]
                end
                
                # Calculate database_coverage if we have positions
                if subject_start != "" && subject_end != ""
                    if haskey(v_lengths, subject_id)
                        v_len = v_lengths[subject_id]
                        s_start = tryparse(Int, subject_start)
                        s_end = tryparse(Int, subject_end)
                        if s_start !== nothing && s_end !== nothing && v_len > 0
                            # Handle cases where start > end (reverse alignment)
                            aligned_length = abs(s_end - s_start) + 1
                            database_coverage = round(aligned_length / v_len * 100, digits=2)
                            # Debug: print first few calculations
                            if length(results) < 3
                                println("    Debug database_coverage: subject=$subject_id, start=$s_start, end=$s_end, v_len=$v_len, database_coverage=$database_coverage%")
                            end
                        elseif length(results) < 3
                            println("    Debug: Failed to parse subject positions for $subject_id: start='$subject_start', end='$subject_end'")
                        end
                    elseif length(results) < 3
                        println("    Debug: Subject ID '$subject_id' not found in v_lengths. Available keys: $(collect(keys(v_lengths))[1:min(5, length(v_lengths))])")
                    end
                elseif length(results) < 3
                    println("    Debug: Missing subject_start or subject_end for $subject_id")
                end
                
                push!(results, (
                        query_id = query_id,
                        subject_id = subject_id,
                        identity = identity,
                        alignment_length = alignment_length,
                        mismatches = mismatches,
                        gap_opens = gap_opens,
                        query_start = query_start,
                        query_end = query_end,
                        subject_start = subject_start,
                        subject_end = subject_end,
                        e_value = e_value,
                        bit_score = bit_score,
                        query_coverage = query_coverage,
                        database_coverage = database_coverage,
                        database_length = database_length
                    ))
            end
        end
    end
    
    # Convert to DataFrame and save as TSV
    if !isempty(results)
        df_results = DataFrame(results)
        
        # For each query, keep only the best hit (highest bit score, lowest e-value)
        # Group by query and select best match
        best_hits = combine(
            groupby(df_results, :query_id),
            x -> begin
                if nrow(x) == 1
                    return x[1, :]
                end
                
                # Parse scores and e-values, handling empty strings
                scores = Float64[]
                e_vals = Float64[]
                
                for row in eachrow(x)
                    bit_score_str = strip(string(row.bit_score))
                    e_value_str = strip(string(row.e_value))
                    
                    # Parse bit_score - should be a positive number
                    score = try
                        if bit_score_str == "" || bit_score_str == "missing" || bit_score_str == "na"
                            0.0
                        else
                            parsed = parse(Float64, bit_score_str)
                            parsed > 0 ? parsed : 0.0
                        end
                    catch
                        0.0
                    end
                    
                    # Parse e_value - can be scientific notation like "2.24e-58"
                    e_val = try
                        if e_value_str == "" || e_value_str == "missing" || e_value_str == "na"
                            Inf  # Empty e-value means worst match
                        else
                            parsed = parse(Float64, e_value_str)
                            parsed > 0 ? parsed : Inf
                        end
                    catch
                        # Try parsing scientific notation manually if needed
                        try
                            # Handle scientific notation like "2.24e-58"
                            if occursin(r"[eE]", e_value_str)
                                parse(Float64, e_value_str)
                            else
                                Inf
                            end
                        catch
                            Inf
                        end
                    end
                    
                    push!(scores, score)
                    push!(e_vals, e_val)
                end
                
                # Find index with highest bit score, breaking ties with lowest e-value
                # Use tuple comparison: (bit_score, -e_value) so higher bit_score wins,
                # and if tied, lower e_value wins
                # For e_value, lower is better, so we use -e_value in comparison
                best_idx = argmax([(s, isfinite(e) ? -e : -Inf) for (s, e) in zip(scores, e_vals)])
                return x[best_idx, :]
            end
        )
        
        CSV.write(tsv_output, best_hits, delim='\t')
        println("  ✓ Saved $(nrow(best_hits)) best V gene assignments")
    else
        # Create empty file with headers
        empty_df = DataFrame(
            query_id = String[],
            subject_id = String[],
            identity = String[],
            alignment_length = String[],
            mismatches = String[],
            gap_opens = String[],
            query_start = String[],
            query_end = String[],
            subject_start = String[],
            subject_end = String[],
            e_value = String[],
            bit_score = String[],
            query_coverage = Union{Float64, Missing}[],
            database_coverage = Union{Float64, Missing}[],
            database_length = Union{Int, Missing}[]
        )
        CSV.write(tsv_output, empty_df, delim='\t')
        println("  ⚠ Warning: No results found in IgBLAST output")
    end
end

# Run the example if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Run with human databases (default)
    println("Running analysis with human databases...")
    main()
    
    # Run with mouse databases
    println("\n\n" * "=" ^ 60)
    println("Running analysis with mouse (Mus musculus) IMGT databases...")
    println("=" ^ 60)
    main(
        v_db = V_DB_MOUSE,
        output_suffix = "MusMusculusIMGT202530-1",
        organism = "mouse"
    )
end

