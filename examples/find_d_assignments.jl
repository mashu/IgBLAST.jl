"""
Find D gene assignments for antibody sequences.

This script:
1. Loads D gene sequences from data/KI+1KGP-IGHD-SHORT.fasta (nucleotide format)
2. Translates each D sequence in all 3 reading frames with suffixes _frame1, _frame2, _frame3
3. Loads antibody sequences from TheraSAbDab.csv (protein format)
4. Uses local alignment to find best matching D allele(s) with frame information
5. Reports frame number and number of matching amino acids
6. Only reports D if >=50% of D sequence is covered
7. Adds D assignments to combined_assignments.csv
"""

using CSV
using DataFrames
using FASTX
using BioSequences
using ProgressMeter

const OUTPUT_DIR = "thera_sabdab_results"
const D_DB = "data/KI+1KGP-IGHD-SHORT.fasta"
const THERASABDAB_CSV = joinpath(OUTPUT_DIR, "TheraSAbDab.csv")
const COMBINED_CSV = joinpath(OUTPUT_DIR, "igblast_results_combined.csv")
const OUTPUT_CSV = joinpath(OUTPUT_DIR, "igblast_results_combined_with_d.csv")

"""
Translate a DNA sequence in all 3 reading frames.
Returns a vector of 3 translated sequences (as strings).
"""
function translate_all_frames(dna_seq::LongDNA{4})
    # Ensure length is divisible by 3
    trim_length = length(dna_seq) - (length(dna_seq) % 3)
    trimmed_seq = dna_seq[1:trim_length]
    
    # Translate in frame 1 (starting at position 1)
    frame1 = translate(trimmed_seq)
    
    # Translate in frame 2 (starting at position 2)
    if length(trimmed_seq) >= 2
        frame2_seq = trimmed_seq[2:end]
        trim_frame2 = length(frame2_seq) - (length(frame2_seq) % 3)
        frame2 = translate(frame2_seq[1:trim_frame2])
    else
        frame2 = LongAA()
    end
    
    # Translate in frame 3 (starting at position 3)
    if length(trimmed_seq) >= 3
        frame3_seq = trimmed_seq[3:end]
        trim_frame3 = length(frame3_seq) - (length(frame3_seq) % 3)
        frame3 = translate(frame3_seq[1:trim_frame3])
    else
        frame3 = LongAA()
    end
    
    return [string(frame1), string(frame2), string(frame3)]
end


"""
Find the longest common substring between two strings.
This function finds the ABSOLUTE longest exact match by trying all possible
starting positions in both strings. It guarantees finding the longest possible
exact substring match.
Returns the length and the substring itself.
"""
function longest_common_substring(s1::String, s2::String)
    m, n = length(s1), length(s2)
    max_len = 0
    longest_substr = ""
    
    # Try all starting positions in s1 (antibody sequence)
    for i in 1:m
        # Try all starting positions in s2 (D sequence)
        for j in 1:n
            # Find the longest common substring starting at these positions
            # This finds ALL possible matches, ensuring we get the absolute longest
            k = 0
            while i + k <= m && j + k <= n && s1[i+k] == s2[j+k]
                k += 1
            end
            
            if k > max_len
                max_len = k
                longest_substr = s1[i:i+k-1]
            end
        end
    end
    
    return (max_len, longest_substr)
end

"""
Find the best matching D allele(s) for an antibody sequence.
Antibody sequence is protein (amino acids), D sequences are already translated with frame suffixes.
Uses exact substring matching to find longest common substrings.
Returns a tuple: (assignment_string, matching_d_sequence) or (missing, missing) if no good match.
- assignment_string: comma-separated D names with frame and match info (only best/tied matches)
- matching_d_sequence: the exact matching substring from the antibody sequence
Format: "D_name_frameX_matchN" where X is frame (1,2,3) and N is number of matching amino acids.
"""
function find_d_assignment(ab_seq_protein::String, d_sequences::Dict{String, String})
    # Antibody sequence is already protein - use it directly
    ab_seq_protein = uppercase(ab_seq_protein)
    
    # Store all matches: D_name_frame => (match_length, coverage, match_count, frame_num, d_seq, ab_match_seq)
    # ab_match_seq is the exact matching substring from the antibody sequence
    all_matches = Dict{String, Tuple{Int, Float64, Int, Int, String, String}}()
    
    # Try each D gene sequence (already translated with frame suffixes)
    # We check ALL D sequences in ALL 3 frames to find the absolute longest match
    for (d_name_frame, d_seq_protein) in d_sequences
        if length(d_seq_protein) == 0
            continue
        end
        
        d_protein_length = length(d_seq_protein)
        
        # Extract frame number from name (format: "D_name_frame1", "D_name_frame2", etc.)
        frame_num = 0
        if endswith(d_name_frame, "_frame1")
            frame_num = 1
        elseif endswith(d_name_frame, "_frame2")
            frame_num = 2
        elseif endswith(d_name_frame, "_frame3")
            frame_num = 3
        end
        
        # Find longest common substring (exact match only)
        # This finds the ABSOLUTE longest exact match between antibody and this D sequence
        # by checking all possible alignments. We then compare across all D sequences
        # and frames to find the globally longest match.
        match_length, matching_substring = longest_common_substring(ab_seq_protein, d_seq_protein)
                
        # Only consider matches that are significant
        # Require at least 50% coverage of the D sequence, OR minimum 3 amino acids
        # This ensures we don't miss longer matches just because they don't cover enough of a very long D sequence
        min_match_length = max(3, ceil(Int, d_protein_length * 0.5))
        
        if match_length >= min_match_length
            # Calculate coverage: how much of the D sequence is matched
            coverage = match_length / d_protein_length
            
            # Always store the longest match for this D (we find absolute longest first, then filter)
            if !haskey(all_matches, d_name_frame) || match_length > all_matches[d_name_frame][1]
                all_matches[d_name_frame] = (match_length, coverage, match_length, frame_num, d_seq_protein, matching_substring)
            end
        end
    end
    
    if isempty(all_matches)
        return (missing, missing)
    end
    
    # Find the best match length (longest exact match)
    max_match_length = maximum([v[1] for v in values(all_matches)])
    
    # Get all D alleles with the best match length (ties)
    # Keep D alleles paired with their matching sequences
    # Format: (formatted_name, matching_sequence)
    best_matches = Tuple{String, String}[]
    
    for (d_name_frame, (match_len, coverage, match_count, frame_num, d_seq, ab_match_seq)) in all_matches
        if match_len == max_match_length
            # Extract base D name (remove _frameX suffix)
            base_name = replace(d_name_frame, r"_frame[123]$" => "")
            formatted_name = "$(base_name)_frame$(frame_num)_match$(match_count)"
            # Keep the matching sequence paired with its D allele
            push!(best_matches, (formatted_name, ab_match_seq))
        end
    end
    
    # Sort by D name for consistent output
    sort!(best_matches, by = x -> x[1])
    
    # Extract the D allele names
    best_d_alleles = [m[1] for m in best_matches]
    
    # Get all unique matching sequences (there may be multiple different substrings of the same length)
    unique_match_seqs = unique([m[2] for m in best_matches])
    sort!(unique_match_seqs)  # Sort for consistent output
    
    # Return comma-separated strings
    assignment_str = join(best_d_alleles, ", ")
    matching_d_seq = join(unique_match_seqs, ", ")
    
    return (assignment_str, matching_d_seq)
end

function main()
    println("=" ^ 60)
    println("Finding D Gene Assignments")
    println("=" ^ 60)
    
    # Step 1: Load D gene sequences and translate in all 3 frames
    println("\n[1/4] Loading D gene sequences and translating in all 3 frames...")
    if !isfile(D_DB)
        error("D gene database not found: $D_DB")
    end
    
    d_sequences = Dict{String, String}()  # D_name_frameX => translated_protein_sequence
    reader = open(FASTA.Reader, D_DB)
    for record in reader
        d_name = FASTA.identifier(record)
        d_seq_dna = FASTA.sequence(record)
        
        # Convert to LongDNA{4}
        d_dna = try
            LongDNA{4}(uppercase(d_seq_dna))
        catch
            continue
        end
        
        # Translate in all 3 frames
        d_frames = translate_all_frames(d_dna)
        
        # Store each frame with suffix
        for (frame_idx, frame_protein) in enumerate(d_frames)
            if length(frame_protein) > 0
                frame_name = "$(d_name)_frame$(frame_idx)"
                d_sequences[frame_name] = string(frame_protein)
            end
        end
    end
    close(reader)
    println("  ✓ Loaded $(length(d_sequences)) D gene sequences (with frame variants)")
    
    # Step 2: Load antibody sequences from TheraSAbDab.csv
    println("\n[2/4] Loading antibody sequences from TheraSAbDab.csv...")
    if !isfile(THERASABDAB_CSV)
        error("TheraSAbDab CSV file not found: $THERASABDAB_CSV")
    end
    
    df_thera = CSV.read(THERASABDAB_CSV, DataFrame)
    ab_sequences = Dict{String, String}()
    
    for row in eachrow(df_thera)
        therapeutic = get(row, :Therapeutic, missing)
        heavy_seq = get(row, :HeavySequence, missing)
        
        if !ismissing(therapeutic) && !ismissing(heavy_seq) && heavy_seq != "" && heavy_seq != "na"
            # Clean sequence (remove whitespace)
            heavy_seq_clean = replace(string(heavy_seq), r"\s" => "")
            if length(heavy_seq_clean) > 0
                ab_sequences[string(therapeutic)] = heavy_seq_clean
            end
        end
    end
    
    println("  ✓ Loaded $(length(ab_sequences)) antibody sequences")
    
    # Step 3: Load combined assignments
    println("\n[3/4] Loading combined assignments...")
    if !isfile(COMBINED_CSV)
        error("Combined assignments file not found: $COMBINED_CSV")
    end
    
    df = CSV.read(COMBINED_CSV, DataFrame)
    println("  ✓ Loaded $(nrow(df)) assignments")
    
    # Remove old subject_id_D column if it exists
    if :subject_id_D in names(df)
        select!(df, Not(:subject_id_D))
        println("  ✓ Removed old subject_id_D column")
    end
    
    # Step 4: Find D assignments for each antibody
    println("\n[4/4] Finding D assignments...")
    println("  This may take a while...")
    
    d_assignments = Vector{Union{String, Missing}}(missing, nrow(df))
    d_matching_sequences = Vector{Union{String, Missing}}(missing, nrow(df))
    ab_full_sequences = Vector{Union{String, Missing}}(missing, nrow(df))
    
    pbar = Progress(nrow(df), desc="Processing antibodies...")
    for (idx, row) in enumerate(eachrow(df))
        query_id = string(row.query_id)
        
        # Get antibody sequence
        if haskey(ab_sequences, query_id)
            ab_seq = ab_sequences[query_id]
            ab_full_sequences[idx] = ab_seq
            d_assignment, d_matching_seq = find_d_assignment(ab_seq, d_sequences)
            d_assignments[idx] = d_assignment
            d_matching_sequences[idx] = d_matching_seq
        else
            d_assignments[idx] = missing
            d_matching_sequences[idx] = missing
            ab_full_sequences[idx] = missing
            if idx <= 5  # Debug: show first few missing
                println("  ⚠ Warning: No sequence found for $query_id")
            end
        end
        
        next!(pbar)
    end
    
    # Add columns with intuitive names
    df.heavy_chain_sequence = ab_full_sequences
    df.d_gene_assignment = d_assignments
    df.matching_d_sequence = d_matching_sequences
    
    # Count assignments
    num_assigned = sum(.!ismissing.(d_assignments))
    println("\n  ✓ Found D assignments for $num_assigned out of $(nrow(df)) antibodies")
    
    # Step 5: Save results
    println("\n[5/5] Saving results...")
    CSV.write(OUTPUT_CSV, df, delim=',')
    println("  ✓ Saved to: $OUTPUT_CSV")
    
    println("\n" * "=" ^ 60)
    println("D assignment complete!")
    println("=" ^ 60)
    println("\nOutput file: $OUTPUT_CSV")
    println("\nNew columns added:")
    println("  - heavy_chain_sequence: Full amino acid sequence of the antibody heavy chain")
    println("  - d_gene_assignment: D gene assignments with frame and match count (best/tied only)")
    println("  - matching_d_sequence: The actual D gene protein sequence that matched")
end

# Run if executed directly
main()
