"""
Post-processing script to combine IgBLAST assignments from different databases.

This script:
1. Reads results from human (KIAVRA) and mouse (IMGT) databases
2. Selects query_id, subject_id, identity, and e_value from each
3. Suffixes column names appropriately (KIAVRA and IMGT)
4. Left joins on both query_id and subject_id
5. Saves the combined results
"""

using CSV
using DataFrames
using FASTX

const OUTPUT_DIR = "thera_sabdab_results"
const KIAVRA_TSV = joinpath(OUTPUT_DIR, "igblast_results.tsv")
const IMGT_TSV = joinpath(OUTPUT_DIR, "igblast_results_MusMusculusIMGT202530-1.tsv")
const THERASABDAB_CSV = joinpath(OUTPUT_DIR, "TheraSAbDab.csv")
const COMBINED_CSV = joinpath(OUTPUT_DIR, "igblast_results_combined.csv")

function main()
    println("=" ^ 60)
    println("Combining IgBLAST Assignments")
    println("=" ^ 60)
    
    # Step 1: Read KIAVRA (human) results
    println("\n[1/3] Reading KIAVRA (human) assignments...")
    if !isfile(KIAVRA_TSV)
        error("KIAVRA results file not found: $KIAVRA_TSV")
    end
    
    df_kiavra = CSV.read(KIAVRA_TSV, DataFrame)
    println("  ✓ Loaded $(nrow(df_kiavra)) assignments from KIAVRA")
    println("  ✓ Columns: $(names(df_kiavra))")
    
    # Select and rename columns for KIAVRA (suffix with _KIAVRA)
    # Coverage column should already be in the TSV file
    df_kiavra_selected = select(
        df_kiavra,
        :query_id => :query_id,
        :subject_id => :subject_id_KIAVRA,
        :identity => :identity_KIAVRA,
        :coverage => :coverage_KIAVRA,
        :e_value => :e_value_KIAVRA
    )
    
    # Step 2: Read IMGT (mouse) results
    println("\n[2/3] Reading IMGT (mouse) assignments...")
    if !isfile(IMGT_TSV)
        error("IMGT results file not found: $IMGT_TSV")
    end
    
    df_imgt = CSV.read(IMGT_TSV, DataFrame)
    println("  ✓ Loaded $(nrow(df_imgt)) assignments from IMGT")
    println("  ✓ Columns: $(names(df_imgt))")
    
    # Select and rename columns for IMGT (suffix with _IMGT)
    # Coverage column should already be in the TSV file
    df_imgt_selected = select(
        df_imgt,
        :query_id => :query_id,
        :subject_id => :subject_id_IMGT,
        :identity => :identity_IMGT,
        :coverage => :coverage_IMGT,
        :e_value => :e_value_IMGT
    )
    
    # Step 3: Left join on query_id only
    println("\n[3/4] Combining assignments...")
    
    # Join on query_id only
    # This will match each query from KIAVRA with its corresponding assignment from IMGT
    # Since subject_ids will be different (human vs mouse), we join only on query_id
    # Coverage is already calculated in each TSV file during the IgBLAST run
    
    # Debug: Show sample coverage values before join
    println("  Sample KIAVRA coverage values:")
    for i in 1:min(3, nrow(df_kiavra_selected))
        row = df_kiavra_selected[i, :]
        println("    Query: $(row.query_id), Coverage: $(row.coverage_KIAVRA)")
    end
    
    println("  Sample IMGT coverage values:")
    for i in 1:min(3, nrow(df_imgt_selected))
        row = df_imgt_selected[i, :]
        println("    Query: $(row.query_id), Coverage: $(row.coverage_IMGT)")
    end
    
    # Left join on query_id
    df_combined = leftjoin(
        df_kiavra_selected,
        df_imgt_selected,
        on = :query_id
    )
    
    # Debug: Show sample combined coverage values
    println("  Sample combined coverage values after join:")
    for i in 1:min(5, nrow(df_combined))
        row = df_combined[i, :]
        println("    Query: $(row.query_id)")
        println("      KIAVRA coverage: $(row.coverage_KIAVRA)")
        println("      IMGT coverage: $(row.coverage_IMGT)")
    end
    
    # Ensure proper column order
    df_combined = select(
        df_combined,
        :query_id,
        :subject_id_KIAVRA,
        :identity_KIAVRA,
        :coverage_KIAVRA,
        :e_value_KIAVRA,
        :subject_id_IMGT,
        :identity_IMGT,
        :coverage_IMGT,
        :e_value_IMGT
    )
    
    println("  ✓ Combined $(nrow(df_combined)) rows")
    println("  ✓ KIAVRA assignments: $(sum(.!ismissing.(df_combined.subject_id_KIAVRA)))")
    println("  ✓ IMGT assignments: $(sum(.!ismissing.(df_combined.subject_id_IMGT)))")
    println("  ✓ Queries with both assignments: $(sum(.!ismissing.(df_combined.subject_id_KIAVRA) .& .!ismissing.(df_combined.subject_id_IMGT)))")
    
    # Step 5: Merge with TheraSAbDab.csv
    println("\n[5/5] Merging with TheraSAbDab metadata...")
    if !isfile(THERASABDAB_CSV)
        error("TheraSAbDab CSV file not found: $THERASABDAB_CSV")
    end
    
    df_thera = CSV.read(THERASABDAB_CSV, DataFrame)
    println("  ✓ Loaded $(nrow(df_thera)) rows from TheraSAbDab")
    
    # Find the Genetics column name (handle potential variations in column name)
    genetics_col = nothing
    for col in names(df_thera)
        if occursin("Genetics", string(col)) && occursin("Bispecifics", string(col))
            genetics_col = col
            break
        end
    end
    
    if genetics_col === nothing
        error("Genetics column not found in TheraSAbDab CSV")
    end
    
    println("  ✓ Found Genetics column: $genetics_col")
    
    # Select Therapeutic, Format, and Genetics columns
    df_thera_selected = select(
        df_thera,
        :Therapeutic => :Therapeutic,
        :Format => :Format,
        genetics_col => :Genetics
    )
    
    # Join on query_id = Therapeutic
    df_final = leftjoin(
        df_combined,
        df_thera_selected,
        on = :query_id => :Therapeutic
    )
    
    println("  ✓ Merged with TheraSAbDab metadata")
    println("  ✓ Final rows: $(nrow(df_final))")
    println("  ✓ Rows with Format: $(sum(.!ismissing.(df_final.Format)))")
    println("  ✓ Rows with Genetics: $(sum(.!ismissing.(df_final.Genetics)))")
    
    # Reorder columns: query_id, KIAVRA columns, IMGT columns, Format, Genetics (last)
    df_final = select(
        df_final,
        :query_id,
        :subject_id_KIAVRA,
        :identity_KIAVRA,
        :coverage_KIAVRA,
        :e_value_KIAVRA,
        :subject_id_IMGT,
        :identity_IMGT,
        :coverage_IMGT,
        :e_value_IMGT,
        :Format,
        :Genetics
    )
    
    # Step 6: Save combined results as CSV
    CSV.write(COMBINED_CSV, df_final, delim=',')
    println("\n  ✓ Saved combined results to: $COMBINED_CSV")
    
    println("\n" * "=" ^ 60)
    println("Combination complete!")
    println("=" ^ 60)
    println("\nOutput file: $COMBINED_CSV")
    println("\nColumn structure:")
    println("  - query_id: Query sequence identifier (Therapeutic name)")
    println("  - subject_id_KIAVRA: V gene assignment from KIAVRA (human)")
    println("  - identity_KIAVRA: Identity percentage from KIAVRA")
    println("  - coverage_KIAVRA: Query coverage percentage from KIAVRA")
    println("  - e_value_KIAVRA: E-value from KIAVRA")
    println("  - subject_id_IMGT: V gene assignment from IMGT (mouse)")
    println("  - identity_IMGT: Identity percentage from IMGT")
    println("  - coverage_IMGT: Query coverage percentage from IMGT")
    println("  - e_value_IMGT: E-value from IMGT")
    println("  - Format: Format from TheraSAbDab (e.g., Whole mAb, Fab, etc.)")
    println("  - Genetics: Genetics from TheraSAbDab (Bispecifics delimited with semicolon)")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



