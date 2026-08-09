"""
Plot IGHD gene usage from therapeutic antibody assignments.

This script:
1. Reads the combined IgBLAST results with D assignments
2. Filters for high-confidence assignments (identity > 90%)
3. Counts how many therapeutics use each D gene
4. Handles comma-separated ambiguous matches by counting each separately
5. Creates a barplot sorted by usage frequency
"""

using CSV
using DataFrames
using CairoMakie
using CategoricalArrays

const OUTPUT_DIR = "thera_sabdab_results"
const COMBINED_CSV = joinpath(OUTPUT_DIR, "igblast_results_combined_with_d.csv")
const OUTPUT_PLOT = joinpath(OUTPUT_DIR, "d_gene_usage_plot.pdf")

println("=" ^ 60)
println("Plotting IGHD Gene Usage")
println("=" ^ 60)

# Step 1: Read combined results
println("\n[1/3] Reading combined results...")
if !isfile(COMBINED_CSV)
    error("Combined results file with D assignments not found: $COMBINED_CSV")
end

df = CSV.read(COMBINED_CSV, DataFrame)
total_initial = nrow(df)
filter!(row -> occursin("Genetically human", string(row.Genetics)), df)
total_after_genetics = nrow(df)
println("  ✓ Loaded $total_initial assignments (after Genetics filter: $total_after_genetics)")

# Step 2: Filter for high-confidence assignments (identity > 90%)
println("\n[2/3] Filtering for assignments (Genetically human)...")

# Parse identity columns
df_parsed = transform(df,
    :identity_KIAVRA => (x -> parse.(Float64, string.(x))) => :identity_KIAVRA_parsed,
    :identity_IMGT => (x -> parse.(Float64, string.(x))) => :identity_IMGT_parsed
)

# Drop missing
df_before = nrow(df_parsed)
df_parsed = dropmissing(df_parsed, [:identity_KIAVRA_parsed, :identity_IMGT_parsed])
df_dropped = df_before - nrow(df_parsed)
total_after_missing = nrow(df_parsed)
if df_dropped > 0
    println("  ⚠ Dropped $df_dropped rows with missing identities")
end

# Keep if KIAVRA > 85% because we want reliable V assignments
df_filtered = filter(row ->
    row.identity_KIAVRA_parsed > 85.0,
    df_parsed
)
total_after_identity = nrow(df_filtered)

# Filter for mAbs only
filter!(x -> match(r"(?i)whole\s+mab", string(x.Format)) !== nothing && 
                   match(r"(?i)(mouse|feline|canine|bispecific)", x.Format) === nothing, 
                   df_filtered)
total_final = nrow(df_filtered)

# Filter out rows with missing D assignments
df_filtered = filter(row -> !ismissing(row.d_gene_assignment) && 
                              string(row.d_gene_assignment) != "", 
                     df_filtered)
total_after_d_filter = nrow(df_filtered)

println("  ✓ Filtered to $total_final high-confidence assignments")
println("    (Total: $total_after_genetics → After missing: $total_after_missing → After identity>85%: $total_after_identity → After Format filter: $total_final → After D assignment filter: $total_after_d_filter)")

# Step 3: Count D gene usage by Format (gene level, not allele level)
println("\n[3/3] Counting D gene usage by Format (gene level)...")

# Helper function to extract gene name from D assignment
# Format: "IGHD1-1*01_frame2_match3" -> "IGHD1-1"
# Also handles: "IGHD5-18*01/IGHD5-5*01_frame1_match3" -> "IGHD5-18" or "IGHD5-5"
function extract_d_gene_name(assignment_str)
    assignment_str = string(assignment_str)  # Convert to String if needed
    # Remove everything after the first underscore (frame/match info)
    gene_part = split(assignment_str, "_")[1]
    # Handle cases like "IGHD5-18*01/IGHD5-5*01" - split by "/" and take first
    if occursin("/", gene_part)
        gene_part = split(gene_part, "/")[1]
    end
    # Extract gene name (before asterisk)
    gene_name = strip(split(gene_part, "*")[1])
    return gene_name
end

# Expand rows with comma-separated D assignments
# Each comma-separated assignment should be counted separately
expanded_data = Dict[]
for row in eachrow(df_filtered)
    d_assignment = string(row.d_gene_assignment)
    
    # Split by comma to get individual assignments
    assignments = split(d_assignment, ",")
    
    for assignment in assignments
        assignment = string(strip(assignment))
        if length(assignment) > 0
            # Extract gene name
            gene_name = extract_d_gene_name(assignment)
            
            # Create a new row dictionary
            new_row = Dict(pairs(row))
            new_row[:d_gene_assignment] = assignment
            new_row[:d_gene_name] = gene_name
            push!(expanded_data, new_row)
        end
    end
end

# Create expanded dataframe
df_expanded = DataFrame(expanded_data)

# Count gene usage by both gene and Format
df_gene_format_counts = combine(groupby(df_expanded, [:d_gene_name, :Format]), nrow => :count)

# Get total counts per gene for sorting
df_gene_totals = combine(groupby(df_expanded, :d_gene_name), nrow => :total_count)
df_gene_totals = sort(df_gene_totals, :total_count, rev=false)

# Get unique Format values for stacking
format_values = sort(unique(df_gene_format_counts.Format))

println("  ✓ Found $(nrow(df_gene_totals)) unique D genes")
println("  ✓ Found $(length(format_values)) Format types: $(join(format_values, ", "))")
println("  ✓ Top 20 D genes:")
for i in 1:min(20, nrow(df_gene_totals))
    println("    $(df_gene_totals[i, :d_gene_name]): $(df_gene_totals[i, :total_count]) assignments")
end

# Step 4: Create plot using Makie
println("\n[4/4] Creating stacked barplot...")

# Get sorted gene names (by total count, descending)
gene_names = df_gene_totals.d_gene_name
num_genes = length(gene_names)

# Prepare data for stacking: for each gene, get counts for each Format
# Create arrays where each position corresponds to a gene-format combination
gene_positions = Int[]
heights = Float64[]
format_indices = Int[]

for (gene_idx, gene) in enumerate(gene_names)
    for (format_idx, format) in enumerate(format_values)
        # Find count for this gene-format combination
        matching_rows = df_gene_format_counts[
            (df_gene_format_counts.d_gene_name .== gene) .& 
            (df_gene_format_counts.Format .== format), :]
        
        count_val = nrow(matching_rows) > 0 ? matching_rows[1, :count] : 0
        
        # Include all values (including zeros) for proper stacking
        push!(gene_positions, gene_idx)
        push!(heights, Float64(count_val))
        push!(format_indices, format_idx)
    end
end

# Calculate figure height based on number of D genes (taller for more genes)
fig_height = max(2000, 10 + num_genes * 10)
fig_width = 1200  # Wider to accommodate legend

# Create figure with layout for legend
fig = Figure(size=(fig_width, fig_height))

# Create informative title with filtering details
main_title = "IGHD Gene Usage in Therapeutic Antibodies (Gene Level)"
subtitle_text = "$total_initial→ $total_after_genetics (Genetically human)→ $total_after_identity (identity>85%)→ $total_final (Whole mAb, excl. Mouse/Feline/Canine/Bispecific) → $total_after_d_filter"

ax = Axis(fig[1, 1],
    title=main_title,
    subtitle=subtitle_text,
    xlabel="Number of Assignments",
    ylabel="D Gene",
    yticks=(1:num_genes, string.(gene_names))
)

# Create stacked horizontal bar plot
# For direction=:x, positions are on y-axis, values on x-axis
# Use a categorical colormap with enough colors
num_formats = length(format_values)
colormap_name = num_formats <= 10 ? :Set3_10 : :tab20
bp = barplot!(ax, gene_positions, heights,
    direction=:x,
    stack=format_indices,
    color=format_indices,
    colormap=colormap_name
)

# Create legend with proper color mapping
format_labels = string.(format_values)
# Get colors from the colormap - normalize indices to [0, 1] range
cmap = Makie.cgrad(colormap_name)
colors_list = [Makie.to_color(cmap[(i-1)/(max(1, num_formats-1))]) 
               for i in 1:num_formats]
legend_elements = [PolyElement(color=colors_list[i]) for i in 1:num_formats]
Legend(fig[1, 2], legend_elements, format_labels, "Format", 
    tellwidth=false, tellheight=false,
    valign=:top, halign=:left)

# Add text labels showing the total count values at the end of each bar
# Position labels slightly to the right of the bar end for better visibility
max_count = maximum(df_gene_totals.total_count)
offset = max_count * 0.02  # 2% offset from bar end
for i in 1:num_genes
    total = df_gene_totals[i, :total_count]
    text!(ax, total + offset, i,
        text=string(total),
        align=(:left, :center),
        fontsize=10,
        color=:black
    )
end

# Save plot as PDF
save(OUTPUT_PLOT, fig)
println("  ✓ Saved plot to: $OUTPUT_PLOT")
println("  ✓ Plot size: $(fig_width)x$(fig_height) pixels")

println("\n" * "=" ^ 60)
println("Plotting complete!")
println("=" ^ 60)
println("\nOutput file: $OUTPUT_PLOT")
println("\nNote: Ambiguous matches (comma-separated) are counted separately,")
println("so the total number of assignments may exceed the number of therapeutics.")
