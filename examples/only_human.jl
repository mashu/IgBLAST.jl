using CSV
using DataFrames
using FASTX

const OUTPUT_DIR = "thera_sabdab_results"
const KIAVRA_TSV = joinpath(OUTPUT_DIR, "igblast_results.tsv")
const THERASABDAB_CSV = joinpath(OUTPUT_DIR, "TheraSAbDab.csv")
const OUTPUT_CSV = joinpath(OUTPUT_DIR, "igblast_results_combined.csv")
const OUTPUT_CSV_FILTERED = joinpath(OUTPUT_DIR, "igblast_results_only_human_filtered.csv")

thera_df = CSV.read(THERASABDAB_CSV, DataFrame)
igblast_df = CSV.read(KIAVRA_TSV, DataFrame)
combined = leftjoin(thera_df, igblast_df, on=:Therapeutic=>:query_id)
dropmissing!(combined, :identity)

# Filter for identity > 85
filtered = filter(row->row.identity > 85, combined)
# Filter for mAbs only
filter!(x -> match(r"(?i)whole\s+mab", string(x.Format)) !== nothing && 
                   match(r"(?i)(mouse|feline|canine|bispecific)", x.Format) === nothing, 
                   filtered)
filtered = filter(row -> occursin("Genetically human", row.var"Genetics (Bispecifics delimited with semicolon)"), filtered)

# Add gene column by splitting subject_id after *
filtered.gene = [occursin("*", string(sid)) ? split(string(sid), "*")[1] : string(sid) for sid in filtered.subject_id]

# Compute SHM (Somatic Hypermutation) as 100 - identity (since identity is a percentage)
combined.shm = 100.0 .- combined.identity
filtered.shm = 100.0 .- filtered.identity

# Select only the specified columns in the requested order
filtered = select(filtered, [:Therapeutic, :subject_id, :gene, :identity, :shm, :query_coverage, :database_coverage, :Format, :Target, :Companies, :HeavySequence])

CSV.write(OUTPUT_CSV, combined, delim=',')
CSV.write(OUTPUT_CSV_FILTERED, filtered, delim=',')

using CairoMakie
using AlgebraOfGraphics

# Make a plot of the SHM histogram with both combined and filtered data
# Compute bin edges based on the combined range of both datasets
all_shm = vcat(combined.shm, filtered.shm)
min_shm = minimum(all_shm)
max_shm = maximum(all_shm)
n_bins = 20
bin_edges = range(min_shm, max_shm, length=n_bins+1)

fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel="SHM (%)", 
    ylabel="Frequency",
    title="SHM Distribution: Combined vs Filtered")
hist!(ax, combined.shm, bins=bin_edges, color=(:blue, 0.5), label="Combined")
hist!(ax, filtered.shm, bins=bin_edges, color=(:red, 0.5), label="Filtered")
vlines!(ax, 15.0, linestyle=:dash, color=:black, linewidth=2)
axislegend(ax, position=:rt)
save(joinpath(OUTPUT_DIR, "shm_histogram.pdf"), fig)
