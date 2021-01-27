### A Pluto.jl notebook ###
# v0.12.19

using Markdown
using InteractiveUtils

# ╔═╡ 656ecfa8-6066-11eb-0e06-23b2166bad85
using Pkg, DrWatson, CSV, DataFrames, HypothesisTests, Plotly, Plots, StatsPlots, Gadfly

# ╔═╡ e298dea2-6065-11eb-19f5-cd362ecd26e7
md"
# Test whether K-S test is affected by sample size
"

# ╔═╡ 44b248f0-6066-11eb-3cb1-93de8055bf5a
md"
## Setup
"

# ╔═╡ 102de98e-6098-11eb-2609-eb951f2f391d
Pkg.add("Gadfly")

# ╔═╡ d7128e60-608e-11eb-2327-b58c4d24f609
md"
### Parameters
"

# ╔═╡ dffa22ae-608e-11eb-1a88-1176d911fd14
trait_order = ["Height", "BMI", "Educational attainment", "Intelligence", "IBD", "Pigmentation"]

# ╔═╡ 0d078b2a-6066-11eb-0318-f7538d2069a4
md"
## Read in data
"

# ╔═╡ 218f2f12-6066-11eb-0aa2-85ea2dcbf4ac
data_path = "/Users/brettell/Documents/Repositories/human_traits_fst/data/20210127_results/20210127_fst.csv"

# ╔═╡ 59a5462e-6067-11eb-1cf9-e3c630e3955b
df = DataFrame(CSV.File(data_path))

# ╔═╡ 2b38df6c-6079-11eb-3cb2-f1eba304453d
md"
## Run KS test over each dataset and phenotype
"

# ╔═╡ 7c054434-607a-11eb-2d79-f99b3fd1fc78
# Filter for phenotype
begin
	traits = unique(df[!, "phenotype"])
	output = Dict(trait => [] for trait in traits)

	for (i, trait_a) in enumerate(traits)
		# filter for `trait_a`
		df_a = filter(row -> row.phenotype in [trait_a], df)
		for (j, trait_b) in enumerate(traits)
			# Filter for `trait_b`
			df_b = filter(row -> row.phenotype in [trait_b], df)
			# Run KS Test
			p_out = pvalue(ApproximateTwoSampleKSTest(df_a[!, "Fst"], df_b[!, "Fst"]))
			# Append to array
			push!(output[trait_a], p_out)
		end	
	end	
end	

# ╔═╡ c622e7ec-609a-11eb-2975-291569841606
output

# ╔═╡ 62b49708-608d-11eb-0d9f-29630c19ef7f
begin
	new = DataFrame(output)
	# reorder columns
	#new = new[!, indexin(names(new), trait_order)]
end	

# ╔═╡ 7c2a7cd8-6098-11eb-3580-7b2432ee3f53
names(new)

# ╔═╡ b3d46eb2-609a-11eb-1f6f-891036efc662
new

# ╔═╡ 2d9bbfa0-6090-11eb-1d3c-4fbe7487a928
md"
## Plot
"

# ╔═╡ 7e9a5e08-6097-11eb-2ddf-f9a88ac3a105
begin
	Plots.heatmap(1:size(data,1),
		1:size(data,2), data,
		c=cgrad([:blue, :white,:red, :yellow]),
		xlabel="x values", ylabel="y values",
		title="My title")
end

# ╔═╡ d85945c8-6097-11eb-3baf-45b17d093aa5
@df new Plots.scatter(
	:Height,
	:BMI
)

# ╔═╡ a681d820-6099-11eb-38fc-7f485e4df8c0
new_mat = convert(Matrix, new)

# ╔═╡ 3001694c-6099-11eb-0e64-470dab1f11e1
Gadfly.spy(new_mat)

# ╔═╡ 73b92fd4-609a-11eb-1814-8f261bd62c10
Gadfly.plot(new_mat, x = "IBD", Geom.histogram)

# ╔═╡ Cell order:
# ╟─e298dea2-6065-11eb-19f5-cd362ecd26e7
# ╟─44b248f0-6066-11eb-3cb1-93de8055bf5a
# ╠═656ecfa8-6066-11eb-0e06-23b2166bad85
# ╠═102de98e-6098-11eb-2609-eb951f2f391d
# ╠═d7128e60-608e-11eb-2327-b58c4d24f609
# ╠═dffa22ae-608e-11eb-1a88-1176d911fd14
# ╟─0d078b2a-6066-11eb-0318-f7538d2069a4
# ╠═218f2f12-6066-11eb-0aa2-85ea2dcbf4ac
# ╠═59a5462e-6067-11eb-1cf9-e3c630e3955b
# ╟─2b38df6c-6079-11eb-3cb2-f1eba304453d
# ╠═7c054434-607a-11eb-2d79-f99b3fd1fc78
# ╠═c622e7ec-609a-11eb-2975-291569841606
# ╠═62b49708-608d-11eb-0d9f-29630c19ef7f
# ╠═7c2a7cd8-6098-11eb-3580-7b2432ee3f53
# ╠═b3d46eb2-609a-11eb-1f6f-891036efc662
# ╠═2d9bbfa0-6090-11eb-1d3c-4fbe7487a928
# ╠═7e9a5e08-6097-11eb-2ddf-f9a88ac3a105
# ╠═d85945c8-6097-11eb-3baf-45b17d093aa5
# ╠═a681d820-6099-11eb-38fc-7f485e4df8c0
# ╠═3001694c-6099-11eb-0e64-470dab1f11e1
# ╠═73b92fd4-609a-11eb-1814-8f261bd62c10
