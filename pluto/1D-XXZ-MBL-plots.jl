### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ d05f72ed-31ac-4983-b510-ce6d6527879a
html"""
<style>main { max-width: 40%;}</style>
"""

# ╔═╡ a35a6fb5-dd1e-4998-8c2a-097acb9bf7d4
tempdir = mkdir(tempname())

# ╔═╡ 1516ba16-a8ef-11ec-3979-1f7edd32dc6d
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate(tempdir)
	Pkg.add(["Plots", "PyPlot", "Statistics", "LinearAlgebra","LaTeXStrings","LsqFit","SpinSymmetry","PlutoUI"])
	Pkg.add(url="https://github.com/pablosanjose/FilteredMatrices.jl")
	Pkg.add(path=expanduser("~/julia/XXZNumerics"))
	Pkg.add(path=expanduser("~/julia/SimLib"))
	using Plots,Statistics,LinearAlgebra,SimLib,LaTeXStrings,LsqFit,XXZNumerics,SpinSymmetry,PlutoUI
	pyplot()
	TableOfContents()
end

# ╔═╡ fc5ab3d2-acdd-47a7-bbdc-c40b28976907
begin
	DATAPATH = expanduser("~/results/xxzpaper")
	HIGH_DENSITY = joinpath(DATAPATH, "high_density-avg/merged")
	LOW_DENSITY = joinpath(DATAPATH, "low_density-avg/merged")
	SAVEPATH = expanduser("~/Documents/PhD/my-papers/1D-XXZ-MBL/gfx2")
end

# ╔═╡ f92f65b3-59f0-4067-a744-3105c9411d0a
begin
	function center_region(arr, percent; dim=1)
		N = size(arr, dim)
		cutoff = div(round(Int, N*(1-percent)), 2)
		region = 1+cutoff:N-cutoff
		return view(arr, axes(arr)[1:dim-1]..., region, axes(arr)[1+dim:end]...)
	end

	maskedmean(A, mask; dims) = sum(A .* mask;dims) ./ sum(mask; dims)
	maskedvar(A, mask; dims) =
		maskedmean((A .- maskedmean(A, mask; dims)) .^ 2, mask; dims)
	maskedstd(A, mask; dims) = sqrt.(maskedvar(A, mask; dims))
	masked_errorofmean(A, mask; dims) = maskedstd(A, mask; dims) ./ sqrt.(sum(mask;dims))

	meandrop(A, mask::BitArray;dims) = dropdims(maskedmean(A, mask; dims); dims)
	vardrop(A, mask::BitArray;dims) = dropdims(maskedvar(A, mask; dims); dims)
	stddrop(A, mask::BitArray;dims) = dropdims(maskedstd(A, mask; dims); dims)
	errorofmeandrop(A, mask::BitArray;dims) = dropdims(masked_errorofmean(A, mask; dims); dims)

	function _prep_vals(data, center_percent=1.0; dims, func=meandrop)
		infs = findall(isinf, data.data)
		NaNs = findall(isnan, data.data)
		length(infs)+length(NaNs) == 0 || @warn "Found $(length(NaNs)) NaNs and $(length(infs)) infs for α=$(data.interaction.α), N=$(data.system_size)"
		mask = trues(size(data))
		mask[infs] .= 0
		mask[NaNs] .= 0
		return func(center_region(data.data, center_percent), mask; dims)
	end
	function bin_count(data, nbins=50, (min,max)=extrema(data); normalize=true)
		sorted_data = sort(vec(data))
		bins = range(min,max,nbins+1)
		bin_vals = [(a+b)/2 for (a,b) in zip(bins,bins[2:end])]
		indices = searchsortedfirst.(Ref(sorted_data), bins)
		normalize && (indices = indices ./ length(data))
		return bin_vals, [a-b for (a,b) in zip(indices[2:end], indices)]
	end
	function normalize_energy(vals)
		min,max = extrema(vals)
		return (vals .- min) ./ (max-min)
	end
end

# ╔═╡ 48544d55-446d-4024-8bf9-fa428de11c0b
Ns = [10,11,12,13,14,15,16]
#Ns = [7,7,7,8,9,10,11,12,13,14]

# ╔═╡ 87c93f63-927f-4348-8534-bac286bb2809
alpha6_high = EDDataDescriptor.(RandomPositionsXXZWithXField.(PositionDataDescriptor.(:noisy_chain_pbc, 1, Ns), Ref(PowerLaw(6))); prefix=HIGH_DENSITY)

# ╔═╡ b9e129d6-de8d-48f4-b07d-48046b64148f
alpha6_low = EDDataDescriptor.(RandomPositionsXXZWithXField.(PositionDataDescriptor.(:box_pbc, 1, Ns), Ref(PowerLaw(6))); prefix=LOW_DENSITY)

# ╔═╡ b428793f-ba0b-4d29-8f3d-8c28e6518334
begin
	lsr_low_descriptor = load_lsr(alpha6_low[1]).descriptor
	lsr_high_descriptor = load_lsr(alpha6_high[1]).descriptor
end#

# ╔═╡ 885244f9-c824-492b-8e41-cef5bebc32aa
ρsortperm = let ρlow = copy(lsr_low_descriptor.ρs), ρhigh = copy(lsr_high_descriptor.ρs),
	overlap = ρhigh .< maximum(ρlow)
	ρhigh[overlap] .= typemax(Float64)
	ρsortperm_full = sortperm(vcat(ρlow, ρhigh))
	ρsortperm_full[1:end-count(overlap)]
end

# ╔═╡ 899fcdd8-8cad-4c6c-a730-3ca7b558a84d
ρvals = vcat(lsr_low_descriptor.ρs, lsr_high_descriptor.ρs)[ρsortperm]

# ╔═╡ 9e513bf3-583c-4e9e-b1bc-47e01844f22a
lsr_high_descriptor.ρs

# ╔═╡ aa0c3da5-3bb2-4a3a-983d-78fd34a99ad0
lsr_low_descriptor.ρs

# ╔═╡ 9bef52d6-3de3-4566-be89-6f229749ff6a
@. W(ρ) = 1 / (ρ)

# ╔═╡ a27dbedb-3fb3-443c-b7ce-23a9aab51e83
W(lsr_high_descriptor.ρs)

# ╔═╡ 18df533c-ca08-4a27-b2af-d5d75d335d49
W(lsr_low_descriptor.ρs)

# ╔═╡ ac75eaf5-d4ef-4d2a-8447-35a9124ee0fd
@. δW(ρ, δρ) = W(ρ)*δρ/ρ

# ╔═╡ f080ed0a-f02b-4dd1-b73e-da2d498f528d
δW(transition::AbstractMatrix) = δW(transition[:,1],transition[:,2])

# ╔═╡ 0438606f-810b-4423-a1de-c416dce708f3
W_xlim = (0.46,2.08)

# ╔═╡ e144f706-23aa-4143-adcb-d8e40b558b8b
W_xticks = [0.5; 0.6:0.2:2.0]

# ╔═╡ bc11f188-7bd3-49f0-bb29-044d4f27de3c
md"""
## LSR
"""

# ╔═╡ ac92e5e6-171e-4b8e-aca3-344c0d8ea1a9
nancounterdrop(unused,mask;dims) = dropdims(mean(mask; dims);dims)

# ╔═╡ e34fa01d-466b-48e4-9a10-5f0e592909cd
lsr_nans = let lsr_means_high = hcat(_prep_vals.(load_lsr.(alpha6_high); dims=(1,2,3,4), func=nancounterdrop)...),
	lsr_means_low  = hcat(_prep_vals.(load_lsr.(alpha6_low ); dims=(1,2,3,4),func=nancounterdrop)...)
	vcat(lsr_means_low, lsr_means_high)[ρsortperm, :]
end;

# ╔═╡ 9d67a245-2d8f-48f7-81da-09b2993cda37
plot(ρvals, lsr_nans; xaxis=:log, label=Ns')

# ╔═╡ dadb5ab4-0e95-4366-84b1-ee2809313927
lsr_means = let lsr_means_high = hcat(_prep_vals.(load_lsr.(alpha6_high); dims=(1,2,3,4))...),
	lsr_means_low  = hcat(_prep_vals.(load_lsr.(alpha6_low ); dims=(1,2,3,4))...)
	vcat(lsr_means_low, lsr_means_high)[ρsortperm, :]
end;

# ╔═╡ f0665d91-1f67-4716-87ce-a3ea3b9039bf
lsr_main_W = let p = plot(;xlabel="W", ylabel=L"\langle r \rangle", legend=:topright, xaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	hline!(p, [0.5295]; label="GOE", ls=:dash, width=2, color=:limegreen)
	hline!(p, [2 * log(2)-1]; label="Poisson", ls=:dot, width=2)
	color = palette(:ice, length(Ns)+3, rev=true)#Plots.colormap("RdBu", length(Ninds)+3)[4:end]'
	plot!(p,  W(ρvals), lsr_means; palette=color, label=["N=$N" for N in Ns'])
	p
end

# ╔═╡ 0a4e0811-036d-46f5-ae89-8ffa7fb7ab1a
section_delims = [0.515,0.58,0.95]

# ╔═╡ 8b8bafa7-e325-4d36-9624-1e1efb6c972d
lsr_full_W = let p = plot(lsr_main_W),
	lfont = ("serif",15,:hcenter,:vcenter),
	ly = 0.21
	vline!(p, section_delims; label=nothing, color=:black, alpha=0.7, width=1.5)
	annotate!(p, [(0.488,ly,text("I", lfont...)),
				((section_delims[1]+section_delims[2])/2,ly,text("II", lfont...)),
				((section_delims[1]+section_delims[3])/2,ly,text("III", lfont...)),
				(1.4,ly,text("IV", lfont...))])
	p
end

# ╔═╡ 1301fc69-9d92-43c6-8274-eb3c11995b8b
savefig(lsr_full_W, joinpath(SAVEPATH, "lsr_W.pdf"))

# ╔═╡ 7f082cf1-12e8-4134-a6e3-371af39164c7
md"""
## PR
"""

# ╔═╡ 80e82ba2-8a72-4ba9-9828-42aea93ba3bd
pr_edds_high = filter(edd->iseven(edd.system_size), alpha6_high)

# ╔═╡ 4721cdab-edc9-42dc-9656-64534000ed35
pr_edds_low = filter(edd->iseven(edd.system_size), alpha6_low)

# ╔═╡ 5b30d6d4-6b12-4da3-b5f5-1e1170abaac5
pr_Ns = getproperty.(pr_edds_high, :system_size)

# ╔═╡ a35f6334-1d05-4925-8c16-07b3c1c75130
pr_means_zbasis = let pr_means_high =
	[_prep_vals(load(PRDataDescriptor(ZBasis(), edd)); dims=(1,2,3,4)) for edd in alpha6_high]
	pr_means_low =
	[_prep_vals(load(PRDataDescriptor(ZBasis(), edd)); dims=(1,2,3,4)) for edd in alpha6_low]
	vcat(hcat(pr_means_low...), hcat(pr_means_high...))[ρsortperm, :]
end;

# ╔═╡ 7e66c883-eb17-410f-b915-3d9b273af71b
pr_means_pairs = let pr_means_high =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size,PetersMethod()), edd)); dims=(1,2,3,4)) for edd in pr_edds_high]
	pr_means_low =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size,PetersMethod()), edd)); dims=(1,2,3,4)) for edd in pr_edds_low]
	vcat(hcat(pr_means_low...), hcat(pr_means_high...))[ρsortperm, :]
end;

# ╔═╡ bdd0e16f-5615-477e-acf6-670ea51b2381
pr_means_naive_shots = let pr_means_high =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size, NaivePairing(false)), edd)); dims=(1,2,4)) for edd in pr_edds_high]
	pr_means_low =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size, NaivePairing(false)), edd)); dims=(1,2,4)) for edd in pr_edds_low]
	[pr_means_low, pr_means_high]#[high_low][Ns][shots, ρ]
end;

# ╔═╡ 2eb28070-691f-4860-a7ac-87c2d9fb81d0
pr_means_naive2_shots = let pr_means_high =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size, NaivePairing(true)), edd)); dims=(1,2,4)) for edd in pr_edds_high]
	pr_means_low =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size, NaivePairing(true)), edd)); dims=(1,2,4)) for edd in pr_edds_low]
	[pr_means_low, pr_means_high]#[high_low][Ns][shots, ρ]
end;

# ╔═╡ 8a4d457c-f788-4c04-b84a-3734a4e3e865
pr_means_naive = mapreduce(Ndata->hcat(SimLib.meandrop.(Ndata;dims=1)...), vcat, pr_means_naive_shots)[ρsortperm, :];

# ╔═╡ 21752b81-f146-4ffb-8cfd-627e63266610
pr_means_naive2 = mapreduce(Ndata->hcat(SimLib.meandrop.(Ndata;dims=1)...), vcat, pr_means_naive2_shots)[ρsortperm, :];

# ╔═╡ aff5c4e6-8ed1-4b12-88b3-43e9ebdc617a
pr_means_naive_lowest =
	mapreduce(vcat, zip(pr_means_naive_shots,pr_means_naive2_shots)) do rhoblock #[[N1,N2,...],...]
		mapreduce(arrs->SimLib.meandrop(min.(arrs...);dims=1),hcat, zip(rhoblock...))
	end


# ╔═╡ 5e57efdc-5ae3-4a86-ae6d-9f0a49c4d5c7
begin
	pr_lines_in_overview = 3
	pr_linfit_index = 15
	pr_zbasis_color = palette(:ice, 5, rev=true)[2:end]
	pr_pairs_color = :Reds
	pr_naive_color = :PRGn_10
	pr_palette = [palette(:Blues,5)[3:5];
		palette(:Reds,5)[3:5];
		palette(:Greens,7)[4:6]]
end

# ╔═╡ 21352df0-3493-48ef-865e-2f6b0614b0f4
pr_main_W = let Ns = pr_Ns[end-pr_lines_in_overview+1:end],
	pr_means_zbasis = pr_means_zbasis[:, end-pr_lines_in_overview*2+2:2:end],
	pr_means_pairs = pr_means_pairs[:, end-pr_lines_in_overview+1:end],
	pr_means_naive = pr_means_naive[:,end-pr_lines_in_overview+1:end],
	palette = pr_palette,
	p = plot(;xlabel="W", ylabel=L"\mathrm{PR}/|\mathcal{H}|", legend=:topright, yaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	plot!(p, W(ρvals), pr_means_zbasis ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["z-basis N=$N" for N in Ns'])
	plot!(p, W(ρvals), pr_means_pairs ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["pair-basis N=$N" for N in Ns'], ls=:dash)
	# plot!(p, W(ρvals), pr_means_naive ./ binomial.(Ns, div.(Ns .- 1, 2))';
	# 	palette, label=["N=$N" for N in Ns'])
	p
end

# ╔═╡ c4eae1f2-95c1-4114-a1dd-3b437414c189
function pr_pairbasis(N)
	symm = symmetrized_basis(N,div(N-1,2))
	pairbasis = SimLib.PRModule.construct_basis(PairBasis(N,PetersMethod()))
	return mean(participation_ratio(symmetrize_operator(pairbasis,symm)))
end

# ╔═╡ aef7cc15-7c6a-4f17-b6aa-5d6398ba8487
PR_zbasis_vs_pairs = pr_pairbasis.(pr_Ns)

# ╔═╡ 97bf13c5-cfc2-4554-925d-759444bf3a45
md"""
## Entropy
"""

# ╔═╡ 15646ef5-cf1d-4a14-8fa9-4bc5e0842272
hce_means = let hcedds_high = HCEDataDescriptor.(div.(Ns, 2), alpha6_high),
	hcedds_low = HCEDataDescriptor.(div.(Ns, 2), alpha6_low),
	hce_means_high = hcat(_prep_vals.(load.(hcedds_high; prefix=HIGH_DENSITY); dims=(1,2,3,4,5))...),
	hce_means_low = hcat(_prep_vals.(load.(hcedds_low; prefix=LOW_DENSITY); dims=(1,2,3,4,5))...)
	vcat(hce_means_low, hce_means_high)[ρsortperm, :]
end;

# ╔═╡ d542da72-4c4a-4cfa-b490-4ae881e8538a
hce_std = let hcedds_high = HCEDataDescriptor.(div.(Ns, 2), alpha6_high),
	hcedds_low = HCEDataDescriptor.(div.(Ns, 2), alpha6_low),
	hce_means_high = hcat(_prep_vals.(load.(hcedds_high; prefix=HIGH_DENSITY); dims=(1,2,3,4,5), func=stddrop)...),
	hce_means_low = hcat(_prep_vals.(load.(hcedds_low; prefix=LOW_DENSITY); dims=(1,2,3,4,5), func=stddrop)...)
	vcat(hce_means_low, hce_means_high)[ρsortperm, :]
end;

# ╔═╡ f9491876-b0cd-4f2d-9b9f-ca24283dfb44
hce_errorofmean = let hcedds_high = HCEDataDescriptor.(div.(Ns, 2), alpha6_high),
	hcedds_low = HCEDataDescriptor.(div.(Ns, 2), alpha6_low),
	hce_means_high = hcat(_prep_vals.(load.(hcedds_high; prefix=HIGH_DENSITY); dims=(1,2,3,4,5), func=errorofmeandrop)...),
	hce_means_low = hcat(_prep_vals.(load.(hcedds_low; prefix=LOW_DENSITY); dims=(1,2,3,4,5), func=errorofmeandrop)...)
	vcat(hce_means_low, hce_means_high)[ρsortperm, :]
end;

# ╔═╡ 99d387f0-ae7e-4ce4-8ad8-6c8bf95e3819
entropy_main_W = let p = plot(;xlabel="W", ylabel="Half-chain entropy", legend=:topright, xaxis=:log,
		xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	#color = Plots.colormap("RdBu", length(Ns)+3)[4:end]'
	color = palette(:ice, length(Ns)+2, rev=true)[2:end]
	plot!(p, W(ρvals), hce_means; palette=color, label=["N=$N" for N in Ns'])
	p#
end

# ╔═╡ 5d930be7-5df7-47bf-be84-2a055106920b
pol1(x, p) = @. p[1] + x * p[2]

# ╔═╡ 241ddd4b-f075-4785-b161-2192982f2137
# ╠═╡ show_logs = false
pr_side = let index = pr_linfit_index,
	zbasis_color = pr_palette[2],
	pairbasis_color = pr_palette[5],
	p0 = [0.0,0.1],#[0.5,-0.25],
	norms = binomial.(pr_Ns, div.(pr_Ns .- 1, 2)),
	zbasis_data = pr_means_zbasis[index,:],# ./ norms,
	pairs_data = pr_means_pairs[index,:],# ./ norms,
	zbasis_fit = curve_fit(pol1, Ns, log10.(zbasis_data), p0),
	pairs_fit = curve_fit(pol1, pr_Ns, log10.(pairs_data), p0),
	p = plot(; yaxis=:log, legend=:topleft, xlabel="N",ylabel="PR")

	scatter!(p, Ns, zbasis_data; color=zbasis_color, marker=:circle, label=false)
	scatter!(p, pr_Ns, pairs_data; color=pairbasis_color, marker=:rect, label=false)
	plot!(p, Ns, 10 .^ pol1(Ns, zbasis_fit.param); color=zbasis_color, label=false)
	plot!(p, Ns, 10 .^ pol1(Ns, pairs_fit.param); color=pairbasis_color, label=false)
	plot!(p, [],[]; color=zbasis_color, marker=:circle, label="PR eigenstate vs. z-basis")
	plot!(p, [],[]; color=pairbasis_color, marker=:rect, label="PR eigenstate vs. pair-basis")
	plot!(p, Ns, 1.5 .^ (Ns./2); color=:black, ls=:dash, label=L"1.5^{N/2}")
	plot!(p, Ns, 1.5 .^ (Ns./4); color=:black, ls=:dot, label=L"1.5^{N/4}")
	plot!(p, pr_Ns, PR_zbasis_vs_pairs; color=:green, label="PR pair-basis vs. z-basis", ls=:dashdot)
	yticks!(p, 1.5 .^ (2:9), ["\$1.5^{$i}\$" for i in 2:9])
	p
end

# ╔═╡ fc0a8c34-1bc8-4003-8165-13a362257648
pr_full_W = let main = plot(pr_main_W; title="(a)", legend=:topright)
	lens!(main, [0.5,0.55], [0.14,0.32] ; inset = (1, bbox(0.36, 0, 0.3, 0.3)))
	plot!(main, [0.55,0.93], [0.32,0.5]; label=false, color=:gray, alpha=0.35)
	plot!(main, [0.55,0.93], [0.17,0.06]; label=false, color=:gray, alpha=0.35)
	vline!(main, [W(ρvals[pr_linfit_index])]; color=:black, ls=:dashdot, label=false)
	p = plot(main, plot(pr_side;title="(b)"); layout=@layout([a;b]), size=(600,800), dpi=200, titlelocation = :left)

	p
end

# ╔═╡ 81fcb0e6-ad67-4313-8af5-a6d79d1e941f
savefig(pr_full_W, joinpath(SAVEPATH, "pr_W.pdf"))

# ╔═╡ f9abd6af-f125-41ef-915d-a73b041f9690
pol1_error(x,p,perr) = @. sqrt(x^2 * perr[2]^2 + perr[1]^2)

# ╔═╡ a8b925f1-1763-43eb-9162-219b85ede984
function _cut_prediction(posdata)
	interaction = PowerLaw(6)
	N = posdata.system_size
	pbc_pairsize(x) = min(x[2]-x[1], N-x[2]+x[1])

	res = zeros(posdata.shots, length(posdata.ρs))
	for (i,ρ) in enumerate(posdata.ρs)
		geom = geometry_from_density(posdata.geometry, ρ, posdata.system_size,1)
		for j in 1:posdata.shots
			clust = SimLib.PRModule.peters_clustering!(interaction_matrix(interaction, geom, posdata.data[:,:,j,i]))
			res[j,i] = sum(pbc_pairsize, clust)
		end
	end
	return res
end

# ╔═╡ 895cace4-ae8b-4d39-8b76-3cdf09815940
entropy_per_cut(N,r) = 2* (N^2-r^2)/(4N^2-2N)

# ╔═╡ 9f8782a6-8c0b-45cb-94ec-64a07addb55b
function pair_entropy_prediction(posdata)
	N = posdata.system_size
	number_of_pair_cuts = 2*_cut_prediction(posdata)
	value_of_cut_pair = entropy_per_cut(N÷2,1)
	return number_of_pair_cuts .* value_of_cut_pair ./ N
end

# ╔═╡ f73cc7af-6f60-48c2-bdf1-3de67260e9b6
entropy_plot_full_W = let detailsat = reverse([9,40,45]),
	details_color = palette(:Set2_3),
	entropy_main = plot(entropy_main_W; legend=:bottomleft),
	prediction_posdata = load_positions(:box_pbc,1,16;prefix=LOW_DENSITY)
	# ,detail_plots = [entropy_subplot(1:7,at) for at in detailsat]

	for (at, color) in zip(detailsat, details_color)
		vline!(entropy_main, [W(ρvals[at])]; color, alpha = 0.9, ls=:dash, label=nothing)
	end

	fits = [curve_fit(pol1, Ns, hce_means[at,:], [1.0,0.5]) for at in detailsat]

	fitdata = reduce(hcat, pol1(Ns,fit.param) for fit in fits)
	fiterrors = reduce(hcat, pol1_error(Ns,fit.param,stderror(fit)) for fit in fits)

	plot!(entropy_main, Ns, fitdata; ribbon=fiterrors,
		xlabel="N", label=false, palette=details_color,
		inset=(1,bbox(0.38,0.05,0.55,0.6)), subplot=2)
	scatter!(entropy_main, Ns, hce_means[detailsat,:]'; marker=:x, color=:black,
		subplot=2, label=false, legend=:right)

	# inset legend
	scatter!(entropy_main, [], []; marker=:x, color=:black, subplot=2, label="data", legend=:topleft)
	plot!(entropy_main,[],[];ribbon=[], label="fit", subplot=2, color=:gray)
	fit_to_annotation(p,color) = text("S=$(round(p[1];sigdigits=2)) + $(round(p[2];sigdigits=2))N"; pointsize=8, color,halign=:left)
	texts = [fit_to_annotation(fit.param, details_color[i]) for (i,fit) in enumerate(fits)]
	annotate!([(14,4.6,texts[1]), (13,2.2,texts[2]),(11,1.15,texts[3])];subplot=2, palette=details_color)

	### Prediction
	plot!(W(prediction_posdata.ρs)[1:23],
		vec(mean(pair_entropy_prediction(prediction_posdata);dims=1))[1:23];
		ls = :dot, width=2, color=:black, label="prediction")
end

# ╔═╡ 9f63456f-9da1-4203-ad01-8d159ebd1bee
savefig(entropy_plot_full_W,joinpath(SAVEPATH, "hce_W.pdf"))

# ╔═╡ 300ae731-3551-43b8-ab40-924d521fdb7b
md"""
### Entropy Variance
"""

# ╔═╡ 941cc59d-c4f3-45bf-885b-290a2a7634f2
entropy_var_main_W = let p = plot(;xlabel="W", ylabel="Half-chain entropy std.dev", legend=:topleft, xaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	color = palette(:ice, length(Ns)+2, rev=true)[2:end]
	plot!(p, W(ρvals), hce_std; palette=color, label=["N=$N" for N in Ns'])
	p
end

# ╔═╡ 67bf03ac-bc0c-4fce-afc9-edd8938925ee
entropy_var_full_parabola_W = let fitrange = 38:42,
	fitplotrange = fitrange,
	p = plot(entropy_var_main_W; legend=:outside, legendposition=(1.0,0.0), dpi=200)
	@. model(x, p) = p[1] + p[2]*(x-p[3])^2

	fits = [curve_fit(model, W(ρvals[fitrange]), col[fitrange], [1.0, -30, 0.7]) for col in eachcol(hce_std)]
	for fit in fits
		@show fit.param
	end
	fitplots = [model(W(ρvals[fitplotrange]), fit.param) for fit in fits]
	maxima_loc = [fit.param[3] for fit in fits]
	maxima_val = [fit.param[1] for fit in fits]
	maxima_std = [stderror(fit)[3] for fit in fits]

	plot!(p, W(ρvals[fitplotrange]), fitplots; color=:black, ls=:dash, alpha=0.9, label=false)
	plot!(p, [], []; color=:black, ls=:dash, alpha=0.9, label="Quad. fit")
		#label=reshape(["$(round(fit.param[3];sigdigits=4)) ± $(round(stderror(fit)[3];sigdigits=1))" for fit in fits], 1, :))
	scatter!(p, maxima_loc, maxima_val; label="Maxima", marker=:x, color=:black)

	plot!(p, Ns, maxima_loc; yerror=maxima_std,
	 	legend=nothing, xlabel="N",ylabel=L"W_{max}", inset=(1,bbox(0.65,0.05,0.5,0.35)), subplot=2)
	hline!(p, [1/(2*0.7475979202)];subplot=2,color=:black, ls=:dot,width=1,label="Reyni limit"	)
	annotate!(p, [(15.2, 0.666, ("Rényi limit",8))]; subplot=2)
	# inset = plot!(p, Ns, maxima_loc; yerror=maxima_std,
	# 	inset=(1,bbox(0.14,0.05,0.4,0.3)), subplot=2, legend=nothing, xlabel="N",ylabel="Maximum")

	p
end

# ╔═╡ c60ded9a-a363-4f99-96a9-9f073065a8c4
savefig(entropy_var_full_parabola_W, joinpath(SAVEPATH, "entropy_variance_W.pdf"))

# ╔═╡ c2f1cdc5-a6e8-4fbf-909f-bfb966709ddb
md"""
## Thouless Parameter
"""

# ╔═╡ 81f27a5a-9765-4f16-9174-cd7e217e7a41
tp_means_sz = let pr_means_high = hcat(_prep_vals.(load_el.("sz", alpha6_high); dims=(1,2,3,4))...),
	pr_means_low  = hcat(_prep_vals.(load_el.("sz", alpha6_low ); dims=(1,2,3,4))...)
	vcat(pr_means_low, pr_means_high)[ρsortperm, :]
end;

# ╔═╡ 20721250-f8e6-4c55-8b8d-d2d20db10602
tp_means_szsz = let pr_means_high = hcat(_prep_vals.(load_el.("szsz", alpha6_high); dims=(1,2,3,4))...),
	pr_means_low  = hcat(_prep_vals.(load_el.("szsz", alpha6_low ); dims=(1,2,3,4))...)
	vcat(pr_means_low, pr_means_high)[ρsortperm, :]
end;

# ╔═╡ e30bae01-4922-4dc6-82ce-28f5d8b4e368
tp_means_hop = let pr_means_high = hcat(_prep_vals.(load_el.("hopping", alpha6_high); dims=(1,2,3,4))...),
	pr_means_low  = hcat(_prep_vals.(load_el.("hopping", alpha6_low ); dims=(1,2,3,4))...)
	vcat(pr_means_low, pr_means_high)[ρsortperm, :]
end;

# ╔═╡ b5d37d9f-7c44-4218-9e0b-1cb6e4915ad3
interpolate(x1,x2, steps=201) = x1 == x2 ? [x1] : collect(range(x1,x2; length=steps))

# ╔═╡ 1376a743-1355-49a5-977d-264e63bc0881
function tp_min_variance(data,ρs; iterations=5)
	for i in 1:iterations
		means_spread = SimLib.stddrop(data, dims=2)
		minind = argmin(means_spread)
		i == iterations && return ρs[minind]
		left = max(1,minind-1)
		right = min(minind+1,length(ρs))
		data = hcat(interpolate(data[left,:], data[minind,:])...,
			interpolate(data[minind, :], data[right, :])...)'
		ρs = vcat(interpolate(ρs[left,:], ρs[minind,:])...,
			interpolate(ρs[minind, :], ρs[right, :])...)
	end
end

# ╔═╡ f69d1bb5-20cf-416f-a581-4acb56676783
function tp_find_transition(data, ρs)
	transition = tp_min_variance(data,ρs)
	two_line_intersections = [tp_min_variance(data[:,[i,j]],ρs) for i in 1:size(data,2) for j in i+1:size(data,2)]
	return [transition, std(two_line_intersections)]
end

# ╔═╡ 5d1d0e06-5269-4b65-ba36-c71ef6afaa52
begin
	tp_transition_window = 20:length(ρvals)-3
	tp_sz_transition = tp_find_transition(tp_means_sz[tp_transition_window, :], ρvals[tp_transition_window])
	tp_szsz_transition = tp_find_transition(tp_means_szsz[tp_transition_window, :], ρvals[tp_transition_window])
	tp_hop_transition = tp_find_transition(tp_means_hop[tp_transition_window, :], ρvals[tp_transition_window])
end

# ╔═╡ 99dbe4d7-2919-4677-8f97-0566ff66b146
tp_N_window = 1

# ╔═╡ 01cb22ef-1d44-4910-8e0c-77da2411081d
begin
	tp_sz_transitions = hcat([tp_find_transition(tp_means_sz[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
	tp_szsz_transitions = hcat([tp_find_transition(tp_means_szsz[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
	tp_hop_transitions = hcat([tp_find_transition(tp_means_hop[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
end

# ╔═╡ 6a05e029-0b39-4ea5-a7ef-e52defbe771e
tp_nans = let pr_means_high = hcat(_prep_vals.(load_el.("hopping", alpha6_high); dims=(1,2,3,4), func=nancounterdrop)...),
	pr_means_low  = hcat(_prep_vals.(load_el.("hopping", alpha6_low ); dims=(1,2,3,4), func=nancounterdrop)...)
	vcat(pr_means_low, pr_means_high)[ρsortperm, :]
end;

# ╔═╡ 1b3e8683-02b7-45a1-a65e-cd74083bb2c5
tp_sz_transition

# ╔═╡ 30100b77-d704-4f7e-9aea-9a2953b86c3b
thouless_main_W = let p = plot(;xlabel="W", ylabel=L"\langle \mathcal{G} \rangle", legend=:topright, xaxis=log, xlim=W_xlim)
	xticks!(p, W_xticks, string.(W_xticks))
	color = palette(:ice, length(Ns)+3; rev=true)[2:end]
	plot!(p, W(ρvals), tp_means_sz; palette=color, label=["N=$N" for N in Ns'])
	p
end

# ╔═╡ 5a366799-c8a1-4fe5-afe2-909524e30e5e
thouless_full_W = let p = plot(thouless_main_W; dpi=200, legend=:topright, legend_columns=-1)
	loc, width = tp_sz_transition
	#vline!([W(tp_sz_transition)]; label="Transition", color=:darkred, alpha=0.7)
	vspan!(W([loc-width,loc+width]); label="Crossover", color=:darkred, alpha=0.7)

	# y = W(tp_sz_transitions[:,1])
	# δy = δW(tp_sz_transitions)
	# plot!(p, Ns[(1+tp_N_window):(end-tp_N_window)], y; subplot=2, inset=(1, bbox(0.55,0.05,0.4,0.35)), legend=false, xlabel="System size N", ylabel="Crossover point", yerror=δy)
	# lens!(p, [0.65,0.75],[-1.9,-0.3] ; inset = (1, bbox(0.05, 0.7, 0.25, 0.25)))
	p
end

# ╔═╡ c3bb699c-e941-4a46-b430-bc1222b2ed09
savefig(thouless_full_W, joinpath(SAVEPATH, "thouless_W.pdf"))

# ╔═╡ da20058c-4073-4ac6-badd-00adfa0860eb
hilberspacesize(N) = binomial(N,div(N-1,2))

# ╔═╡ f37acd5b-b8a4-4caf-91a0-9cba909e5a78
md"""
## Entropy per cut
"""

# ╔═╡ 5e099ce4-8f01-4087-bdb6-ba1d6ab23165
let p = plot(;xlabel="Number of pairs N",ylabel="Entropy per cut", dpi=200)
	plot!(p, 3:20, entropy_per_cut.(3:20,0); label="r=0")
	plot!(p, 4:20, entropy_per_cut.(4:20,1); label="r=1")
	plot!(p, 5:20, entropy_per_cut.(5:20,2); label="r=2")
	plot!(p, 6:20, entropy_per_cut.(6:20,3); label="r=3")
	plot!(p, 7:20, entropy_per_cut.(7:20,4); label="r=4")
	hline!([0.5];label=false, color=:black, width=2,ls=:dot)
	savefig(p, joinpath(SAVEPATH, "entropy_per_cut.pdf"))
	p
end

# ╔═╡ Cell order:
# ╠═d05f72ed-31ac-4983-b510-ce6d6527879a
# ╠═a35a6fb5-dd1e-4998-8c2a-097acb9bf7d4
# ╠═1516ba16-a8ef-11ec-3979-1f7edd32dc6d
# ╠═fc5ab3d2-acdd-47a7-bbdc-c40b28976907
# ╠═f92f65b3-59f0-4067-a744-3105c9411d0a
# ╠═48544d55-446d-4024-8bf9-fa428de11c0b
# ╠═87c93f63-927f-4348-8534-bac286bb2809
# ╠═b9e129d6-de8d-48f4-b07d-48046b64148f
# ╠═b428793f-ba0b-4d29-8f3d-8c28e6518334
# ╠═885244f9-c824-492b-8e41-cef5bebc32aa
# ╠═a27dbedb-3fb3-443c-b7ce-23a9aab51e83
# ╠═18df533c-ca08-4a27-b2af-d5d75d335d49
# ╠═899fcdd8-8cad-4c6c-a730-3ca7b558a84d
# ╠═9e513bf3-583c-4e9e-b1bc-47e01844f22a
# ╠═aa0c3da5-3bb2-4a3a-983d-78fd34a99ad0
# ╠═9bef52d6-3de3-4566-be89-6f229749ff6a
# ╠═ac75eaf5-d4ef-4d2a-8447-35a9124ee0fd
# ╠═f080ed0a-f02b-4dd1-b73e-da2d498f528d
# ╠═0438606f-810b-4423-a1de-c416dce708f3
# ╠═e144f706-23aa-4143-adcb-d8e40b558b8b
# ╟─bc11f188-7bd3-49f0-bb29-044d4f27de3c
# ╠═ac92e5e6-171e-4b8e-aca3-344c0d8ea1a9
# ╠═e34fa01d-466b-48e4-9a10-5f0e592909cd
# ╠═9d67a245-2d8f-48f7-81da-09b2993cda37
# ╠═dadb5ab4-0e95-4366-84b1-ee2809313927
# ╠═f0665d91-1f67-4716-87ce-a3ea3b9039bf
# ╠═0a4e0811-036d-46f5-ae89-8ffa7fb7ab1a
# ╠═8b8bafa7-e325-4d36-9624-1e1efb6c972d
# ╠═1301fc69-9d92-43c6-8274-eb3c11995b8b
# ╟─7f082cf1-12e8-4134-a6e3-371af39164c7
# ╠═80e82ba2-8a72-4ba9-9828-42aea93ba3bd
# ╠═4721cdab-edc9-42dc-9656-64534000ed35
# ╠═5b30d6d4-6b12-4da3-b5f5-1e1170abaac5
# ╠═a35f6334-1d05-4925-8c16-07b3c1c75130
# ╠═7e66c883-eb17-410f-b915-3d9b273af71b
# ╠═bdd0e16f-5615-477e-acf6-670ea51b2381
# ╠═2eb28070-691f-4860-a7ac-87c2d9fb81d0
# ╠═8a4d457c-f788-4c04-b84a-3734a4e3e865
# ╠═21752b81-f146-4ffb-8cfd-627e63266610
# ╠═aff5c4e6-8ed1-4b12-88b3-43e9ebdc617a
# ╠═5e57efdc-5ae3-4a86-ae6d-9f0a49c4d5c7
# ╠═21352df0-3493-48ef-865e-2f6b0614b0f4
# ╠═241ddd4b-f075-4785-b161-2192982f2137
# ╠═fc0a8c34-1bc8-4003-8165-13a362257648
# ╠═c4eae1f2-95c1-4114-a1dd-3b437414c189
# ╠═aef7cc15-7c6a-4f17-b6aa-5d6398ba8487
# ╠═81fcb0e6-ad67-4313-8af5-a6d79d1e941f
# ╟─97bf13c5-cfc2-4554-925d-759444bf3a45
# ╠═15646ef5-cf1d-4a14-8fa9-4bc5e0842272
# ╠═d542da72-4c4a-4cfa-b490-4ae881e8538a
# ╠═f9491876-b0cd-4f2d-9b9f-ca24283dfb44
# ╠═99d387f0-ae7e-4ce4-8ad8-6c8bf95e3819
# ╠═5d930be7-5df7-47bf-be84-2a055106920b
# ╠═f9abd6af-f125-41ef-915d-a73b041f9690
# ╠═f73cc7af-6f60-48c2-bdf1-3de67260e9b6
# ╠═9f8782a6-8c0b-45cb-94ec-64a07addb55b
# ╠═a8b925f1-1763-43eb-9162-219b85ede984
# ╠═895cace4-ae8b-4d39-8b76-3cdf09815940
# ╠═9f63456f-9da1-4203-ad01-8d159ebd1bee
# ╟─300ae731-3551-43b8-ab40-924d521fdb7b
# ╠═941cc59d-c4f3-45bf-885b-290a2a7634f2
# ╠═67bf03ac-bc0c-4fce-afc9-edd8938925ee
# ╠═c60ded9a-a363-4f99-96a9-9f073065a8c4
# ╟─c2f1cdc5-a6e8-4fbf-909f-bfb966709ddb
# ╠═81f27a5a-9765-4f16-9174-cd7e217e7a41
# ╠═20721250-f8e6-4c55-8b8d-d2d20db10602
# ╠═e30bae01-4922-4dc6-82ce-28f5d8b4e368
# ╠═b5d37d9f-7c44-4218-9e0b-1cb6e4915ad3
# ╠═f69d1bb5-20cf-416f-a581-4acb56676783
# ╠═1376a743-1355-49a5-977d-264e63bc0881
# ╠═5d1d0e06-5269-4b65-ba36-c71ef6afaa52
# ╠═99dbe4d7-2919-4677-8f97-0566ff66b146
# ╠═01cb22ef-1d44-4910-8e0c-77da2411081d
# ╠═6a05e029-0b39-4ea5-a7ef-e52defbe771e
# ╠═1b3e8683-02b7-45a1-a65e-cd74083bb2c5
# ╠═30100b77-d704-4f7e-9aea-9a2953b86c3b
# ╠═5a366799-c8a1-4fe5-afe2-909524e30e5e
# ╠═c3bb699c-e941-4a46-b430-bc1222b2ed09
# ╠═da20058c-4073-4ac6-badd-00adfa0860eb
# ╟─f37acd5b-b8a4-4caf-91a0-9cba909e5a78
# ╠═5e099ce4-8f01-4087-bdb6-ba1d6ab23165
