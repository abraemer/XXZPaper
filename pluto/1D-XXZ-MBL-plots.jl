### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 2a544bad-35cb-4229-b75e-8d3283b8c288
md"""
# Setup and base definitions
"""

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
	Pkg.add(["Plots", "PyPlot", "Statistics", "LinearAlgebra","LaTeXStrings","LsqFit","SpinSymmetry","PlutoUI","Bootstrap"])
	Pkg.add(url="https://github.com/pablosanjose/FilteredMatrices.jl")
	Pkg.add(path=expanduser("~/julia/XXZNumerics"))
	Pkg.add(path=expanduser("~/julia/SimLib"))
	using Plots,Statistics,LinearAlgebra,SimLib,LaTeXStrings,LsqFit,XXZNumerics,SpinSymmetry,PlutoUI,Bootstrap
	pyplot()
end

# ╔═╡ 4d171def-b37a-44ea-a192-03a88debe63c
TableOfContents(indent=true)

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

# ╔═╡ ed2305f5-7bc9-43a8-b5eb-75b5a8c00eca
begin
	function nan_mask(data)
		mask = trues(size(data))
		mask[isnan.(data)] .= 0
		return mask
	end
	
	_ndim(arr) = length(size(arr))
	_default_dims(arr) = Tuple(1:_ndim(arr))
	mean_ignore_nan(data; dims=_default_dims(data)) = maskedmean(data, nan_mask(data); dims)
	std_ignore_nan(data; dims=_default_dims(data)) = maskedstd(data, nan_mask(data); dims)
	function error_of_mean_ignore_nan(data; dims=_default_dims(data))
		mask = nan_mask(data)
		return maskedstd(data, mask; dims) ./ sqrt.(sum(mask; dims))
	end
end

# ╔═╡ da20058c-4073-4ac6-badd-00adfa0860eb
hilberspacesize(N) = binomial(N,div(N-1,2))

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
# LSR
"""

# ╔═╡ 799f75b3-6bbd-4325-89b2-6a24f780f51a
md"""
## Load data
"""

# ╔═╡ ac92e5e6-171e-4b8e-aca3-344c0d8ea1a9
nancounterdrop(unused,mask;dims) = dropdims(mean(mask; dims);dims)

# ╔═╡ 33b9ecab-6984-416a-8127-6367da4d04ef
lsr_means, lsr_errorofmean = let lsr_shot_data_high = _prep_vals.(load_lsr.(alpha6_high); dims=(1,2,4)), #[N][shot,rho]
	lsr_shot_data_low =
	_prep_vals.(load_lsr.(alpha6_low ); dims=(1,2,4)) #[N][shot,rho]
	lsr_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(lsr_shot_data_low,lsr_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean_ignore_nan.(lsr_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(error_of_mean_ignore_nan.(lsr_shot_data; dims=1))...)
	means, errorofmean
end;

# ╔═╡ 336c5f98-71a9-40e8-9821-7c56033a40f9
md"""
## Make and save plots
"""

# ╔═╡ f0665d91-1f67-4716-87ce-a3ea3b9039bf
lsr_main_W = let p = plot(;xlabel="W", ylabel=L"\langle r \rangle", legend=:topright, xaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	hline!(p, [0.5295]; label="GOE", ls=:dash, width=2, color=:limegreen)
	hline!(p, [2 * log(2)-1]; label="Poisson", ls=:dot, width=2)
	color = palette(:ice, length(Ns)+3, rev=true)#Plots.colormap("RdBu", length(Ninds)+3)[4:end]'
	plot!(p,  W(ρvals), lsr_means; palette=color, label=["N=$N" for N in Ns'])
	p
end

# ╔═╡ a5c221a6-c6c2-4301-b93b-7178859bcf68
## error of mean as ribbon
## not visible
let p = plot(;xlabel="W", ylabel=L"\langle r \rangle", legend=:topright, xaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	hline!(p, [0.5295]; label="GOE", ls=:dash, width=2, color=:limegreen)
	hline!(p, [2 * log(2)-1]; label="Poisson", ls=:dot, width=2)
	color = palette(:ice, length(Ns)+3, rev=true)#Plots.colormap("RdBu", length(Ninds)+3)[4:end]'
	plot!(p,  W(ρvals), lsr_means; palette=color, label=["N=$N" for N in Ns'], ribbon=lsr_errorofmean)
	p
end

# ╔═╡ 0a4e0811-036d-46f5-ae89-8ffa7fb7ab1a
section_delims = [0.515,0.58,0.95]

# ╔═╡ 8b8bafa7-e325-4d36-9624-1e1efb6c972d
lsr_full_W = let p = plot(lsr_main_W; dpi=200,tickfontsize=10, labelfontsize=12),
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
# PR
"""

# ╔═╡ 81071ead-3233-42bc-a7cb-c0aa4fc33ffe
md"""
## Load data
"""

# ╔═╡ 80e82ba2-8a72-4ba9-9828-42aea93ba3bd
pr_edds_high = filter(edd->iseven(edd.system_size), alpha6_high)

# ╔═╡ 4721cdab-edc9-42dc-9656-64534000ed35
pr_edds_low = filter(edd->iseven(edd.system_size), alpha6_low)

# ╔═╡ 5b30d6d4-6b12-4da3-b5f5-1e1170abaac5
pr_Ns = getproperty.(pr_edds_high, :system_size)

# ╔═╡ bbea8da3-c3a5-440a-875d-99ea1293e546
pr_means_zbasis, pr_errorofmean_zbasis = let pr_shot_data_high = [_prep_vals(load(PRDataDescriptor(ZBasis(), edd)); dims=(1,2,4)) for edd in alpha6_high], #[N][shot,rho]
	pr_shot_data_low =
	[_prep_vals(load(PRDataDescriptor(ZBasis(), edd)); dims=(1,2,4)) for edd in alpha6_low] #[N][shot,rho]
	pr_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(pr_shot_data_low,pr_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean_ignore_nan.(pr_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(error_of_mean_ignore_nan.(pr_shot_data; dims=1))...)
	means, errorofmean
end;

# ╔═╡ b1e54703-202c-47a1-8a15-c8f4e86e9387
pr_means_pairs, pr_errorofmean_pairs = let pr_shot_data_high = [_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size,PetersMethod()), edd)); dims=(1,2,4)) for edd in pr_edds_high], #[N][shot,rho]
	pr_shot_data_low =
	[_prep_vals(load(PRDataDescriptor(PairBasis(edd.system_size,PetersMethod()), edd)); dims=(1,2,4)) for edd in pr_edds_low] #[N][shot,rho]
	pr_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(pr_shot_data_low,pr_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean.(pr_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(std.(pr_shot_data; dims=1)) ./ sqrt.(size.(pr_shot_data, 1))...)
	means, errorofmean
end;

# ╔═╡ c4eae1f2-95c1-4114-a1dd-3b437414c189
function pr_pairbasis(N)
	symm = symmetrized_basis(N,div(N-1,2))
	pairbasis = SimLib.PRModule.construct_basis(PairBasis(N,PetersMethod()))
	return mean(participation_ratio(symmetrize_operator(pairbasis,symm)))
end

# ╔═╡ aef7cc15-7c6a-4f17-b6aa-5d6398ba8487
PR_zbasis_vs_pairs = pr_pairbasis.(pr_Ns)

# ╔═╡ 9db9970d-436b-46d5-bc9b-65477dc529f1
md"""
## Make and save plots
"""

# ╔═╡ 5e57efdc-5ae3-4a86-ae6d-9f0a49c4d5c7
begin
	pr_lines_in_overview = 3
	pr_linfit_index = 9
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
	palette = pr_palette,
	p = plot(;xlabel="W", ylabel=L"\mathrm{PR}/|\mathcal{H}|", legend=:topright, yaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	plot!(p, W(ρvals), pr_means_zbasis ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["z-basis N=$N" for N in Ns'])
	plot!(p, W(ρvals), pr_means_pairs ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["pair-basis N=$N" for N in Ns'], ls=:dash)
	p
end

# ╔═╡ aad4eeb4-8fdb-40e3-9301-49bd9a4de439
## ribbon indicates error of mean
## -> line width is of same order
let Ns = pr_Ns[end-pr_lines_in_overview+1:end],
	pr_means_zbasis = pr_means_zbasis[:, end-pr_lines_in_overview*2+2:2:end],
	pr_eom_zbasis = pr_errorofmean_zbasis[:, end-pr_lines_in_overview*2+2:2:end],
	pr_means_pairs = pr_means_pairs[:, end-pr_lines_in_overview+1:end],
	pr_eom_pairs = pr_errorofmean_pairs[:, end-pr_lines_in_overview+1:end],
	palette = pr_palette,
	p = plot(;xlabel="W", ylabel=L"\mathrm{PR}/|\mathcal{H}|", legend=:topright, yaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	plot!(p, W(ρvals), pr_means_zbasis ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["z-basis N=$N" for N in Ns'],
		ribbon = pr_eom_zbasis ./ binomial.(Ns, div.(Ns .- 1, 2))')
	plot!(p, W(ρvals), pr_means_pairs ./ binomial.(Ns, div.(Ns .- 1, 2))';
		palette, label=["pair-basis N=$N" for N in Ns'], ls=:dash,
		ribbon = pr_eom_pairs ./ binomial.(Ns, div.(Ns .- 1, 2))')
	p
end

# ╔═╡ 97bf13c5-cfc2-4554-925d-759444bf3a45
md"""
# Entropy
"""

# ╔═╡ 0f802771-29b5-472d-bcee-c6987bd97aea
md"""
## Load data
"""

# ╔═╡ eb5b02be-5195-47b0-a300-f4071e0f56c3
hce_shot_data = let hcedds_high = HCEDataDescriptor.(div.(Ns, 2), alpha6_high),
	hcedds_low = HCEDataDescriptor.(div.(Ns, 2), alpha6_low),
	hce_shots_high = _prep_vals.(load.(hcedds_high; prefix=HIGH_DENSITY); dims=(1,2,3,5)),
	hce_shots_low = _prep_vals.(load.(hcedds_low; prefix=LOW_DENSITY); dims=(1,2,3,5)) #[N][shot,rho_low]
	[hcat(low, high)[:,ρsortperm] for (low,high) in zip(hce_shots_low,hce_shots_high)]
end; # [N][shot,rho]

# ╔═╡ b3ca4a97-d758-4fa5-b881-fc34063b5a28
hce_means = hcat(vec.(mean_ignore_nan.(hce_shot_data;dims=1))...); # [rho,N]

# ╔═╡ 6ea0c192-202f-472a-a7f0-51a1a8e72256
hce_std, hce_errorofmean, hce_std_std = let bs_std=zeros(size(hce_shot_data[1],2), length(hce_shot_data)), bs_eom = zeros(size(bs_std)), bs_std_std = zeros(size(bs_std))
	for (i, data) in enumerate(hce_shot_data)
		for (j, col) in enumerate(eachcol(data))
			mask = nan_mask(col)
			bs = bootstrap(std, col[mask], BasicSampling(10000))
			bs_std[j,i] = original(bs)[1]
			bs_eom[j,i] = original(bs)[1] / sqrt(sum(mask))
			bs_std_std[j,i] = stderror(bs)[1]
		end
		@info "$i done"
	end
	bs_std, bs_eom, bs_std_std
end; # [rho,N], [rho,N]

# ╔═╡ 05ccfc6e-92c2-45b5-b49d-ead83bf1ad4c
md"""
## Pair prediction code
"""

# ╔═╡ 895cace4-ae8b-4d39-8b76-3cdf09815940
entropy_per_cut(N,r) = 2* (N^2-r^2)/(4N^2-2N)

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

# ╔═╡ 9f8782a6-8c0b-45cb-94ec-64a07addb55b
function pair_entropy_prediction(posdata)
	N = posdata.system_size
	number_of_pair_cuts = 2*_cut_prediction(posdata)
	value_of_cut_pair = entropy_per_cut(N÷2,1)
	return number_of_pair_cuts .* value_of_cut_pair ./ N
end

# ╔═╡ 0504d150-3841-4a44-9d06-f5d27b3f9b39
md"""
## Make and save plots
"""

# ╔═╡ 99d387f0-ae7e-4ce4-8ad8-6c8bf95e3819
entropy_main_W = let p = plot(;xlabel="W", ylabel="Half-chain entropy", legend=:topright, xaxis=:log,
		xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	#color = Plots.colormap("RdBu", length(Ns)+3)[4:end]'
	color = palette(:ice, length(Ns)+2, rev=true)[2:end]
	plot!(p, W(ρvals), hce_means; palette=color, label=["N=$N" for N in Ns'])
	p#
end

# ╔═╡ bbf4bd03-23eb-4dd9-b91c-d92eb3c4ff6a
## include error of mean as ribbon
## -> invisible
entropy_main_W_errors = let p = plot(;xlabel="W", ylabel="Half-chain entropy", legend=:topright, xaxis=:log,
		xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	#color = Plots.colormap("RdBu", length(Ns)+3)[4:end]'
	color = palette(:ice, length(Ns)+2, rev=true)[2:end]
	plot!(p, W(ρvals), hce_means; palette=color, label=["N=$N" for N in Ns'], 
		ribbon=hce_errorofmean, fillalpha=0.2) 
	# only consider statistical error of disorder -> _prep_vals also divided by #states so multiply back on
	# this only happend for entropy because for the other indicators average is over eigenstate
	# but for entropy it's actually over the cut location
	p
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
	zbasis_data = pr_means_zbasis[index,:],
	pairs_data = pr_means_pairs[index,:],
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
	p = plot(main, plot(pr_side;title="(b)"); layout=@layout([a;b]), size=(600,800), dpi=200, titlelocation = :left, tickfontsize=10, labelfontsize=12)

	p
end

# ╔═╡ 81fcb0e6-ad67-4313-8af5-a6d79d1e941f
savefig(pr_full_W, joinpath(SAVEPATH, "pr_W.pdf"))

# ╔═╡ f9abd6af-f125-41ef-915d-a73b041f9690
pol1_error(x,p,perr) = @. sqrt(x^2 * perr[2]^2 + perr[1]^2)

# ╔═╡ f73cc7af-6f60-48c2-bdf1-3de67260e9b6
entropy_plot_full_W = let detailsat = reverse([9,40,45]),
	details_color = palette(:Set2_3),
	entropy_main = plot(entropy_main_W; legend=:bottomleft, dpi=200),
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
	annotate!(entropy_main, [(14,4.6,texts[1]), (13,2.2,texts[2]),(11,1.15,texts[3])];subplot=2, palette=details_color)
			
	### Prediction
	plot!(entropy_main, W(prediction_posdata.ρs)[1:23],
		vec(mean(pair_entropy_prediction(prediction_posdata);dims=1))[1:23];
		ls = :dot, width=2, color=:black, label="prediction")
end

# ╔═╡ 6b5ad674-b541-45ec-8bbf-76543f77297f
# The final plot was assembled in an image editor to arrange subplots correctly.

# ╔═╡ bfc60db7-04c1-4934-a3de-c6b18d841a2f
entropy_plot_full_W_new = let detailsat = reverse([9,40,45]),
	details_color = palette(:Set2_3),
	entropy_main = plot(entropy_main_W_errors; legend=false, labelfontsize=12, tickfontsize=10),
	prediction_posdata = load_positions(:box_pbc,1,16;prefix=LOW_DENSITY)
	# ,detail_plots = [entropy_subplot(1:7,at) for at in detailsat]

	for (at, color) in zip(detailsat, details_color)
		vline!(entropy_main, [W(ρvals[at])]; color, alpha = 0.9, ls=:dash, label=nothing)
	end
			
	### Prediction
	preddata_full = pair_entropy_prediction(prediction_posdata)
	preddata_mean = vec(mean(preddata_full;dims=1))[1:end]
	preddata_std = vec(std(preddata_full;dims=1))[1:end] ./ sqrt(size(preddata_full,1))
	plot!(entropy_main, W(prediction_posdata.ρs)[1:end], preddata_mean;
		ls = :dash, width=1, color=:red, label="prediction", ribbon=preddata_std, fillalpha=0.2)
	# lens
	p_lens = plot(entropy_main)
	plot!(p_lens; legend=:outerleft, ylabel="Half-chain entropy", xlim=[1.05,2.75], ylim=[0.7,1.2], ticks=:native, tickfontsize=10, labelfontsize=12)
	#p_lens = plot(;xlabel="W", legend=false, xaxis=:log, xlim=[1.5,2.0], ylim=[0.7,1.0])
	xticks!(p_lens, 1.1:0.2:2.7, string.(1.1:0.2:2.7))
	# plot!(p_lens, W(ρvals), hce_means; palette=palette(:ice, length(Ns)+2, rev=true)[2:end])
	# plot!(p_lens, W(prediction_posdata.ρs)[1:23],
	# 	vec(mean(pair_entropy_prediction(prediction_posdata);dims=1))[1:23];
	# 	ls = :dot, width=2, color=:black, label="prediction")

	# curve fits
	fits = [curve_fit(pol1, Ns, hce_means[at,:], [1.0,0.5]) for at in detailsat]

	fitdata = reduce(hcat, pol1(Ns,fit.param) for fit in fits)
	fiterrors = reduce(hcat, pol1_error(Ns,fit.param,stderror(fit)) for fit in fits)

	p_linfit = entropy_main
	
	plot!(p_linfit, Ns, fitdata; ribbon=fiterrors,
		xlabel="N", label=false, palette=details_color,
		inset=(1,bbox(0.38,0.05,0.55,0.6)), subplot=2)
	scatter!(p_linfit, Ns, hce_means[detailsat,:]'; marker=:x, color=:black,label=false, legend=:right, subplot=2)

	# inset legend
	scatter!(p_linfit, [], []; marker=:x, color=:black, label="data", legend=:topleft, subplot=2)
	plot!(p_linfit,[],[];ribbon=[], label="fit", color=:gray, subplot=2)
	fit_to_annotation(p,color) = text("S=$(round(p[1];sigdigits=2)) + $(round(p[2];sigdigits=2))N"; pointsize=8, color,halign=:left)
	texts = [fit_to_annotation(fit.param, details_color[i]) for (i,fit) in enumerate(fits)]
	annotate!(p_linfit, [(14,4.6,texts[1]), (13,2.2,texts[2]),(11,1.15,texts[3])]; palette=details_color, subplot=2)

	# draw lens effect
	plot!(entropy_main, [2.5,1.05,1.05,2.5],[0.7,0.7,1.2,1.2]; color=:black, alpha=0.35, label=false)
	#plot!(entropy_main, [0.0, 1.0],[0.0,1.0]; inset=(2,bbox(0.6,0.9,0.4,0.3)), label=false, color=:gray, axis=false, ticks=false, background=false, background_color=:transparent, subplot=3)	

	#plot!(entropy_main, p_linfit; inset=(1,bbox(0.38,0.05,0.55,0.6)), subplot=2)
	layout = grid(2,1; heights=(0.666,0.334))#[1;1]#@layout [a{0.5h};b]
	#layout = Any[1;1]
	savefig(plot!(entropy_main; title="(a)", titlelocation=:left, size=(600,400),subplot=1, dpi=200), joinpath(SAVEPATH, "entropy_new_part1.png"))
	savefig(plot!(p_lens; title="(b)", titlelocation=:left, size=(550,200), dpi=200), joinpath(SAVEPATH, "entropy_new_part2.png"))
	plot(plot!(entropy_main; title="(a)", subplot=1), 
		 plot!(p_lens,title="(b)"); 
		layout, dpi=200, size=(600,600), titlelocation=:left, margin_bottom=20Plots.mm)
	# layout = @layout [a 
	# 				  [b c]]
	# plot(entropy_main, p_linfit, p_lens; layout, dpi=200)
end

# ╔═╡ 9f63456f-9da1-4203-ad01-8d159ebd1bee
savefig(entropy_plot_full_W,joinpath(SAVEPATH, "hce_W.pdf"))

# ╔═╡ 9b4520ca-472e-4e1d-b612-62064f65889d
savefig(entropy_plot_full_W_new,joinpath(SAVEPATH, "hce_W_new.pdf"))

# ╔═╡ 300ae731-3551-43b8-ab40-924d521fdb7b
md"""
# Entropy Variance
"""

# ╔═╡ 33ec0a03-1d39-4b7f-91f2-55cc7c35b7fb
entropy_var_main_W_errors = let p = plot(;xlabel="W", ylabel="Half-chain entropy std.dev", legend=:topleft, xaxis=:log, xlim=W_xlim)
	xticks!(W_xticks, string.(W_xticks))
	color = palette(:ice, length(Ns)+2, rev=true)[2:end]
	plot!(p, W(ρvals), hce_std; palette=color, label=["N=$N" for N in Ns'], ribbon=hce_std_std)
	p
end

# ╔═╡ 0fe73110-2508-4e83-a959-fd40b41e9ecc
# use errors for fitting
entropy_var_full_parabola_W = let fitrange = 37:42,
	fitplotrange = fitrange,
	p = plot(entropy_var_main_W_errors; legend=:outside, legendposition=(1.0,0.0), dpi=200, tickfontsize=10, labelfontsize=12)
	@. model(x, p) = p[1] + p[2]*(x-p[3])^2

	fits = [curve_fit(model, W(ρvals[fitrange]), ycol[fitrange], 1 ./ yerrcol[fitrange], [1.0, -30, 0.7]) for (ycol,yerrcol) in zip(eachcol(hce_std), eachcol(hce_std_std))]
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
	#scatter!(p, maxima_loc, maxima_val; label="Maxima", marker=:x, color=:black)
	scatter!(p, maxima_loc, maxima_val; label="Maxima", marker=:vline, color=:black, xerror=maxima_std)

	plot!(p, Ns, maxima_loc; yerror=maxima_std,
	 	legend=nothing, xlabel="N",ylabel=L"W_{max}", inset=(1,bbox(0.65,0.05,0.5,0.35)), subplot=2)
	hline!(p, [1/(2*0.7475979202)];subplot=2,color=:black, ls=:dot,width=1,label="Reyni limit"	)
	annotate!(p, [(15.2, 0.658, ("Rényi limit",8))]; subplot=2)
	# inset = plot!(p, Ns, maxima_loc; yerror=maxima_std,
	# 	inset=(1,bbox(0.14,0.05,0.4,0.3)), subplot=2, legend=nothing, xlabel="N",ylabel="Maximum")

	p
end

# ╔═╡ c60ded9a-a363-4f99-96a9-9f073065a8c4
savefig(entropy_var_full_parabola_W, joinpath(SAVEPATH, "entropy_variance_W.pdf"))

# ╔═╡ c2f1cdc5-a6e8-4fbf-909f-bfb966709ddb
md"""
# Thouless Parameter
"""

# ╔═╡ ab2d8ca2-5fc8-4f47-ba17-9141b709eef5
md"""
## Load data
"""

# ╔═╡ 3f676374-ed14-4844-98fd-6a23eee6dfc4
tp_means_sz, tp_errorofmean_sz = let tp_shot_data_high = _prep_vals.(load_el.("sz", alpha6_high); dims=(1,2,4)), #[N][shot,rho]
	tp_shot_data_low =
	_prep_vals.(load_el.("sz", alpha6_low); dims=(1,2,4)) #[N][shot,rho]
	tp_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(tp_shot_data_low,tp_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean_ignore_nan.(tp_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(error_of_mean_ignore_nan.(tp_shot_data; dims=1))...)
	means, errorofmean
end;

# ╔═╡ 386d2101-cc03-4cbd-b986-8c5b1cd055b2
tp_means_szsz, tp_errorofmean_szsz = let tp_shot_data_high = _prep_vals.(load_el.("szsz", alpha6_high); dims=(1,2,4)), #[N][shot,rho]
	tp_shot_data_low =
	_prep_vals.(load_el.("szsz", alpha6_low); dims=(1,2,4)) #[N][shot,rho]
	tp_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(tp_shot_data_low,tp_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean_ignore_nan.(tp_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(error_of_mean_ignore_nan.(tp_shot_data; dims=1))...)
	means, errorofmean
end;

# ╔═╡ 0513fe1d-d908-447a-9ff7-88b03e9addf8
tp_means_hop, tp_errorofmean_hop = let tp_shot_data_high = _prep_vals.(load_el.("hopping", alpha6_high); dims=(1,2,4)), #[N][shot,rho]
	tp_shot_data_low =
	_prep_vals.(load_el.("hopping", alpha6_low); dims=(1,2,4)) #[N][shot,rho]
	tp_shot_data = [hcat(low, high)[:,ρsortperm] for (low,high) in zip(tp_shot_data_low,tp_shot_data_high)] #[N][shot,rho]
	means = hcat(vec.(mean_ignore_nan.(tp_shot_data; dims=1))...) #[rho,N]
	errorofmean = hcat(vec.(error_of_mean_ignore_nan.(tp_shot_data; dims=1))...)
	means, errorofmean
end;

# ╔═╡ e91c4684-7530-47b8-9ce4-12b06e235726
md"""
## Find transition location
"""

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
	(tp_sz_transition, tp_szsz_transition, tp_hop_transition)
end

# ╔═╡ 99dbe4d7-2919-4677-8f97-0566ff66b146
tp_N_window = 1

# ╔═╡ 01cb22ef-1d44-4910-8e0c-77da2411081d
begin
	tp_sz_transitions = hcat([tp_find_transition(tp_means_sz[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
	tp_szsz_transitions = hcat([tp_find_transition(tp_means_szsz[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
	tp_hop_transitions = hcat([tp_find_transition(tp_means_hop[tp_transition_window, (i-tp_N_window):(i+tp_N_window)], ρvals[tp_transition_window]) for i in (1+tp_N_window):(length(Ns)-tp_N_window)]...)'
end

# ╔═╡ 1b3e8683-02b7-45a1-a65e-cd74083bb2c5
tp_sz_transition

# ╔═╡ 30100b77-d704-4f7e-9aea-9a2953b86c3b
thouless_main_W = let p = plot(;xlabel="W", ylabel="Thouless parameter \$\\mathcal{G}\$", legend=:topright, xaxis=log, xlim=W_xlim)
	xticks!(p, W_xticks, string.(W_xticks))
	color = palette(:ice, length(Ns)+3; rev=true)[2:end]
	plot!(p, W(ρvals), tp_means_sz; palette=color, label=["N=$N" for N in Ns'])
	loc, width = tp_sz_transition
	vspan!(W([loc-width,loc+width]); label="Crossover", color=:darkred, alpha=0.7)
	p
end

# ╔═╡ 7122e488-d1cc-44c8-ac39-c92fa3b031fc
## errors only barely visible
thouless_main_W_errors = let p = plot(;xlabel="W", ylabel="Thouless parameter \$\\mathcal{G}\$", legend=:topright, xaxis=log, xlim=W_xlim)
	xticks!(p, W_xticks, string.(W_xticks))
	color = palette(:ice, length(Ns)+3; rev=true)[2:end]
	plot!(p, W(ρvals), tp_means_sz; palette=color, label=["N=$N" for N in Ns'], ribbon=tp_errorofmean_sz)
	loc, width = tp_sz_transition
	vspan!(W([loc-width,loc+width]); label="Crossover", color=:darkred, alpha=0.7)
	p
end

# ╔═╡ a50b42e4-2001-4dbc-8c67-9b914e543cfb
md"""
## Comparison of the three operators
"""

# ╔═╡ 41ce7ca3-6ea2-4132-b931-bfd4db9fa763
let sz_plot = plot(thouless_main_W_errors; title="\$W_1 = \\sigma_z\$"),
	szsz_plot = let p = plot(;title="\$W_2 = \\sigma_z \\sigma_z\$", xlabel="W", ylabel="Thouless parameter \$\\mathcal{G}\$", legend=:topright, xaxis=log, xlim=W_xlim)
		xticks!(p, W_xticks, string.(W_xticks))
		color = palette(:ice, length(Ns)+3; rev=true)[2:end]
		plot!(p, W(ρvals), tp_means_szsz; palette=color, label=["N=$N" for N in Ns'], ribbon=tp_errorofmean_szsz)
		loc, width = tp_szsz_transition
		vspan!(W([loc-width,loc+width]); label="Crossover", color=:darkred, alpha=0.7)
		p
	end
	hop_plot = let p = plot(;title="\$W_3 = \\sigma_+ \\sigma_- + \\mathrm{h.c.}\$", xlabel="W", ylabel="Thouless parameter \$\\mathcal{G}\$", legend=:topright, xaxis=log, xlim=W_xlim)
		xticks!(p, W_xticks, string.(W_xticks))
		color = palette(:ice, length(Ns)+3; rev=true)[2:end]
		plot!(p, W(ρvals), tp_means_hop; palette=color, label=["N=$N" for N in Ns'], ribbon=tp_errorofmean_hop)
		loc, width = tp_hop_transition
		vspan!(W([loc-width,loc+width]); label="Crossover", color=:darkred, alpha=0.7)
		p
	end
	plot(sz_plot,szsz_plot,hop_plot; layout=[1;1;1], size=(500,1000))
end

# ╔═╡ 83ec7e23-8f38-4523-b068-1e4a19aa02e1
md"""
Virtually identical
"""

# ╔═╡ 0e173da7-feb0-4c4e-aa34-084a756a168f
md"""
## Make and save plot
"""

# ╔═╡ 5a366799-c8a1-4fe5-afe2-909524e30e5e
thouless_full_W = plot(thouless_main_W; dpi=200, legend=:topright, legend_columns=-1, tickfontsize=10)

# ╔═╡ c3bb699c-e941-4a46-b430-bc1222b2ed09
savefig(thouless_full_W, joinpath(SAVEPATH, "thouless_W.pdf"))

# ╔═╡ f37acd5b-b8a4-4caf-91a0-9cba909e5a78
md"""
# Entropy per cut
"""

# ╔═╡ 5e099ce4-8f01-4087-bdb6-ba1d6ab23165
let p = plot(;xlabel="Number of pairs N",ylabel="Entropy per cut", dpi=200, labelfontsize=12, tickfontsize=10)
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
# ╟─2a544bad-35cb-4229-b75e-8d3283b8c288
# ╠═d05f72ed-31ac-4983-b510-ce6d6527879a
# ╠═a35a6fb5-dd1e-4998-8c2a-097acb9bf7d4
# ╠═1516ba16-a8ef-11ec-3979-1f7edd32dc6d
# ╠═4d171def-b37a-44ea-a192-03a88debe63c
# ╠═fc5ab3d2-acdd-47a7-bbdc-c40b28976907
# ╠═f92f65b3-59f0-4067-a744-3105c9411d0a
# ╠═ed2305f5-7bc9-43a8-b5eb-75b5a8c00eca
# ╠═da20058c-4073-4ac6-badd-00adfa0860eb
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
# ╟─799f75b3-6bbd-4325-89b2-6a24f780f51a
# ╠═ac92e5e6-171e-4b8e-aca3-344c0d8ea1a9
# ╠═33b9ecab-6984-416a-8127-6367da4d04ef
# ╟─336c5f98-71a9-40e8-9821-7c56033a40f9
# ╠═f0665d91-1f67-4716-87ce-a3ea3b9039bf
# ╠═a5c221a6-c6c2-4301-b93b-7178859bcf68
# ╠═0a4e0811-036d-46f5-ae89-8ffa7fb7ab1a
# ╠═8b8bafa7-e325-4d36-9624-1e1efb6c972d
# ╠═1301fc69-9d92-43c6-8274-eb3c11995b8b
# ╟─7f082cf1-12e8-4134-a6e3-371af39164c7
# ╟─81071ead-3233-42bc-a7cb-c0aa4fc33ffe
# ╠═80e82ba2-8a72-4ba9-9828-42aea93ba3bd
# ╠═4721cdab-edc9-42dc-9656-64534000ed35
# ╠═5b30d6d4-6b12-4da3-b5f5-1e1170abaac5
# ╠═bbea8da3-c3a5-440a-875d-99ea1293e546
# ╠═b1e54703-202c-47a1-8a15-c8f4e86e9387
# ╠═c4eae1f2-95c1-4114-a1dd-3b437414c189
# ╠═aef7cc15-7c6a-4f17-b6aa-5d6398ba8487
# ╟─9db9970d-436b-46d5-bc9b-65477dc529f1
# ╠═5e57efdc-5ae3-4a86-ae6d-9f0a49c4d5c7
# ╠═21352df0-3493-48ef-865e-2f6b0614b0f4
# ╠═aad4eeb4-8fdb-40e3-9301-49bd9a4de439
# ╠═241ddd4b-f075-4785-b161-2192982f2137
# ╠═fc0a8c34-1bc8-4003-8165-13a362257648
# ╠═81fcb0e6-ad67-4313-8af5-a6d79d1e941f
# ╟─97bf13c5-cfc2-4554-925d-759444bf3a45
# ╟─0f802771-29b5-472d-bcee-c6987bd97aea
# ╠═eb5b02be-5195-47b0-a300-f4071e0f56c3
# ╠═b3ca4a97-d758-4fa5-b881-fc34063b5a28
# ╠═6ea0c192-202f-472a-a7f0-51a1a8e72256
# ╠═05ccfc6e-92c2-45b5-b49d-ead83bf1ad4c
# ╠═9f8782a6-8c0b-45cb-94ec-64a07addb55b
# ╠═895cace4-ae8b-4d39-8b76-3cdf09815940
# ╠═a8b925f1-1763-43eb-9162-219b85ede984
# ╟─0504d150-3841-4a44-9d06-f5d27b3f9b39
# ╠═99d387f0-ae7e-4ce4-8ad8-6c8bf95e3819
# ╠═bbf4bd03-23eb-4dd9-b91c-d92eb3c4ff6a
# ╠═5d930be7-5df7-47bf-be84-2a055106920b
# ╠═f9abd6af-f125-41ef-915d-a73b041f9690
# ╠═f73cc7af-6f60-48c2-bdf1-3de67260e9b6
# ╠═6b5ad674-b541-45ec-8bbf-76543f77297f
# ╠═bfc60db7-04c1-4934-a3de-c6b18d841a2f
# ╠═9f63456f-9da1-4203-ad01-8d159ebd1bee
# ╠═9b4520ca-472e-4e1d-b612-62064f65889d
# ╟─300ae731-3551-43b8-ab40-924d521fdb7b
# ╠═33ec0a03-1d39-4b7f-91f2-55cc7c35b7fb
# ╠═0fe73110-2508-4e83-a959-fd40b41e9ecc
# ╠═c60ded9a-a363-4f99-96a9-9f073065a8c4
# ╟─c2f1cdc5-a6e8-4fbf-909f-bfb966709ddb
# ╟─ab2d8ca2-5fc8-4f47-ba17-9141b709eef5
# ╠═3f676374-ed14-4844-98fd-6a23eee6dfc4
# ╠═386d2101-cc03-4cbd-b986-8c5b1cd055b2
# ╠═0513fe1d-d908-447a-9ff7-88b03e9addf8
# ╟─e91c4684-7530-47b8-9ce4-12b06e235726
# ╠═b5d37d9f-7c44-4218-9e0b-1cb6e4915ad3
# ╠═f69d1bb5-20cf-416f-a581-4acb56676783
# ╠═1376a743-1355-49a5-977d-264e63bc0881
# ╠═5d1d0e06-5269-4b65-ba36-c71ef6afaa52
# ╠═99dbe4d7-2919-4677-8f97-0566ff66b146
# ╠═01cb22ef-1d44-4910-8e0c-77da2411081d
# ╠═1b3e8683-02b7-45a1-a65e-cd74083bb2c5
# ╠═30100b77-d704-4f7e-9aea-9a2953b86c3b
# ╠═7122e488-d1cc-44c8-ac39-c92fa3b031fc
# ╟─a50b42e4-2001-4dbc-8c67-9b914e543cfb
# ╠═41ce7ca3-6ea2-4132-b931-bfd4db9fa763
# ╟─83ec7e23-8f38-4523-b068-1e4a19aa02e1
# ╟─0e173da7-feb0-4c4e-aa34-084a756a168f
# ╠═5a366799-c8a1-4fe5-afe2-909524e30e5e
# ╠═c3bb699c-e941-4a46-b430-bc1222b2ed09
# ╟─f37acd5b-b8a4-4caf-91a0-9cba909e5a78
# ╠═5e099ce4-8f01-4087-bdb6-ba1d6ab23165
