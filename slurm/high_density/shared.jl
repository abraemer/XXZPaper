println("zero-field with alpha=6")

## environment
import Pkg
Pkg.activate(".")

if haskey(ENV, "ON_CLUSTER")
    @eval using MKL
    println("Added MKL.jl")
end

using LinearAlgebra # for BLAS threads

println("Working Directory:          $(pwd())" )
println("Running on host:            $(gethostname())" )
println("Job id:                     $(get(ENV, "SLURM_JOB_ID", ""))" )
println("Job name:                   $(get(ENV, "SLURM_JOB_NAME", ""))" )
println("Number of nodes allocated:  $(get(ENV, "SLURM_JOB_NUM_MODES", ""))" )
println("Number of cores allocated:  $(get(ENV, "SLURM_NTASKS", ""))" )
println("#threads of Julia:          $(Threads.nthreads())")
println("#threads of BLAS:           $(BLAS.get_num_threads())")
println("#BLAS config:               $(BLAS.get_config())")

@show ARGS

Pkg.instantiate(; io=stdout)
Pkg.status(; io=stdout)

using SimLib, XXZNumerics, SpinSymmetry, Random

## constants and ARGS
const RUNNUMBER = parse(Int, ARGS[3])
const PREFIX = joinpath(path_prefix("xxzpaper"), "high_density/run_$RUNNUMBER")
const GEOMETRY = :noisy_chain_pbc
const DIM = 1
const ALPHA = 6
const N = parse(Int, ARGS[1])
const ρs = 1 ./ (0.5025:0.025:1.2) # ρ=1.99...0.83
const SHOTS = parse(Int, ARGS[2])
const FIELDS = [0.0]
const BLOCK = div(N-1,2)

const BASIS = SymmetrizedBasis(zbasis(N, BLOCK), [], [])
const DIAGTYPE = Full() # TODO
const RUNMODE = haskey(ENV, "ON_CLUSTER") ? Parallel(48, true, 2, true) : Parallel(4, true, 2, false)

@show PREFIX
@show GEOMETRY
@show DIM
@show N
@show ALPHA
@show ρs
@show SHOTS
@show FIELDS
@show BLOCK
@show BASIS
@show DIAGTYPE
@show RUNMODE

const LOCATION = SaveLocation(;prefix=PREFIX)




pdd = PositionDataDescriptor(GEOMETRY, DIM, N, SHOTS, ρs, LOCATION)
if !isfile(datapath(pdd))
    logmsg("No positions generated yet! Generating with seed=$(RUNNUMBER).")
    Random.seed!(RUNNUMBER)
    save(create(pdd))
end

model = RandomPositionsXXZWithXField(pdd, PowerLaw(ALPHA), FIELDS, :none, BASIS)

edd = EDDataDescriptor(model, DIAGTYPE, LOCATION)

hopping_operator = 1/2 * real.(σx⊗σx + σy⊗σy) ⊗ identity_op(N-2)

prtasks = isodd(N) ? [] :  [ParticipationRatio(BASIS, PairBasis(N,PetersMethod())),
                            ParticipationRatio(BASIS, PairBasis(N,NaivePairing(true))),
                            ParticipationRatio(BASIS, PairBasis(N,NaivePairing(false)))]

tasks = [Energies(),
        OperatorDiagonal("hopping", symmetrize_operator(hopping_operator, BASIS)),
        EigenstateOccupation("xpol", symmetrize_state(normalize!(ones(2^N)), BASIS)),
        LevelSpacingRatio(),
        EigenstateLocality("sz", symmetrize_operator(single_spin_op(σz, 1, N), BASIS)),
        EigenstateLocality("hopping", symmetrize_operator(hopping_operator, BASIS)),
        EigenstateLocality("szsz", symmetrize_operator(correlator(σz, 1, 2, N), BASIS)),
        HalfChainEntropyZBlock(L=1, basis=BASIS),
        HalfChainEntropyZBlock(L=2, basis=BASIS),
        HalfChainEntropyZBlock(L=div(N, 3), basis=BASIS),
        HalfChainEntropyZBlock(L=div(N, 2), basis=BASIS),
        ParticipationRatio(),
        prtasks...]

@show typeof.(tasks)

logmsg("*"^10 * "run ed" * "*"^10; doflush=true)
@time edata = run_ed(edd, tasks, RUNMODE)

logmsg("*"^10 * "Saving" * "*"^10)
save.(edata)

using Statistics
const AVG_PREFIX = joinpath(path_prefix("xxzpaper"), "high_density-avg/run_$RUNNUMBER")
for (data, datatype) in zip(edata[4:end], [LSRData, ELData,ELData,ELData, HCEData,HCEData,HCEData,HCEData, PRData,PRData,PRData,PRData])
    avg = datatype(data.descriptor, mean(data.data; dims=1))
    avg.descriptor.derivedfrom.pathdata.prefix = AVG_PREFIX
    save(avg)
end
