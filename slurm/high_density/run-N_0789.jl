#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=paper-N_789-hd
#SBATCH --output="logs/high_density-N_789-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
julia --color=no --threads=1 --startup-file=no slurm/high_density/shared.jl 7 2000 $1
julia --color=no --threads=1 --startup-file=no slurm/high_density/shared.jl 8 2000 $1
julia --color=no --threads=1 --startup-file=no slurm/high_density/shared.jl 9 2000 $1
=#
