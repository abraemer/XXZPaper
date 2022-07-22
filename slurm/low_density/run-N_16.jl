#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=61:00:00
#SBATCH --mem=350gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=paper-N_16-ld
#SBATCH --output="logs/low_density-N_16-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=1 --startup-file=no slurm/low_density/shared.jl 16 400 $1
=#
