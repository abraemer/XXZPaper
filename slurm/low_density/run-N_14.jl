#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:30:00
#SBATCH --mem=50gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=paper-N_14-ld
#SBATCH --output="logs/low_density-N_14-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=1 --startup-file=no slurm/low_density/shared.jl 14 2000 $1
=#
