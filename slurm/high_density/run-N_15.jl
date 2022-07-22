#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=29:10:00
#SBATCH --mem=140gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=paper-N_15-hd
#SBATCH --output="logs/high_density-N_15-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=1 --startup-file=no slurm/high_density/shared.jl 15 800 $1
=#
