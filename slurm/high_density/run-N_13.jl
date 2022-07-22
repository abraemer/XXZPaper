#!/bin/sh
# ########## Begin Slurm header ##########
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:13:00
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=48
#SBATCH --job-name=paper-N_13-hd
#SBATCH --output="logs/high_density-N_13-%j.out"
########### End Slurm header ##########
#=
# load modules
# not needed - julia installed locally

# export JULIA_DEPOT_PATH=$SCRATCH
export ON_CLUSTER=1
exec julia --color=no --threads=1 --startup-file=no slurm/high_density/shared.jl 13 2000 $1
=#
