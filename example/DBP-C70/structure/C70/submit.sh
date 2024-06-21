#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 40
#SBATCH --time=24:00:00
#SBATCH --mem=100gb
#SBATCH --job-name="testgs"
#SBATCH --output=slurm_%j.o
#SBATCH --error=slurm_%j.e
#SBATCH --partition=parallel

module use /gpfsnyu/xspace/sungroup/modules
module load gaussian16

g16 resp.com
#/gpfsnyu/xspace/sungroup/modules
