#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --mem=100gb
#SBATCH --job-name="PCT-triad"
#SBATCH --mail-type=END
#SBATCH --mail-user=netid@nyu.edu
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --partition=argon
#SBATCH -c 96
##SBATCH --qos=argon
##SBATCH --constraint=g6132

TEMPDIR=/tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
mkdir -p $TEMPDIR
case=conf3
module use /gpfsnyu/xspace/sungroup/modules
module load qchem
module load openmpi/gnu/1.10.7 gcc/7.3
source /gpfsnyu/packages/qchem/4.4/qcenv.sh
export QC=/gpfsnyu/packages/qchem/4.4
export QCSCRATCH=$TEMPDIR
export QCAUX=/gpfsnyu/packages/qchem/4.4/qcaux

qchem -nt 96 triad_conf3.inp > triad_conf3.out


#remove TEMPDIR
rm -rf $TEMPDIR

cd /gpfsnyu/xspace/projects/PCT-example/triad/results/triad_conf3/QC/
