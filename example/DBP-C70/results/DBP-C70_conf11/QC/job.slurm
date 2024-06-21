#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --mem=180gb
#SBATCH --job-name="DBPC70"
#SBATCH --mail-type=END
#SBATCH --mail-user=xs13@nyu.edu
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --partition=parallel
#SBATCH -c 40
##SBATCH --qos=argon
##SBATCH --constraint=g6132

TEMPDIR=/tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
mkdir -p $TEMPDIR
case=conf11
module use /gpfsnyu/xspace/sungroup/modules
module load qchem/cpu/6.0.1
source /gpfsnyu/xspace/sungroup/software/qchem6cpu/qcenv.sh


export QCSCRATCH=$TEMPDIR


qchem -nt 40 DBP-C70_conf11.inp > DBP-C70_conf11.out


#remove TEMPDIR
rm -rf $TEMPDIR

cd /gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/QC/
