
#Global Parameters
TRAJ_STATE=NVE
STATE=PI
SOLVENT=DBPC70
d_resname=DBP
a_resname=C70
solvent=${SOLVENT}
platform=cpu


#Job parameter
system=DBP-C70_
stripsolvent=no
run_cpu_min=no
run_gpu_MD=no
run_cpu_rec=no
run_energy_analysis=yes
START=
END=
total_traj=1
continue_from=0
inpdir=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/MD/
GIVEN_STRUC=conf11
DEP=9697759,9697760,9697761
state_list='GR LE CT1 CT2 CT3 CT4 CT5 CT6 '
for structure in $GIVEN_STRUC #(( structure = ${START}; structure <= ${END}; structure++))
do

PROJECT_DIR=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/MD/
DATA_DIR=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/results/DBP-C70_conf11/QC/
cpu_script=cpujob_control.slurm
gpu_script=gpujob_control.slurm
SRCDIR=/gpfsnyu/xspace/projects/PCT-example/DBP-C70/template
JOBDIR=${PROJECT_DIR} #/${system}${structure}

mkdir -p ${JOBDIR}

#Setup slurm script for the run
# job option : min , heat , eq_npt, eq_nvt, sample_nvt, relax_nve, prod_nve, recalculate 
# in cpu script default is all no 

# other option : for getting more trajectory chagnge NTRAJ as needed, default NTRAJ=1
#inputfile is define at the section min_inputfile, heat_inputfile and so on , full documentation should follow

if [ "$stripsolvent" == "yes" ];then
module use /gpfsnyu/xspace/sungroup/modules
module load amber/22z

# Go to DATA_DIR and create a .prmtop without THF for all state
cd ${PROJECT_DIR}
for Charge_State in $state_list #GR CT PI
do
cat > processtop.cpptraj <<EOF
parm ${system}${structure}_${solvent}_${Charge_State}.prmtop
parmstrip :${d_resname}
parmstrip :${a_resname}
parmwrite out ${system}${structure}_no${solvent}_${Charge_State}.prmtop
EOF

cpptraj -i processtop.cpptraj
echo .prmtop for state ${Charge_State} is successfully strip from THFs
done
fi



# get all necessary files to start a JOBDIR
#cp ${SRCDIR}/${cpu_script} ${SRCDIR}/${gpu_script} ${JOBDIR}
cp ${DATA_DIR}/*.prmtop ${DATA_DIR}/*.inpcrd ${JOBDIR}
cd ${JOBDIR}




if [ "$run_cpu_min" == "yes" ];then
echo Setting up cpu job run for structure no $structure ...
sed_param="s#THE_DIR=.*#THE_DIR="$JOBDIR"#"
sed -i $sed_param ${cpu_script}
sed_param="s|HQDIR=.*|HQDIR="$inpdir"|"
sed -i "$sed_param"  ${cpu_script}
# For minimization 
sed_param=s/SYSTEM=.*/SYSTEM="${system}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/ID_NO=.*/ID_NO="${structure}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/min=.*/min="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/NTRAJ=.*/NTRAJ="${total_traj}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/CONTINUE_FROM=.*/CONTINUE_FROM="${continue_from}"/
sed -i "$sed_param" ${cpu_script}

sed_param=s/heat=.*/heat="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/eq_npt=.*/eq_npt="yes"/
sed -i "$sed_param" ${cpu_script}
echo "slurm script successfully edited, submitting now"

if [ -z "$DEP" ];
then
jobid=$(sbatch --parsable ${cpu_script})
echo Submitted job id $jobid
else
jobid=$(sbatch --parsable --dependency=afterok:$DEP ${cpu_script})
echo Submitted job id $jobid
fi

#sbatch ${cpu_script}

fi

#########################################################
# for GPU based MD run
if [ "$run_gpu_MD" == "yes" ];then
echo Setting up ${platform}  MD run for structure no $structure ...
sed_param="s#THE_DIR=.*#THE_DIR="$JOBDIR"#"


script=${platform}job_control.slurm

sed -i $sed_param ${script}
sed_param="s|HQDIR=.*|HQDIR="$inpdir"|"
sed -i "$sed_param"  ${script}
# For minimization 
sed_param=s/ID_NO=.*/ID_NO="${structure}"/
sed -i "$sed_param" ${script}
sed_param=s/NTRAJ=.*/NTRAJ="${total_traj}"/
sed -i "$sed_param" ${script}
sed_param=s/CONTINUE_FROM=.*/CONTINUE_FROM="${continue_from}"/
sed -i "$sed_param" ${script}
sed_param=s/SYSTEM=.*/SYSTEM="${system}"/
sed -i "$sed_param" ${script}


echo "slurm script successfully edited, submitting now"

#sbatch ${gpu_script}

if [ -z "$DEP" ];
then
jobid=$(sbatch --parsable ${script})
echo Submitted job id $jobid
else
jobid=$(sbatch --parsable --dependency=afterok:$DEP ${script})
echo Submitted job id $jobid
fi

fi
########################################################
# for recalculate
if [ "$run_cpu_rec" == "yes" ];then
echo Setting up cpu recalculation job run for structure no $structure ...
sed_param="s#THE_DIR=.*#THE_DIR="$JOBDIR"#"
sed -i $sed_param ${cpu_script}
sed_param="s|HQDIR=.*|HQDIR="$inpdir"|"
sed -i "$sed_param"  ${cpu_script}
# For minimization 
sed_param=s/recalculate=.*/recalculate="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/SYSTEM=.*/SYSTEM="${system}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/ID_NO=.*/ID_NO="${structure}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/NTRAJ=.*/NTRAJ="${total_traj}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/CONTINUE_FROM=.*/CONTINUE_FROM="${continue_from}"/
sed -i "$sed_param" ${cpu_script}

echo "slurm script successfully edited, submitting now"

if [ -z "$DEP" ];
then
jobid=$(sbatch --parsable ${cpu_script})
echo Submitted job id $jobid
else
jobid=$(sbatch --parsable --dependency=afterok:$DEP ${cpu_script})
echo Submitted job id $jobid
fi

#sbatch ${cpu_script}

fi
#######################################################
# for energy_analysis
if [ "$run_energy_analysis" == "yes" ];then
echo Setting up cpu energy_analysis job run for structure no $structure ...
sed_param="s#THE_DIR=.*#THE_DIR="$JOBDIR"#"
sed -i $sed_param ${cpu_script}
sed_param="s|HQDIR=.*|HQDIR="$inpdir"|"
sed -i "$sed_param"  ${cpu_script}
# For minimization
sed_param=s/ID_NO=.*/ID_NO="${structure}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/SYSTEM=.*/SYSTEM="${system}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/energy_analysis=.*/energy_analysis="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/get_correction=.*/get_correction="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/run_on_tmp=.*/run_on_tmp="yes"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/NTRAJ=.*/NTRAJ="${total_traj}"/
sed -i "$sed_param" ${cpu_script}
sed_param=s/CONTINUE_FROM=.*/CONTINUE_FROM="${continue_from}"/
sed -i "$sed_param" ${cpu_script}


echo "slurm script successfully edited, submitting now"

#sbatch ${cpu_script}
if [ -z "$DEP" ];
then
jobid=$(sbatch --parsable ${cpu_script})
echo Submitted job id $jobid
else
jobid=$(sbatch --parsable --dependency=afterok:$DEP ${cpu_script})
echo Submitted job id $jobid
fi

fi



done

