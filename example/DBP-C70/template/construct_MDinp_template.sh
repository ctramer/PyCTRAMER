#!/bin/bash 

Template_Dir=
Target_Dir=

# Warning! Warning! Warning!
# one need to use different template dirs with different simulation parameters in order to avoid misuse of MD input files 

MD_steps=      
MD_dt=          
sample_steps= 
MD_temperature=

# 

cd ${Template_Dir}


for simu in  heat  equil_NPT  equil_NVT  nvt_eq_fixed_solute  nvt_sampling_fixed_solute  nve_prod_fixed_solute  nve_relax_fixed_solute
do

ori_file=${Template_Dir}/${simu}.temp 
tar_file=${simu}.in 

sed     "s/MYDT/${MD_dt}/"  $ori_file >    $tar_file
sed -i  "s/TARTEMP/${MD_temperature}/g"     $tar_file
sed -i  "s/MYSAMPLESTEP/${sample_steps}/g"  $tar_file
sed -i  "s/MYMDSTEP/${MD_steps}/"          $tar_file

cp $tar_file $Target_Dir

done  

cp ${Template_Dir}/min.in  $Target_Dir
