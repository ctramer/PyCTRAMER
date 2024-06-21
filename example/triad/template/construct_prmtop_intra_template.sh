export MODULEPATH=$MODULEPATH:/gpfsnyu/xspace/sungroup/modules
# for e5 series cpu node (compute1-27), load avx2 version
if [ -n ${SLURM_JOB_NODELIST:7} ] && [ ${SLURM_JOB_NODELIST:7} -lt 28 ];then
module load amber/20z_avx2
else
# for gold/platinum series cpu node, load avx512 version
module load amber/20z
fi
project_dir=
system=
chargefile= #=charge0.txt
#code=TRIvvi
code=
solvent= #DDDAAA
solute=
solvent_mol=2700
N_layer=
solvent_nameofsolute=TRI
solu_resname=
D_resname=
A_resname=
d_resname=
a_resname=
box_side=100
com=50
forcefield=gaff
solvent_dir=
D_dir=
A_dir=
d_dir=
a_dir=
functional=BNL
SRCDIR=
START=
END=
GIVEN_STRUC=
STATE_NUM=
omega=114514
N_layer=25
#-----------------------------------------
# this is to execute looping directly here 
 for structure in $GIVEN_STRUC  #(( structure = ${START}; structure <= ${END}; structure++))
 do

 for states in $STATE_NUM
 do 

 for w in $omega 
 do 
#-----------------------------------------
#WORKDIR=$project_dir/${system}${structure}_${functional}_${states}states_w${w}
WORKDIR=
#1. Getting Necessary files for solvent and zero charge
cd ${SRCDIR}
cp ${chargefile} ${WORKDIR}

#cd ${solvent_dir}
#cp *THF* ${WORKDIR}
cd ${WORKDIR}


# 2. from pdb to mol2
# antechamber -i ${system}${structure}.pdb -fi pdb -o ${system}${structure}.mol2 -fo mol2 -nc 0 -cf ${chargefile} -rn ${code} -pl 6

# break pdb into 2 pdbs 


# REMINDER: name of the residue in the D, A, a, and d moiety should be the consistent with chemical structure;
# for example, D and d are DBP in different positions, then we should name the residue in the same name such as DBP; 
# similarly, A and a in this case is Fullerene, C70, and we can name it as Ful. 

# 3. from mol2 to frcmod
# parmchk2 -i ${system}${structure}.mol2 -f mol2 -o ${system}${structure}.frcmod

# 4. packmol solvation

cat > packmol.inp <<EOF
## A DDDAAA heterojunction box with Box length of 60A, user can construct 
# there own heterojunction structure at nanoscale 

## All the atoms from different molecules will be separated at least 2.0

tolerance 2.0

# The file type and input type

output ${system}${structure}_${solvent}.pdb
filetype pdb
structure ${D_dir}/${system}${structure}.pdb
  number 1
  centerofmass
  fixed ${com}. ${com}. ${com}. 0. 0. 0.
end structure

structure ${d_dir}/${solvent}.pdb
  number ${solvent_mol}
  inside box 0. 0. 0. ${box_side}. ${box_side}. ${box_side}.
end structure


EOF

packmol < packmol.inp

# 5. from frcmod to prmtop

cat > tleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.${forcefield}

${D_resname} = loadmol2  ${D_dir}/${D_resname}.mol2    # ${WORKDIR}/${solvent}.mol2
loadamberparams     ${D_dir}/${D_resname}.frcmod #${WORKDIR}/${solvent}.frcmod

${d_resname} = loadmol2   ${d_dir}/${d_resname}.mol2
loadamberparams   ${d_dir}/${d_resname}.frcmod

check ${d_resname}
check ${D_resname}

com=loadpdb ${WORKDIR}/${system}${structure}_${solvent}.pdb
setbox com centers

savepdb com ${system}${structure}_solvated_in_${solvent}.pdb
saveamberparm com ${system}${structure}_${solvent}.prmtop ${system}${structure}_${solvent}.inpcrd
quit

EOF

tleap -f tleap.in 

done 
done
done 



