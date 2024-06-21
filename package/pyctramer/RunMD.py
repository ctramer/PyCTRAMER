from .Utilities import * 
from shutil import copyfile

def run_step4(i, job_control, dict_of_simulation, job_info):
    #print(arrow_str,"Started Step 4: RunMD")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] #caseid = str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    template_dir =  dict_of_simulation['template_dir']
    # 
    platform = 'cpu'
    jobtype = 'minimization'
    total_traj= dict_of_simulation['N_traj']
    PARM = ''
    MD_dyn_state = dict_of_simulation['MD_traj_state']
    MD_Donor = dict_of_simulation['MD_traj_state']
    inp_dir = md_dir # dict_of_simulation['template_dir']
    Solvent = dict_of_simulation['solvent_resname']
    MD_dt = dict_of_simulation['MD_dt']
    sample_steps = dict_of_simulation['sample_steps']
    MD_temperature  = dict_of_simulation["MD_temperature"] 
    MD_steps  = dict_of_simulation["MD_steps"]

    generate_MDinp(template_dir, MD_dt, MD_steps, sample_steps, MD_temperature, md_dir)

    state_list_str = "'"
    for info in job_info:
        for state_name in info['State_list']:
            state_list_str = state_list_str + state_name + ' ' 
    state_list_str = state_list_str + "'"  
    
    
    jobs = []
    newline = ''
    
    jobtype = 'minimization'
    platform = 'cpu'
    returning_info = setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir, platform, 
                                             jobtype, project, caseid, total_traj, PARM, MD_dyn_state, MD_Donor, 
                                             Solvent, state_list_str,
                                             continue_from=0, dependency=newline)#@newline)
    print(returning_info)

    if len(jobs) == 0:
        dependency_at_this_step = '' 
    else:
        dependency_at_this_step = list_to_dep_line(jobs)
    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
     
    jobtype = 'heat_equil_run'
    platform = 'gpu'
    returning_info = setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir,
                                              platform, jobtype, project, caseid, total_traj, PARM,
                                              MD_dyn_state, MD_Donor, Solvent,state_list_str, 
                                              continue_from=0, dependency=dependency_at_this_step) #newline)
    print(returning_info)
    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
    
    jobtype='recalculate'
    platform = 'cpu'
    returning_info = setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir, platform, 
                                             jobtype, project, caseid, total_traj, PARM, MD_dyn_state,
                                               MD_Donor, Solvent, state_list_str,
                                               continue_from=0, dependency=dependency_at_this_step) #newline)
    print(returning_info)
    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
    # print("RunMD job dependency: ",dependency_at_this_step, jobs)
    
    # jobtype = 'energy_analysis'
    jobtype = 'strip_solvent' 
    returning_info = setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir, 
                                             platform, jobtype, project, caseid, total_traj, 
                                             PARM, MD_dyn_state, MD_Donor, Solvent, 
                                             state_list_str, continue_from=0, dependency=dependency_at_this_step)#newline)
    print(returning_info)
    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
    
    job_control[i]['MD_jobids'] = dependency_at_this_step

    # jobtype = 'strip_solvent' 
    jobtype = 'energy_analysis'
    returning_info = setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir, platform, jobtype, project, caseid, 
                                             total_traj, PARM, MD_dyn_state, MD_Donor, Solvent, 
                                             state_list_str,
                                             continue_from=0, dependency=dependency_at_this_step)
    
    find_job_id(returning_info, jobs)
    dependency_at_this_step = list_to_dep_line(jobs)
    job_control[i]['MD_jobids'] = dependency_at_this_step
    print(returning_info)

    stroutput = "\n" + arrow_str
    stroutput = stroutput + "STEP 4: RunMD\n"
    stroutput = stroutput + "Input: \n" + md_dir + "\n"
    stroutput = stroutput + md_dir + "/*.in\n*=min,heat,equil_NPT,equil_NVE,sample_NVT,relax_NVE,prod_NVE\n" 
    stroutput = stroutput + "Job submit script: \n" + md_dir + "cpujob_control.slurm\n"
    stroutput = stroutput + md_dir + "gpujob_control.slurm\n"
    stroutput = stroutput + md_dir + "master_control.slurm\n"
    stroutput = stroutput + md_dir + "master_control.slurm\n"
    stroutput = stroutput + "Output: \n" + md_dir + "\n"
    stroutput = stroutput + md_dir + "/Traj_"+ MD_dyn_state + "/*.out \n\n" 


    stroutput = stroutput + arrow_str + "MD simulation details \n" 
     
    stroutput = stroutput + 'N_traj:'  + total_traj + '\n'
    stroutput = stroutput + 'MD_traj_state:'  + MD_dyn_state + '\n'
    stroutput = stroutput + 'solvent_resname:'  + Solvent + '\n'
    stroutput = stroutput + 'MD_dt:'  + MD_dt + '\n'
    stroutput = stroutput + 'sample_steps:'  + sample_steps + '\n'
    stroutput = stroutput + "MD_temperature:" + MD_temperature + '\n'
    stroutput = stroutput + "MD_steps:"   + MD_steps + '\n'
    #stroutput =  stroutput + '\n\n'
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 

    return job_control


def setup_MD_master_control(template_dir, qc_dir, md_dir, inp_dir, platform, jobtype, case, testing_string, 
                            total_traj, PARM, MD_dyn_state, MD_Donor, Solvent,state_list, continue_from=0,dependency=''):
    """
    # submitQC(template_dir,running_dir,QCmethod,basis_set,testing_str,Drange,Arange,N_cis,case,basedir):

    To setup a master control for the ctramer MD job.
    
    template_dir : is slurm template directory where the master_control_template.sh and master_control.sh are located.
    
    project_dir = main project data directory , 
        the quantum chemistry data will be read from and written to project_dir/case/QC/
        the molecular dynamics data will be read from and written to project_dir/case/MD/
        data directory will be created in project_dir/case/QC/case_${structure_id}_PARM
    
    MD_dyn_state: state for MD propagation, be a string
    QCDir: where to store the FF file
    
    platform = gpu or cpu 
        define the platform for which the job will run
        
    jobtype = 'strip_solvent','minimization', 'heat_equil_run', 'recalculate', 'energy_analysis'  
        if platform == gpu only heat_equil_run can be performed
        
    case : is the case being considered, will also be used for naming folder under the project_dir 
        as well as many other important file choose wisely
    start,end : int, integer index for structure identification
            [x,y]
        
    PARM = string identification for the job parameter, will be used for naming.
    
    continue_from :  int, default = 0.
        in almost all cases you will start from 0 , if for some reason you want to repeat some of the job, 
        you can start where to start with
        

    """
    #General setup for master control file
    master_control_template =f'{template_dir}/master_control_template.sh'
    master_control_file =f'{md_dir}/master_control.sh'
    copyfile(master_control_template,master_control_file)
    
    JOB_template=f'{template_dir}/{platform}job_control_template.slurm'
    JOB_FILE=f'{md_dir}/{platform}job_control.slurm'
    copyfile(JOB_template,JOB_FILE)
    # print("master_control_template,master_control_file,JOB_template,JOB_FILE",master_control_template,master_control_file,JOB_template,JOB_FILE)
    # variable setting and naming
    system= case+'_'
    
    #Setting up the job content
    
    # a lof of if functions 
    MD_DIR= md_dir#$ project_dir + '/' + case + '/MD/'
    QC_DIR= qc_dir #project_dir + '/' + case + '/QC/'+system+'${structure}_'+PARM


    if jobtype == 'strip_solvent':
        job = ['yes','no','no','no','no']
        
    elif jobtype == 'minimization':
        job = ['no','yes','no','no','no']
    
    elif jobtype == 'heat_equil_run':
        job = ['no','no','yes','no','no']
        
    elif jobtype == 'recalculate':
        job = ['no','no','no','yes','no']
    
    elif jobtype == 'energy_analysis':
        job = ['no','no','no','no','yes']
    
    else :
        print('Support for other jobtype is under construction')
    print(arrow_str,jobtype)
    sed_inplace(master_control_file,'PROJECT_DIR=','PROJECT_DIR='+MD_DIR)
    #sed_inplace(master_control_file,'JOBDIR=','JOBDIR='+MD_DIR)
    sed_inplace(master_control_file,'SRCDIR=','SRCDIR='+ template_dir)
    sed_inplace(master_control_file,'DATA_DIR=','DATA_DIR='+QC_DIR)
    #sed_inplace(master_control_file,"ID_NO=",   "ID_NO="+testing_string)
    # sed_inplace(master_control_file,'SYSTEM=','SYSTEM='+system)
    print('system  ', system)
    if len(Solvent.split()) > 1:
        print('Solvent.split(',')[0] + Solvent.split(',')[1]', Solvent.split(',')[0] + Solvent.split(',')[1])
    if len(Solvent.split(",")) == 1:
        sed_inplace(master_control_file,'SOLVENT=','SOLVENT='+Solvent)
    else:
        sed_inplace(master_control_file,'d_resname=','d_resname='+Solvent.split(',')[0])
        sed_inplace(master_control_file,'a_resname=','a_resname='+Solvent.split(',')[1])
        sed_inplace(master_control_file,'SOLVENT=','SOLVENT='+ Solvent.split(',')[0] + Solvent.split(',')[1]) 
        print('test solvent ', master_control_file)
    sed_inplace(master_control_file,'system=','system='+system)    
    #clst.sed_inplace(master_control_file,'START=','START='+str(start))
    #clst.sed_inplace(master_control_file,'END=','END='+str(end))
    sed_inplace(master_control_file,'GIVEN_STRUC=','GIVEN_STRUC='+testing_string)
    sed_inplace(master_control_file,'total_traj=','total_traj='+str(total_traj))
    sed_inplace(master_control_file,'continue_from=','continue_from='+str(continue_from))
    sed_inplace(master_control_file,'stripsolvent=','stripsolvent='+job[0])
    print('platform', platform)
    sed_inplace(master_control_file, 'platform=', 'platform='+ platform)
    sed_inplace(master_control_file,'run_'+platform+'_min=','run_'+platform+'_min='+job[1])
    sed_inplace(master_control_file,'run_gpu_MD=','run_gpu_MD='+job[2])
    sed_inplace(master_control_file,'run_cpu_rec=','run_cpu_rec='+job[3])
    sed_inplace(master_control_file,'run_energy_analysis=','run_energy_analysis='+job[4])
    sed_inplace(master_control_file,'inpdir=','inpdir='+inp_dir)
    
    print("md dir ", md_dir, "\ninp_dir", inp_dir, JOB_FILE)
    # learn to use subprocess Popen etc to do this later
    sed_inplace(JOB_FILE,'THE_DIR=','THE_DIR='+md_dir)
    sed_inplace(JOB_FILE,'HQDIR=','HQDIR='+md_dir)
    sed_inplace(JOB_FILE,'DONOR=PI','DONOR='+MD_Donor)
    sed_inplace(JOB_FILE,'STATE=PI','STATE='+MD_dyn_state)
    if len(Solvent.split(",")) == 1:
        sed_inplace(JOB_FILE,'SOLVENT=','SOLVENT='+Solvent)
    else:
        sed_inplace(JOB_FILE,'d_resname=','d_resname='+Solvent.split(',')[0])
        sed_inplace(JOB_FILE,'a_resname=','a_resname='+Solvent.split(',')[1])
        sed_inplace(JOB_FILE,'SOLVENT=','SOLVENT='+ Solvent.split(',')[0] + Solvent.split(',')[1]) 
    sed_inplace(JOB_FILE,'SYSTEM=','SYSTEM='+system)
    sed_inplace(JOB_FILE,'MOLECULE=',"MOLECULE="+system+testing_string)
    sed_inplace(JOB_FILE,"ID_NO=",   "ID_NO="+testing_string)
    sed_inplace(JOB_FILE,'state_list=','state_list='+state_list)
    # sed_inplace(JOB_FILE, "THE_DIR=")

    sed_inplace(master_control_file, 'state_list=','state_list='+state_list)

    sed_inplace(master_control_file,'DEP=','DEP='+dependency)
    #else:
    data2 = os.popen('cat '+JOB_FILE).read()
        
    os.chdir(md_dir)
    data = os.popen('sh master_control.sh').read()
    
    
    
    return data

def generate_MDinp(MD_inp_tempdir, MD_dt, MD_steps, sample_steps, MD_temperature, MD_running_dir) :
    construct_MDinp = shell_copy_template(MD_inp_tempdir, MD_running_dir, filename='construct_MDinp')
     
    os.chdir(MD_running_dir)
    sed_inplace(construct_MDinp,'Template_Dir=',  'Template_Dir=' + MD_inp_tempdir )
    sed_inplace(construct_MDinp,'Target_Dir=',    'Target_Dir=' + MD_running_dir )
    
    sed_inplace(construct_MDinp,'MD_steps=',       'MD_steps='+ MD_steps )
    sed_inplace(construct_MDinp,'MD_dt=',          'MD_dt='+ MD_dt )
    sed_inplace(construct_MDinp,'sample_steps=',   'sample_steps='+ sample_steps )
    sed_inplace(construct_MDinp,'MD_temperature=', 'MD_temperature='+ MD_temperature )
    
    
    data = os.popen('sh construct_MDinp.sh; pwd').read()
    return data 


    
import pickle 
def main_RunMD():
    dict_of_simulation =  parse_input()
    # prepare for quantum chemistry simulation 
    work_dir = dict_of_simulation['work_dir']
    # where we have the structure directory
    structure_dir = dict_of_simulation['structure_dir']
    # load the project folder 
    project = dict_of_simulation['project']
    # get case id for the simulation
    caseid = dict_of_simulation['caseid']
    # get job list
    case_id_list = separate_idlist(caseid)
    # get list of job control token 
    job_control  = init_job_control(case_id_list)

    # job_info = init_job_info(case_id_list)
    startime = time.time ()

    # Read list to memory
    def read_list():
        # for reading also binary mode is important
        with open(work_dir+'/listfile', 'rb') as fp:
            n_list = pickle.load(fp)
            return n_list 
        
    job_info = read_list()
    print("RunMD begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step4(ind, job_control, dict_of_simulation,job_info)

    end_time = time.time  ()
    print("RunMD ends at "+str(end_time))


if __name__ == "__main__":
    main_RunMD(); 
