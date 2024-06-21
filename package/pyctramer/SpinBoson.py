import numpy as np
import pandas as pd
from .Utilities import *


def run_step6(i,job_control,dict_of_simulation,job_info):
    # print('')
    #print(arrow_str,"Started Step 6: SpinBoson")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] # caseid = str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    MD_dyn_state = dict_of_simulation['MD_traj_state']
    spin_boson_modes = dict_of_simulation['spin_boson_modes']
    TCF_corr_length = dict_of_simulation['TCF_corr_length']
    TCF_DT = dict_of_simulation['TCF_DT']
    N_traj = int(dict_of_simulation['N_traj'])
    MD_steps = int(dict_of_simulation['MD_steps'])
    sample_steps = int(dict_of_simulation['sample_steps'])
    traj_N_frames = int(MD_steps / sample_steps)
    temperature  = float(dict_of_simulation['MD_temperature'])
    beta = 1. / (temperature * kT2au)

    TCF_corr_length = float(TCF_corr_length)
    spin_boson_modes = int(spin_boson_modes)
    TCF_DT = float(TCF_DT)
    corr_leng = int(TCF_corr_length/TCF_DT)
    # print('traj_N_frame',traj_N_frames,corr_leng,TCF_corr_length,)
    correctmax = 0
    if job_info[i]['namemax'] == '':
        energy_traj_file =  dict_of_simulation['energy_traj_file'].split(",") 
        namemax = dict_of_simulation['transition_DA_states'].split(",")
        correctmax = float(dict_of_simulation['correction'])
        LE_path = energy_traj_file[0]
        CT_path = energy_traj_file[1]
    else:
        namemax = job_info[i]['namemax'] 
        cpmax = job_info[i]['cpmax']   
        ratemax = job_info[i]['ratemax']  
        correctmax = job_info[i]['correctmax'] 
        max_indpair = job_info[i]['max_indpair'] 
        LE_path = md_dir + "/Traj_"+MD_dyn_state+"/"+project +"_"+ caseid+"_E_LE_TRAJ_"+MD_dyn_state+"_NVE.dat"
        CT_path = md_dir + "/Traj_"+MD_dyn_state+"/"+project +"_"+ caseid+"_E_"+namemax[1]+"_TRAJ_"+MD_dyn_state+"_NVE.dat"
    kcal2eV = 0.0433641043
    print("Transition of the spin-boson model is ", namemax[0], namemax[1])
    print('Energy gap correction of the transition (eV): ', correctmax) #, traj_N_frames)
    Ut = kcal2eV * (np.loadtxt(LE_path) - np.loadtxt(CT_path)) + correctmax
    
    # Get Cuu(t) 
    Cuu = gen_Cuu(Ut, corr_leng, N_traj, traj_N_frames)
    print('C_UU(t=0) = ', Cuu[0])
    np.savetxt(work_dir+"/Cuu"+namemax[0]+namemax[1]+".dat",Cuu)

    # Get discretized SD
    w_SB = np.zeros(spin_boson_modes, dtype=float)
    max_iter_SBM = 1000 
    tol_w_SBM = 1e-8
    for im in range(spin_boson_modes):
        x0 = 0.1 * cm2au 
        x1 = 0.5 * cm2au
        iter_count = 0
        while (abs(x1-x0) > tol_w_SBM):
            f0 = fw(x0, im+1, TCF_DT*ps2au, Cuu, spin_boson_modes) 
            f1 = fw(x1, im+1, TCF_DT*ps2au, Cuu, spin_boson_modes)
            x = (x0 * f1 - x1 * f0) / (f1 - f0) 
            x0 = x1 
            x1 = x 
            iter_count += 1
            if iter_count > max_iter_SBM:
                print('Warning: The frequency of mode '+ str(im+1)+ ' is not converged yet! ')
                # raise RuntimeError("set larger mode number would help")
        w_SB[im] = x1 # get discretized spectral density in au of tiem
    print('frequencies of the bath (au) \n', w_SB) 
    Req_SB = np.zeros(spin_boson_modes)
    c_SB   = np.zeros(spin_boson_modes)
    Er = 0.5 * beta * Cuu[0] * eV2au  * eV2au 
    print('Er = ', Er ,'Hartree')
    for im in range(spin_boson_modes):
        Req_SB[im] = np.sqrt(2. * Er / (1.*spin_boson_modes) ) / w_SB[im]
        c_SB[im] = np.sqrt(Er / (2. * spin_boson_modes)) * w_SB[im]

    SB_data = np.array([np.arange(1,spin_boson_modes+1,dtype=int), w_SB, c_SB, Req_SB, Req_SB]).T
    np.savetxt(work_dir+"/SBM"+namemax[0]+namemax[1]+".dat",SB_data)
    np.savetxt(work_dir+"/SBM"+namemax[0]+"CT.dat",SB_data)
    # Get J(w)
    # max frequency to be 
    max_freq = 1.2 * max(w_SB) * au2cm # 2500 # 2 * np.pi / (TCF_DT * ps2au) / corr_leng
    print('max freq of the bath in cm-1', max_freq)
    dw = 0.2 # in wave number 
    max_w_count = int(max_freq / dw)
    
    cont_w = 0.2 * cm2au * np.arange(max_w_count)
    Jw_cont = np.zeros(max_w_count)

    for ind, w in enumerate(cont_w):
        Jw_cont[ind] = Jw(w, beta, TCF_DT*ps2au, Cuu)
    
    np.savetxt(work_dir+"/Jw"+namemax[0]+namemax[1]+".dat",Jw_cont*eV2au*eV2au)

    stroutput = "\n" 
    stroutput = stroutput + arrow_str  +  "STEP 6: SpinBoson\n"

    stroutput = stroutput + "Cuu data point number " + str(corr_leng) + ",  dt  = " + str(TCF_DT) + " ps\n" 
    stroutput = stroutput + "Jw data point number " + str(max_w_count) + ", dw  = " + str(dw) + "cm^-1 \n"
    
    stroutput = stroutput + "SpinBoson model bath parameters " + work_dir+"/SBM"+namemax[0]+"CT.dat\n"  
    stroutput = stroutput + "SpinBoson bath parameter file headers: w_SB, c_SB, Req_SB, Req_SB \n" 
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 
    # print(job_control, 'step 3 ending')
    return job_control

def gen_Cuu(energy_gaps, corr_leng, Ntraj, MD_steps):
    # energy_gaps = np.array(energy_gaps)
    Cuut = np.zeros(corr_leng,dtype=float)
    Cuu0 =  np.average(energy_gaps)**2
    energy_gaps = energy_gaps - np.average(energy_gaps)
    print("Cuu0",Cuu0,np.std(energy_gaps)**2)
    for i in range(Ntraj):
        for j in range(corr_leng):
            for k in range(MD_steps-j):
                Cuut[j] += energy_gaps[i*MD_steps+k] * energy_gaps[i*MD_steps+k+j]  
            Cuut[j] /= (1.*MD_steps-1.*j-1.)
            #print((1.*MD_steps-1.*j))
    # print('cuu',Cuut[0],Cuu0 ) 
    for j in range(corr_leng):
        Cuut[j] /= (1. * Ntraj) 
        #Cuut[j] = Cuut[j] - Cuu0
    # print(Cuut[0] ,Cuu0)
    return Cuut



def Jw(w, beta, dt_au, Cuu):
    factor = 0.25 * beta * w 
    size = len(Cuu)
    # print('size',size, np.shape(w))
    ft = Cuu * np.cos(w * np.arange(size) * dt_au) # wt 
    return factor * (np.sum(ft)-0.5*(ft[0]+ft[-1])) * dt_au


def fw(w, j, dt, Cuu, N):
    # w in au 
    # j mode order 
    # dt in au 
    # Cuu in au 
    # N: SBM mode num
    factor = 2 * N * w / np.pi / Cuu[0] 
    size =  len(Cuu) 
    ft = np.zeros(size)
    ft[0] = Cuu[0]
    for s in range(1,size):
        wt = s * dt * w 
        ft[s] = Cuu[s] * np.sin(wt) / wt 
    return factor * (np.sum(ft) - 0.5* (ft[0]+ft[-1]))*dt + 0.5 -j 

import pickle 
def main_SB():
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
    print("SpinBoson begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step6(ind, job_control, dict_of_simulation,job_info)

    end_time = time.time  ()
    print("SpinBoson ends at "+str(end_time))


if __name__ == "__main__":
    main_SB()
