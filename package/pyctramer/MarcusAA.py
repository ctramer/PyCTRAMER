import numpy as np
import pandas as pd
from .Utilities import *

#MARCUS All ATOM 
def run_step5(i, job_control, dict_of_simulation,job_info):
    #print(arrow_str,"Started Step 5: MarcusAA")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] # caseid = str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    MD_dyn_state = dict_of_simulation['MD_traj_state']
    if dict_of_simulation['skip_MD'] == 'true' or dict_of_simulation['skip_MD'] == 'True':
        enelst = dict_of_simulation['Ene_list'].split(',')
        statlst = dict_of_simulation['State_list'].split(',')
        cplinlst = dict_of_simulation['coupling_list'].split(",")
        stindlst = dict_of_simulation['State_order'].split(",")
        Ene_list = []
        State_list = {}
        coupling_list = []
        for ind, energ in enumerate(enelst) :
            Ene_list.append(float(enelst[ind]))
            State_list[statlst[ind]] = int(stindlst[ind])
        for cpppp in cplinlst :
            coupling_list.append([float(cpppp)])
        job_info[i] ['Ene_list']  = Ene_list
        job_info[i]['State_list'] = State_list 
        job_info[i]['coupling_list'] = coupling_list
    else:
        Ene_list = job_info[i]['Ene_list'] 
        State_list = job_info[i]['State_list'] 
        coupling_list = job_info[i]['coupling_list'] 
    
    system_name = project +"_"+ caseid
    MD_temperature  = dict_of_simulation["MD_temperature"] 
    temperature = float(MD_temperature)
    charge_list = []
    for info in job_info:
        for state_name in info['State_list']:
            charge_list.append(state_name)
    # charge_list = ['PI','CT1','CT2','GR']
    # job_info[i][]

    print("MD directory: ", md_dir)
    # Trajectory dir 
    traj_dir = md_dir + '/Traj_' + MD_dyn_state
    print("MD trajectory directory: ",traj_dir)
    # get correction for Er
    Correction =  get_energy_correction(traj_dir,charge_list,Ene_list)
    print("Energy corrections (eV): ", Correction)
    
    Coupling_list, rates, transitions_names, Er_all, dE_all = returnMarcusRate(traj_dir, MD_dyn_state,coupling_list, Correction, charge_list,temperature, system_name, ensemble='NVE')
    # print("Marcus rates (Hz) : ", rates)
    print(transitions_names)
    print(arrow_str, 'Summary of Marcus rate constants from all-atom simulation:')
    for indtra, tra in enumerate(transitions_names):
        print("Marcus rate for transition " , tra ," is ", rates[indtra], ' Hz')
    f = open(case_dir + '/'+project+'_'+caseid+'_MarcusAA.out', 'a+')
    #result = {'GR': 0, 'PI': pi_state, 'CT1':CT1_state, 'CT2': CT2_state}
    #energies =[0, PI_energy,CT1_energy, CT2_energy] # list of energy in eV
    for ene in Ene_list: 
        f.write('Excited state order  '+ str(ene)+'\n')
    f.write("\n")
     
    for ind,state in enumerate(charge_list):
        f.write(state+' '+str(Ene_list[ind])+'\n')
    
    f.write("\n") 

    f.write('Excited state energies in eV  \n')
    for ind, name in enumerate(charge_list):
        f.write("Excited energy of the state "+name+' is '+ str(Ene_list[ind]) + ' eV. \n' )
    f.write("\n")

    f.write('diabatic coupling in eV  ')
    for ind, coupling in enumerate(Coupling_list):
        f.write ('Coupling between transition of '+transitions_names[ind] + ' is ' + str(coupling) + ' eV. \n ')
    
    f.write('Rates in Hz  ')
    for ind, rate in enumerate(rates):
        f.write ('Marcus rates at ' + MD_temperature + ' between transition of '+transitions_names[ind] + ' is ' + str(rate/1e12) + ' THz. \n ')
    f.close()

    job_info[i]['transitions_names'] = transitions_names
    job_info[i]['transitions_rates'] = rates
    job_info[i]['Correction'] = Correction
    job_info[i]['Coupling_list'] = Coupling_list

    indmax = rates.index(max(rates))
    namemax = transitions_names[indmax]
    namemax = namemax.split("-")
    cpmax = Coupling_list[indmax]
    correctmax = -Correction['CT'+str(indmax+1)] + Correction['LE']
    ratemax = rates[indmax]
    max_indpair = [1,2+indmax]
    
    job_info[i]['namemax']  = namemax
    job_info[i]['cpmax']    = cpmax
    job_info[i]['ratemax']  = ratemax 
    job_info[i]['correctmax']  = correctmax
    job_info[i]['max_indpair']  = max_indpair

    H_sys = np.array([[1.0,0],[0,1.0]])
    H_sys[0][1] = H_sys[1][0] = cpmax * eV2au
    H_sys[0][0] = 0. 
    H_sys[1][1] = dE_all[indmax] * eV2au 
    print("H_sys (Hartree) \n", H_sys)
    np.savetxt(case_dir+'/H_sys.dat' ,H_sys)

    # print("Correction ", Correction)
    # print(Coupling_list,'Coupling list')

    stroutput = "\n" + arrow_str
    stroutput = stroutput + "STEP 5: MarcusAA  \n"
    # stroutput 
    stroutput = stroutput + arrow_str + "Marcus Rates Summary \n" 
    stroutput = stroutput + "transition,Excitation Energy (eV),Correction (eV),Coupling (eV),Er (eV),DeltaE (eV),all-atom Marcus rates (Hz)\n"
    for ind, name in enumerate(State_list):
        if ind > 1:
            stroutput = stroutput + name + "," + str(Ene_list[ind]) +"," 
            stroutput = stroutput +  str(Correction[name]) + ","  
            stroutput = stroutput + str(Coupling_list[ind-2])  + "," 
            stroutput = stroutput + str(Er_all[ind-2]) + "," 
            stroutput = stroutput + str(dE_all[ind-2]) + "," 
            stroutput = stroutput + str(rates[ind-2]) + "\n" 
    stroutput = stroutput + "\n"

    stroutput = stroutput + arrow_str + "Target Transition \n" 
    stroutput = stroutput + "Transition between: " + transitions_names[indmax] + "\n" 
    stroutput = stroutput + "System Hamiltonian (Hartree) \n "
    for ii in range(len(H_sys)):
        for j in range( len(H_sys)):
            stroutput = stroutput + str(H_sys[ii][j]) + " "
        stroutput = stroutput + "\n"

    
    stroutput = stroutput + "Rate: " + str(ratemax) + " Hz \n" 

    stroutput = stroutput + "Energy correction of the transition: " + str(correctmax) + " eV \n" 
    stroutput = stroutput + "DeltaE = " + str(dE_all[indmax]) + " eV \n" 
    stroutput = stroutput + "LE excition energy = " + str(Ene_list[indmax]-1) + " eV \n" 
    stroutput = stroutput + "Diabatic coupling of the transition = " + str(cpmax) + "eV \n" 
    stroutput = stroutput + "Er = " + str(Er_all[indmax]) + " eV \n" 
     
    #stroutput =  stroutput + '\n\n'
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 
    return job_control

def get_energy_correction(filedir,charge_list,ex_energy_list):
    """
    This functions extract energy correction for replacting solute energy obtained from MD 
    with the value obtained from quantum chemistry calculation
    
    filename for MD energy in gas is assumed to be energy_gas_[STATE]_TRAJ_NVT.dat in kcal/mol
    reported energy_correction dictionary is in eV.
    """
    kcalmol2eV= 0.0433634 #0.0433634
    energy_corr_tmp = []
    energy_corr = {} # Uxg
    final_energy_corr = {}
    
    for ind, i in enumerate (charge_list):
        filename = f'{filedir}/energy_gas_' + i + '_TRAJ_NVT.dat'
        test_file = pd.read_table(filename,header=None)
        energy_in_ev = test_file[0][0] * kcalmol2eV
        energy_corr_tmp.append(energy_in_ev)
        energy_corr[i] = energy_in_ev -    energy_corr_tmp[0]  
        final_energy_corr[i] = ex_energy_list[ind] - energy_corr[i]

    # getting energy gap between excited state and ground (energy_corr_tmp[0])  as obtained from MD  

    print("MD solute energies (eV): ", energy_corr)
    print("QM solute energies (eV): ", ex_energy_list)
    # substracting the excitation energy obtained using MD from the energy calculated by QChem 
    
    return final_energy_corr  

def calc_k_marcus(DeltaE, Er, Gamma, T):
    """calculate charge transfer rate constant based on marcus theory.
    
    Parameters:
    ------------
    Ea : 
    Er :
    Gamma : 
    
    return:
    ------------
    k_M
    
    """
    h = 1 * 4.135667696e-15    # h value in eV s
    # note that hbar = h * 2 * PI
    pi = np.pi 
    kb = 8.617333262e-5
    # T = 300 # in kelvin
    
    #h = 4.135667696...×10−15 eV⋅Hz−1
    term1 = ( 2 * pi * abs(Gamma)**2 / h ) 
    #term1  = abs(Gamma)**2/h
    term2 = np.sqrt(pi /(kb * T * Er) ) 
    term3 = np.exp(- (DeltaE + Er)**2 / (4 * kb * T * Er))
    # print(term1,term2,term1*term2,term3)
    k_M = term1 * term2 * term3
                   
    return k_M

def returnMarcusRate(traj_dir, MD_traj_state, Coupling_list, correction, state_list, temperature,system_name, ensemble='NVE'):
    Ene_list = []
    for name in correction:
        Ene = np.loadtxt(traj_dir+'/'+system_name+'_E_'+name+"_TRAJ_"+MD_traj_state+'_'+ensemble+'.dat')
        Ene *= 0.0433634 # kcal per mol to eV
        Ene_list .append(Ene.copy() + correction[name])
        np.savetxt(traj_dir+'/'+system_name+'_Ecorrected_'+name+"_TRAJ_"+MD_traj_state+'_'+ensemble+'.dat', Ene_list[-1])
        # print(len(Ene_list))
        # print(name, Ene_list)
    kb = 8.617333262e-5
    T = temperature # 300 Kelvin 
    rates_list = []

    transition_pair = []
    transition_indx = []
    for ind, name in enumerate(state_list):
        if ind > 1:
            transition = state_list[1] + "-" + name
            transition_pair.append(transition)
            transition_indx.append([1,ind])

    coupling_list = []
    coup_util_list = []
    cp_count = 0
    # print(Coupling_list, "cplst")
    for i in range(len(state_list)):
        for j in range(len(state_list)):
            
            if i > j:
                #print(i,j)
                #print(cp_count)
                if j == 1:
                    coupling_list.append(Coupling_list[cp_count][0]) # cp_count)
                    coup_util_list.append(cp_count)
                    # print(Coupling_list[cp_count][0])
                cp_count +=1 


    #print("Coupling_list", Coupling_list,'\n  couplist ',  coupling_list)
    inds = transition_indx
    #print(len(inds),'inds')
    #print((2 * kb * T))
    Er_list = []
    dE_list = []
    for ind, coupling in enumerate(coupling_list):
        energy_gap = Ene_list[inds[ind][0] ] -  Ene_list[inds[ind][1]]
        U = np.average(energy_gap) # <U>
        sigma_U = np.std(energy_gap)
        Er = np.std(energy_gap) * np.std(energy_gap) / (2 * kb * T)
        dE = - U - Er 
        rates_list.append(calc_k_marcus(dE,Er,coupling_list[ind], T) )
        Er_list . append(Er) 
        dE_list . append(dE)
        print(arrow_str, 'Transition', transition_pair[ind] ,'Marcus Parameters from all-atom simulation' )
        print("DeltaE = ", dE , " eV" )
        print("<U> = ", U , " eV" )
        print("Er = ", Er , " eV" )
        print("Gamma = ", coupling_list[ind],' eV')
        print("sigma_U = ", sigma_U, ' eV')
        # test Ea
        print("Ea = ", 0.5 * U * U / sigma_U / sigma_U * kb * T , (dE+Er)**2/4./Er, ' eV')
        print("Marcus rate = ", rates_list[ind], ' Hz\n')
        #print("DeltaE = ",dE, " eV, Er = ",Er, " eV, diabatic coupling = ", coupling_list[ind], ' eV')
        #print(rates_list[ind], ' Hz is rate of the transition ' + transition_pair[ind])
    return coupling_list, rates_list, transition_pair, Er_list, dE_list


  
import pickle 
def main_MarcusAA():
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

    def write_list(a_list):
        # store list in binary file so 'wb' mode
        with open(work_dir+'/listfile', 'wb') as fp:
            pickle.dump(a_list, fp)
            print('Done writing list into a binary file')
        
    job_info = read_list()
    print("MarcusAA begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step5(ind, job_control, dict_of_simulation,job_info)

    write_list(job_info)
    end_time = time.time  ()
    print("MarcusAA ends at "+str(end_time))


if __name__ == "__main__":
    main_MarcusAA()
