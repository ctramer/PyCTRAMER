from .Utilities import *  
import pandas as pd
import numpy  as np
import math, time
import io 
def run_step3(i, job_control, dict_of_simulation, job_info):
    #print(arrow_str,"Started Step 3: ConstructFF")
    project  = dict_of_simulation['project']
    work_dir = dict_of_simulation['work_dir']
    caseid =  job_control[i]['caseid'] # caseid = str(i)
    case_dir = work_dir + '/' + project +"_"+ caseid
    qc_dir,md_dir =  get_folder_path(case_dir)
    template_dir = dict_of_simulation['template_dir']

    N_cis = int(dict_of_simulation['N_cis'])

    # charges on donor
    charge_on_frags = get_charge_on_Donor(qc_dir, project +"_"+ caseid, N_cis)
    charges_on_donor = []
    for line in charge_on_frags:
        splitted_charges = line.split()
        charges_on_donor.append(float(splitted_charges[1]))
    print("fragment charge on donor for ", N_cis, " states: \n"  ,charges_on_donor)
    
    
    # find locally excited states
    LE_fragment = dict_of_simulation['LE_fragment'] 
    DA_type = dict_of_simulation['DA_type']
    CT_donor_fragments = dict_of_simulation['CT_donor_fragments']
    CT_acceptor_fragments = dict_of_simulation['CT_acceptor_fragments']
    Max_CT_states = dict_of_simulation['Max_CT_states']
    CT_amount_threshold = dict_of_simulation['CT_amount_threshold']
    OS_threshold = dict_of_simulation['OS_threshold']
    LE_position_threshold = dict_of_simulation['LE_position_threshold']
    CT_type  = dict_of_simulation['CT_type']
    Include_bCT_states = dict_of_simulation['Include_bCT_states']
    
    
    # Load the theodore results
    test_results = load_Theo_summary(qc_dir, DA_type)
    # extract energy and excited state
    if DA_type == 'DBA' and Max_CT_states == "2":
        print(arrow_str, "Charge on Donor Fragments")
        CT_deltaQ_onD = []
        CT_Donor_charge = []
        State_list, Ene_list = Find_CT_states_DBA(test_results, LE_fragment,CT1_thres=CT_amount_threshold,OSthres=OS_threshold)
        for ind, namestate in enumerate(State_list):
            CT_Donor_charge.append(charges_on_donor[State_list[namestate]-1])
            CT_deltaQ_onD.append(charges_on_donor[State_list[namestate]-1] - charges_on_donor[State_list['LE']-1])

    elif DA_type == 'DA':
        # , CT_dQ_D, CT_donor_charge, bdCT_list
        State_list, Ene_list, CT_deltaQ_onD, CT_Donor_charge, bdCT_list = Find_CT_states_DA(test_results, LE_fragment, charges_on_donor, Max_CT_states, Include_bCT_states, POSthres=LE_position_threshold, OSthres=OS_threshold, CT_thres=CT_amount_threshold)
    
    # get total number of atoms
    N_atom = dict_of_simulation['N_atom']
    # get the output name
    out_name = project + '_' + caseid
    # coupling for different states
    coupling_list = []
    N_cis     =  dict_of_simulation['N_cis']
    Solvent_dir =  dict_of_simulation['solvent_dir']
        
    
    # calculate the diabatic coupling
    for ind, name in enumerate(State_list):
        for ind2, name2 in enumerate(State_list):
            if ind > ind2 :
                #print(name, name2)
                coupling_list.append(get_coupling(qc_dir,out_name,State_list[name],State_list[name2],N_cis))
    print("State list = ",State_list)#,"\n")# Diabatic coupling (eV) = ",coupling_list) 
    print("States",end='')
    for ind1, name in enumerate(State_list):
        print(","+name, end='')
    
    fcdcpind = 0
    for ind1, name in enumerate(State_list):
        # stroutput = stroutput + name + 
        if ind1 > 0:
            print(name, end='')
        for ind2, name2 in enumerate(State_list):
            if ind1 > ind2 :
                if ind1 > 0 and ind2 > 0:
                    print(","+str(coupling_list[fcdcpind]),end='')
                    # stroutput = stroutput + name + name2 + ":" + str(coupling_list[fcdcpind]) + "\n"
                fcdcpind += 1
    # get charge 
    charge_name = dict_of_simulation['charge_type'] #'Mulliken'
    # Solvent =  dict_of_simulation['Solvent']
    # Solvent_ad = dict_of_simulation['solvent_dir']
    
    # Solvent_name = dict_of_simulation['solvent_name']
    # Solute_name = dict_of_simulation['solute_name']
    # Solute_name1 = 'TRI'
    
    BoxLength=dict_of_simulation['MD_box_length']
    COM = int(BoxLength)/2
    
    DA_type  = dict_of_simulation['DA_type']
    D_range  =  dict_of_simulation['D_range']
    B_range  =  dict_of_simulation['B_range']
    A_range  =  dict_of_simulation['A_range']
    A1_range = dict_of_simulation['A1_range']
    A2_range = dict_of_simulation['A2_range']
    CT_donor_fragment = dict_of_simulation['CT_donor_fragments']
    CT_acceptor_fragment = dict_of_simulation['CT_acceptor_fragments'] 

    solute_resname = dict_of_simulation['solute_resname']
    solvent_resname = dict_of_simulation['solvent_resname']
    solute_N_atom  = dict_of_simulation['solute_N_atom']
    solvent_N_atom  = dict_of_simulation['solvent_N_atom']
    solvent_dir      = dict_of_simulation['solvent_dir']
    solvent_N_mole   = dict_of_simulation['solvent_N_mole']
    solute_dir       = dict_of_simulation['solute_dir']
     

    # THF_number= dict_of_simulation['solvent_nres']  
    N_cis   = int(dict_of_simulation['N_cis'])
    N_atom = int(dict_of_simulation['N_atom'])  
    force_field_name='gaff'
    charge_file=  qc_dir+ '/charge0.txt'
    
    zeros = np.zeros(N_atom, dtype=int) 
    np.savetxt(charge_file, zeros)
    
    Solvent_name = ""
    if len(solvent_resname.split()) == 2:
        Solvent_name = solvent_resname.split()[0] + solvent_resname.split()[1]
    elif  len(solvent_resname.split()) == 1:
        Solvent_name = solvent_resname

    if CT_type == "inter":
        info = construct_prmtop_inter(template_dir, qc_dir, caseid, project, 
                                      D_range,A_range,B_range,CT_acceptor_fragment,CT_acceptor_fragment,
                                      solute_dir,solvent_dir,solute_resname,solvent_resname,solvent_N_mole,
                                      BoxLength, COM, force_field_name, N_cis, charge_file)
                                    
    else: 
        info = construct_prmtop_intra(template_dir, qc_dir, caseid, project, 
                                      D_range,A_range,B_range,CT_acceptor_fragment,CT_acceptor_fragment,
                                      solute_dir,solvent_dir,solute_resname,solvent_resname,solvent_N_mole,
                                      BoxLength, COM, force_field_name, N_cis, charge_file) 
    
    ## construct FF excited states
    charge_file_prefix = project+'_'+caseid 
    mulliken_log_file_dir = qc_dir +'/mulliken'
    
    Solvent_name = ""
    if len(solvent_resname.split(",")) == 2:
        Solvent_name = solvent_resname.split(',')[0] + solvent_resname.split(',')[1]
    elif  len(solvent_resname.split(',')) == 1:
        Solvent_name = solvent_resname

    old_top = qc_dir + charge_file_prefix + '_' + Solvent_name +'.prmtop' 
    charge_ad = qc_dir  

    new_name_top = charge_file_prefix+'_' +Solvent_name+'_GR'
    new_top_dir = md_dir
    print("State list = ",State_list)
    QChem_version = dict_of_simulation['QChem_version']
    print("N_atom = ", N_atom)
    #replace_charge_to_prmtop(N_atom_triad, new_top_dir, old_top, charge_file, new_name_top)
    for ind, name in enumerate(State_list):
        charge_file = qc_dir+'/charge_'+name+'_'+charge_file_prefix+'.txt'
        get_charges(N_atom, State_list[name], name, mulliken_log_file_dir, charge_ad, 
                    charge_file_prefix,charge_file_prefix, charge_name, D_range, B_range, A_range, QChem_version)
        new_name_top = charge_file_prefix+'_'+Solvent_name+'_'+name
        replace_charge_to_prmtop(N_atom, new_top_dir, old_top, charge_file, new_name_top)
    print("Energies (eV) = ", Ene_list)
    # print('jobinfo=',job_info)

    job_info[i]['Ene_list']= Ene_list 
    job_info[i]['State_list'] = State_list
    job_info[i]['coupling_list'] = coupling_list
    
    stroutput = "\n" + arrow_str
    stroutput = stroutput + "STEP 3: ConstructFF\n"
    stroutput = stroutput + "Input: \n" + qc_dir 
    stroutput = stroutput + "\n" + out_name + ".out\n"
    stroutput = stroutput + out_name + ".fchk\n"
    stroutput = stroutput + qc_dir + "/mulliken/" + out_name + ".out\n"
    stroutput = stroutput + qc_dir + "tden_summ_fchk.txt\n"
    stroutput = stroutput + qc_dir + "tden_summ_tddft.txt\n"
    stroutput = stroutput + qc_dir + "packmol.inp\n"
    stroutput = stroutput + qc_dir + "tleap.in\n\n"
    stroutput = stroutput + "Job submit script: \n" + qc_dir + "construct_prmtop_"+ CT_type +".sh\n\n"
    stroutput = stroutput + "Output: \n" + qc_dir + "\n"
    stroutput = stroutput + qc_dir + "leap.log"
    stroutput = stroutput + qc_dir + "FCD_coupling_raw.txt"
    for ind, name in enumerate(State_list):
        charge_file = qc_dir+'/charge_'+name+'_'+charge_file_prefix+'.txt'
        stroutput = stroutput + charge_file + "\n"
    
    stroutput = stroutput + old_top
    for ind, name in enumerate(State_list):
        new_name_top = charge_file_prefix+'_'+Solvent_name+'_'+name
        stroutput = stroutput + new_name_top + ".txt\n"
    stroutput = stroutput + qc_dir + charge_file_prefix + '_' + Solvent_name +'.inpcrd \n' 
    stroutput = stroutput + "\n"
    stroutput = stroutput + arrow_str  + "Summary of calculation:\n\n"
    stroutput = stroutput + arrow_str + "State characteristic analysis (TheoDORE)\n"
    s = io.StringIO() 
    test_results.to_csv(s)
    stroutput = stroutput + s.getvalue()
    stroutput = stroutput + "\n"

    stroutput = stroutput + arrow_str + "States considered in Marcus rate calculation: \n"
    
    stroutput = stroutput + "State, Excitation Energy (eV) \n"
    e  = np.array(Ene_list)
    for ind, name in enumerate(State_list):
        stroutput = stroutput + name + "," + str(e[ind]) + '\n'
    
    stroutput = stroutput + arrow_str  + "CT States: Charge Transfer Amounts \n"
    
    if CT_type == "DA":
        stroutput = stroutput +   "State,dQ_D (eV),Q_D \n"
        for ind, name in (State_list):
            if ind > 1: # CT_deltaQ_onD, CT_Donor_charge, bdCT_list
                stroutput = stroutput + name + "," + str(CT_deltaQ_onD[ind-2]) +"," +str(CT_Donor_charge[ind-2]) + "\n"
        
        stroutput = stroutput + "\n" 
        stroutput = stroutput + arrow_str + "CT States: Optical Properties \n"
        for line in bdCT_list:
            stroutput = stroutput + "\n"

    elif CT_type == "DBA":
        stroutput = stroutput +   "State,dQ_D (eV),Q_D \n"
        for ind, name in (State_list):
            if ind > 1: # CT_deltaQ_onD, CT_Donor_charge, bdCT_list
                stroutput = stroutput + name + "," + str(CT_deltaQ_onD[ind-2]) +"," +str(CT_Donor_charge[ind-2]) + "\n"
        
        stroutput = stroutput + "\n" 
        
        
    stroutput = stroutput + "\n" + arrow_str  + "FCD Coupling in eV \n"    

    fcdcpind = 0
    for ind1, name in enumerate(State_list):
        # stroutput = stroutput + name + 
        for ind2, name2 in enumerate(State_list):
            if ind1 > ind2 :
                if ind1 > 0 and ind2 > 0:
                    stroutput = stroutput + name + name2 + ":" + str(coupling_list[fcdcpind]) + "\n"
                fcdcpind += 1
    # stroutput = stroutput + '\n\n'
    stroutput = stroutput + "#######################################################\n"
    write_to_output(i, job_control, dict_of_simulation, stroutput) 
    return job_control

def load_Theo_summary(filedir, DA_type):
    """
    load theodore analysis results 
    """
    file_tddft=f'{filedir}/tden_summ_tddft.txt'
    file_fchk=f'{filedir}/tden_summ_fchk.txt'
    file_ehFrag=f'{filedir}/ehFrag.txt'
    theo_tddft=pd.read_table(file_tddft, skiprows= [1],delim_whitespace=True )
    theo_fchk=pd.read_table(file_fchk, skiprows= [1], delim_whitespace=True )
    theo_ehFrag=pd.read_table(file_ehFrag, skiprows= [1], delim_whitespace=True )
    if DA_type == 'DA':
        theo_summ = pd.concat([theo_tddft[['state','dE(eV)','f']],
                           theo_fchk[['Om','POS','PR','CT','COH', 'CTnt', 'PRNTO','Z_HE','RMSeh']],
                           theo_ehFrag[['H_1','E_1','H_2','E_2']]],axis=1)
    elif DA_type == "DBA":
        theo_summ = pd.concat([theo_tddft[['state','dE(eV)','f']],
                            theo_fchk[['Om','POS','PR','CT','COH', 'CTnt', 'PRNTO','Z_HE','RMSeh']],
                            theo_ehFrag[['H_1','E_1','H_2','E_2','H_3','E_3']]],axis=1)
    return theo_summ

def Find_CT_states_DBA(data, LE_fragment, POSthres=0.2, OSthres= 0.1, CT1_thres=0.5, CT2_thres=0.05, verbose=True):
    """
    Algorithm to find the best LE, CT1 and CT2 states
    LE_fragment defines the locally excited fragments in triad or dyad; int; for triad case; LE_fragment=2 for porphyrin excited state
    
    OSthres: threshold of the Oscillator Strength (OS)
    POSthres: thereshold of the POS (averaged position)
    CT1_thres: threshold of the charge transfer from C to PC60
    
    return 
    list of CT states final decision, and the corresponding energies 
    """
    
    ################ Finding LE State ###########################
    energy_list= list(data['dE(eV)'])
    LE_fragment = int(LE_fragment) 
    if LE_fragment != 1 and LE_fragment !=2 :
        print('Error: LE should be on the first or the second fragment in you defined theodore analysis')
        return -1
    else: 
        pass 
    
    POS_list , POS_sorted= sort_and_index(df=data, column_name='POS')
    LE_list=[]
    P_LE_index=[]
    high_OS=[]
    high_OS_index=[]
    PR_NTO_GOOD=[]
    LE_energy_list=[]
    # Criteria 1: POS
    for i in range(len(POS_sorted)):
        diff_POS = abs((POS_sorted[i]-LE_fragment))
        if  diff_POS <= POSthres :
            #print("diff_POS = ", diff_POS)
            # screen for state with local excitation on the LE_fragment 
            LE_list.append(POS_sorted[i])
            POS_candidate_list = list(data['POS'])
            P_LE_index.append(POS_candidate_list.index(POS_sorted[i]))
            #print("PI_candidate", POS_sorted[i])
    
    # Criteria 2: Oscillation strength       
    for j in P_LE_index :
        f_list = list(data['f'])
        if f_list[j] > float(OSthres): # The state is an absorption state
            high_OS.append(f_list[j])
            high_OS_index.append(high_OS.index(f_list[j]))
        #pi_state_candidate = (np.min(LE_list)) # based on lowest excitation porphyrin
    LE_state_candidate = (np.max(high_OS)) # based on oscillator strength
    LE_state = f_list.index(LE_state_candidate) +1 # +1 to convert python indexing to real electronic state

    # Criteria 3: Lowest lying PI-PI state, aborted 
    #for k in high_OS_index :
    #    pi_energy_list.append(energy_list[k])
    #    pi_state_candidate = (np.min(pi_energy_list))   
    #    pi_state = energy_list.index(pi_state_candidate) +1 # +1 for index to real state
        
    
    
    

    LE_energy = energy_list[LE_state -1] # -1 for real state index to index
    if verbose == True:
        print('Potentially LE: ',LE_list)
        print('Potentially LE_index: ',np.array(P_LE_index)+1)
        print('Highest OS: ',high_OS)
        print('Highest OS state index: ',high_OS_index)
        # print('GOOD_PR_NTO: ',PR_NTO_GOOD)
        print('LE energy = ',LE_energy, ' eV')
    else :
        print('LE energy = ',LE_energy, ' eV')
    print('--------------------------------')
    ################ Finding CT1 State ###########################
    
    CT_net_list , CT_net_sorted= sort_and_index(df=data, column_name='CTnt')
    CT_net_candidate=[]
    CT_net_index=[]
    CT1_energy_list=[]

    # Criteria 1: CTnet
    for i in range(len(CT_net_sorted)):
        # diff_CT_net = abs(CT_net_sorted[i]+1) # because CT_net_sorted[i] is all negative value
        diff_CT_net = abs(CT_net_sorted[i]) 
        if diff_CT_net >= float(CT1_thres) :
        # if  diff_CT_net <=CT1_thres :
            #print("diff_CT_net = ", diff_CT_net)
            CT_net_candidate.append(CT_net_sorted[i])
            CT_net_candidate_list = list(data['CTnt'])
            CT_net_index.append(CT_net_candidate_list.index(CT_net_sorted[i]))
            #print("PI_candidate", POS_sorted[i])
    
    # Criteria 2: argmin(CT1_energy-PI_energy)
    LECT1_gap = []
    LECT1_ind = []
    for j in CT_net_index:
        if j != LE_state - 1:
            CT1_energy = abs(energy_list[j] - LE_energy)
            # print('CT1_energy : ', CT1_energy, j)
            LECT1_gap .append(CT1_energy)
            CT1_index = energy_list.index(energy_list[j])
            CT1_energy_list.append(energy_list[CT1_index])
            LECT1_ind . append(j)
    
    # CT1_candidate = np.min(CT1_energy_list) 
    CT1_gap_min = np.min(LECT1_gap)
    # CT1_state=  energy_list.index(CT1_candidate) +1 # +1 for index to real state
    CT1_state = LECT1_ind[LECT1_gap.index(CT1_gap_min)] + 1 # for readable results
    # print('CT1STATE  ',CT1_state)
    CT1_energy = energy_list[CT1_state -1] # -1 for real state to index
    
    if verbose == True:
        print('CT1 index: ', CT1_state)
        print('CT1 candidates: ',CT_net_candidate)
        print('CT1 candidates index: ',CT_net_index)
        print('CT1 candidate energy: ',CT1_energy_list)
        print ('CT1 energy:',CT1_energy, ' eV')
    else :
        print ('CT1 energy:',CT1_energy, ' eV')
    print('--------------------------------')
    ################ Finding CT2 State ###########################    
    CT_list , CT_sorted= sort_and_index(df=data, column_name='CT')
    CT_candidate=[]
    CT_index=[]
    CT2_energy_list=[]

    # # Criteria 1: CT value

    # empty list prevention # do while loop
    target =1 # 1.1 # main target is CT value = 1 , add buffer to start the while loop
    bottom_lim = 0.95
    for i in range(len(CT_sorted)):
        tol1 = target - CT2_thres
        tol2 = target + CT2_thres
        if (target <= CT_sorted[i] and CT_sorted[i] < tol2) or (target >= CT_sorted[i] and CT_sorted[i] > tol1):
            CT_candidate.append(CT_sorted[i])
            CT_candidate_list = list(data['CT'])
            CT_index.append(CT_candidate_list.index(CT_sorted[i]))
    # while (len(CT_index) < 1) and target > bottom_lim:
    #     target = target - 0.01 # refine the search slowly
    #     for i in range(len(CT_sorted)):
    #         diff_CT = abs(CT_sorted[i]-target) # that is CT value closest to 1
    #         print(i,CT_sorted[i],diff_CT)
    #         if  diff_CT <= CT2_thres :
    #             print(CT2_thres,'THRES 2')
    #             #print("diff_CT = ", diff_CT)
    #             CT_candidate.append(CT_sorted[i])
    #             CT_candidate_list = list(data['CT'])
    #             CT_index.append(CT_candidate_list.index(CT_sorted[i]))


    # Criteria 2: argmin(CT2_energy-PI_energy)
    for j in CT_index:
        CT2_energy = abs(energy_list[j] - LE_energy)
        #print('CT2_energy : ', CT2_energy)
        CT2_index = energy_list.index(energy_list[j])  
        CT2_energy_list.append(energy_list[CT2_index])
        # consider selecting maximum CT ???
    
    # Creteria 3: 



    CT2_candidate = np.min(CT2_energy_list) 
    CT2_state=  energy_list.index(CT2_candidate) +1 # +1 for index to real state
    CT2_energy = energy_list[CT2_state -1] # -1 for real state to index
    if verbose == True :
        print('CT2 index: ',CT_candidate)
        print('CT2 candidates index: ',CT_index)
        print('CT2 candidates energy (eV) : ',CT2_energy_list)
        print ('CT2 energy: ', CT2_energy, ' eV')
    else :
        print ('CT2 energy:',CT2_energy , ' eV')

    print('--------------------------------')
    print('FINAL DECISION OF STATES:')
    print("The LE  state index is ", LE_state)
    print("The CT1 state index is ", CT1_state)
    print("The CT2 state index is ", CT2_state)
    
    result = {'GR': 0, 'LE': LE_state, 'CT1':CT1_state, 'CT2': CT2_state}
    energies =[0, LE_energy,CT1_energy, CT2_energy] # list of energy in eV
    return result, energies
 
def Find_CT_states_DA(data, LE_fragment, charge_frag_list, Max_CT_states, Include_bCT_states, POSthres=0.2, OSthres= 0.1, CT_thres=0.5,verbose=True):
    """
    Algorithm to find the LE and CT states 
    LE_fragment defines the locally excited fragments in triad or dyad; int; for triad case; LE_fragment=2 for porphyrin excited state
    
    OSthres: threshold of the Oscillator Strength (OS)
    POSthres: thereshold of the POS (averaged position)
    CT1_thres: threshold of the charge transfer from C to PC60
    
    return 
    list of CT states final decision, and the corresponding energies 
    """
    # Need to load quantum chemistry results 
    
    ################ Finding LE State ###########################
    energy_list= list(data['dE(eV)'])
    LE_fragment = int(LE_fragment) 
    if LE_fragment != 1 and LE_fragment !=2 :
        print('Error: LE should be on the first or the second fragment in you defined theodore analysis')
        return -1
    else: 
        pass 
    
    POS_list , POS_sorted= sort_and_index(df=data, column_name='POS')
    LE_list=[]
    P_LE_index=[]
    high_OS=[]
    high_OS_index=[]
    PR_NTO_GOOD=[]
    LE_energy_list=[]
    POSthres = float(POSthres)


    if Include_bCT_states == "initial":
        print("Consider bCT as an initial state (bCT as an LE state candidate) ")
    elif Include_bCT_states == "final":
        print("Consider bCT as a final state (bCT as a CT state candidate) ")
    elif Include_bCT_states == "initial_or_final":
        print("Consider bCT as an initial or final state (bCT can be the LE or CT states candidate)")

    # Criteria 1: POS
    for i in range(len(POS_sorted)):
        diff_POS = abs((POS_sorted[i]-LE_fragment))
        # print(POSthres, type(POSthres))
        if  diff_POS <= POSthres :
            #print("diff_POS = ", diff_POS)
            # screen for state with local excitation on the LE_fragment 
            LE_list.append(POS_sorted[i])
            POS_candidate_list = list(data['POS'])
            P_LE_index.append(POS_candidate_list.index(POS_sorted[i]))
            # print("PI_candidate", POS_sorted[i])

    
    

    # Criteria 1: Q_D amount

    ####################### CT DECISION ACCORDING TO Q_D ####################################################
    # determine CT states according to Q_D and energy ordering + Max_CT_states to determine the CT states
    ########################################################################################################@
    CT_candidate = []
    bdCT_list = []
    OS_list = list(data['f'])
    OSthres = float(OSthres)
    for i in range(len(charge_frag_list)): 
        # select the CT states: 1. dark (small OS), 2. CT (quite amound of charges on D segment in FCD calculation)
        if charge_frag_list[i] > float(CT_thres): #and OS_list[i] < OSthres :
            #if Include_bCT_states == "initial":
            CT_candidate.append(i)
            #elif Include_bCT_states == "initial_or_final":
            if OS_list[i] < OSthres:
                bdCT_list.append('d')
            else:
                bdCT_list.append('b')
            if Include_bCT_states == "initial" and OS_list[i] >= OSthres :
                bdCT_list = bdCT_list[:-1]
                CT_candidate = CT_candidate[:-1]
                
                #CT_candidate.append(i)

    ######################## BRIGHT STATE DECISION ##################################################################
   
    # Criteria 2: Oscillation strength       
    for j in P_LE_index :
        f_list = list(data['f']) #
        if f_list[j] > OSthres: # The state is an absorption state
            if Include_bCT_states == "initial" or Include_bCT_states == "initial_or_final":
                high_OS.append(f_list[j])
                high_OS_index.append(high_OS.index(f_list[j]))
            elif Include_bCT_states == "final": # bCT as CT states
                if charge_frag_list[high_OS.index(f_list[j])] <= float(CT_thres):
                    high_OS.append(f_list[j])
                    high_OS_index.append(high_OS.index(f_list[j]))
            # high_OS_index.append(high_OS.index(f_list[j]))
    if len(high_OS) == 0:
        raise RuntimeError("bCT is the only bright state within the OSthreshold, revise your fucking input accordingly")
    
    # DECIDE LE State
    LE_state_candidate = (np.max(high_OS)) # based on oscillator strength
    LE_state = f_list.index(LE_state_candidate) +1 # +1 to convert python indexing to real electronic state
    LE_energy = energy_list[LE_state -1] # -1 for real state to index

    if verbose == True:
        print('Potentially LE: ',LE_state_candidate)
        print('Potentially LE index: ',P_LE_index)
        print('Highest OS: ',high_OS)
        print('Highest OS state index: ',high_OS_index)
        print('LE energy = ',LE_energy, ' eV')
        print('bdCT: ', bdCT_list)
    else :
        print('LE energy = ',LE_energy, ' eV')
        print('bdCT: ', bdCT_list)
    print('--------------------------------')



    
    

    # print(type(charge_frag_list[0]),Max_CT_states.type)

    
    
    Max_CT_states = int(Max_CT_states)
    # Criteria 2: energy ordering
    if len(CT_candidate) < Max_CT_states:
        print('The number of CT states is less than the maximum number of CT states '+len(CT_candidate))
        CT_state_order = CT_candidate
    else:
        CT_state_order = CT_candidate[0:Max_CT_states]
    # get homosapiens readable format for CT states
    CT_states = []
    
    for i in CT_state_order:
        CT_states.append(i+1) # +1 to convert python indexing to real electronic state
    # get energy of CT states in eV
    CT_energy_list=[] 
    for i in CT_state_order:
        CT_energy_list.append(energy_list[i])

    # get the charge of CT states on donor segment
    CT_donor_charge = []
    for i in CT_state_order:
        CT_donor_charge.append(charge_frag_list[i])

    # get the charge charge on the donor segment for different CT states
    # dQ_D = Q_D of the i-th CT state - Q_D of the LE state

    CT_dQ_D = []
    for i in CT_state_order:
        CT_dQ_D.append(charge_frag_list[i] - charge_frag_list[LE_state-1])
    



    if verbose == True:
        print('CT states index: ',CT_states)
        print('CT donor charge: ',CT_donor_charge)
        print('CT dQ_D (e): ',CT_dQ_D)
        print('CT energy (eV): ',CT_energy_list)
        
    else :
        print ('CT states index: ', CT_states)
        print ('CT energy (eV): ',CT_energy_list)
    print('--------------------------------')
    

    print('---------------------------------------------------------')
    print('FINAL DECISION OF STATES')
    print("The CT state index: ", CT_states)
    print("The CT state energy (eV): ", CT_energy_list)

    result = {'GR' : 0, 'LE': LE_state}
    for ind, state in enumerate(CT_states):
        result['CT'+str(ind+1)] = state
    energies =[0, LE_energy] # list of energy in eV
    for ct_ene in CT_energy_list:
        energies.append(ct_ene)
    return result, energies,  CT_dQ_D, CT_donor_charge, bdCT_list 



def get_charges(N, STATE, CASE, filedir, charge_dir, name, prefname, charge_type, Drange, Brange, Arange, QChem_version, scaling = None , verbose = True):
    
    """
    hartree2kcalmol = 627.51 
    Angs2Bohr = 0.5291772
    # in Amber.parm E(in Kcal/mol) = ( q1 * q2 ) /r(in Angs)
    # q in electron charge unit
    Scale=np.sqrt(hartree2kcalmol*Angs2Bohr) # scaling for amber formatted charge
    
    """
    ordersD  = np.fromstring(Drange, dtype=int, sep='-')#Drange.split(',')
    ordersA  = np.fromstring(Arange, dtype=int, sep='-') #
    ordersB  = np.fromstring(Brange, dtype=int, sep='-')#D 
    systemSeg = ''
    if len(ordersB) == 0:
        if ordersA[1] + ordersD[1] - ordersA[0] - ordersD[0] + 2 == N:
            systemSeg = systemSeg + 'DA'
        else: 
            systemSeg  
    
        if ordersD[1] > ordersA[1]:
            systemSeg = 'AD'
        else:
            systemSeg = 'DA'
        
    # print('segments in atmoic order: ', systemSeg)
    output_file= f'{filedir}/{name}.out'
    print('QC result file: ' ,output_file)
    #Reading line by line from the output file
    lines = []
    with open(output_file) as f:
        lines = f.readlines()
    #print(lines[-11])
    if QChem_version == '4':
        skip = -11
    elif QChem_version == '6':
        skip = -7
    qchem4int = -11
    qchem6int = -7
    while (lines[skip] != '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n'):
        print("Thank you...Q-Chem... line was not found but we found this line: \n",lines[skip])
        time.sleep(1)
        with open(output_file) as f:
            lines = f.readlines()
        
    #Find section of interest 
    if charge_type == 'Mulliken':
        if STATE == 0:
            keywords_to_find = '          Ground-State Mulliken Net Atomic Charges\n'
        elif STATE < 10:
            keywords_to_find = '         TDA Excited State  '+ str(STATE) +':  Mulliken Net Atomic Charges\n'
        else: 
            keywords_to_find = '         TDA Excited State '+ str(STATE) +':  Mulliken Net Atomic Charges\n'
    else:
        raise Exception(charge_type + ' are NOT defined in CTRAMER' )
        return -1

    # capture the section
    count = 0
    lines_after_title=3 # this is QCHEM 4.4. formatted mulliken charge 
    my_charge_section=[]
    for line in lines:
        count += 1
        if line == keywords_to_find :
            start = count -1 # to include the title too
            end = count + lines_after_title + N
    print('Solute atomic charges from QC result file between lines:')
    print('start line = ', start+4)
    print('end line = ',end)
    contents=lines[start+4:end]
    #print(contents)
    charge_data  = contents
    ############################################################
    savefile = f'{charge_dir}/charge_{CASE}_S{str(STATE)}.txt'
    savefile2 = f'{charge_dir}/charge_{CASE}_'+prefname+'.txt'
    
    with open(savefile, "w") as output:
        #for row in contents:
        #    output.write(str(row) + '\n')
        if systemSeg == "AD":
            for i in range(ordersD[1] - ordersD[0] + 1):
                output.write(str(i+1) + ' ' + str(charge_data[i+ordersA[1]-ordersA[0]+1]) + '\n')
            for i in range(ordersA[1] - ordersA[0] + 1):
                output.write(str(i+1+ordersD[1]-ordersD[0]+1) + ' ' + str(charge_data[i]) + '\n')
        else:
            for ind, row in enumerate(charge_data):
                output.write(str(ind+1) + ' ' + str(row) + '\n')
    
    Qchem_ver = 6
    if Qchem_ver == 6:
        charge_table= pd.read_table(savefile, skiprows= 0, delim_whitespace=True,names=["Order", "No_Atom", "Atoms", '(a.u.)'] )
        print('Headers of charge files: ', charge_table.columns)
    else:
        charge_table= pd.read_table(savefile, skiprows= 0, delim_whitespace=True,names=["Order", "No_Atom", "Atoms", '(a.u.)'] )
        print('Headers of charge files: ', charge_table.columns)
    
    hartree2kcalmol = 627.51 
    Angs2Bohr = 0.5291772
    # in Amber.parm E(in Kcal/mol) = ( q1 * q2 ) /r(in Angs)
    # q in electron charge unit
    Scale=np.sqrt(hartree2kcalmol*Angs2Bohr) # scaling for amber formatted charge
    
    charge_data=charge_table['(a.u.)']#[1:]
    with open(savefile2, "w") as output:
        count = 0
        if systemSeg == "AD":
            for ind, row in enumerate(charge_data):
                output.write(str(ind+1) + ' ' + str(row) + '\n')
            #for i in range(ordersD[1] - ordersD[0] + 1):
            #    output.write(str(i+1) + ' ' + str(charge_data[i+ordersA[1]-ordersA[0]+1]) + '\n')
            #for i in range(ordersA[1] - ordersA[0] + 1):
            #    output.write(str(i+1+ordersD[1]-ordersD[0]+1) + ' ' + str(charge_data[i]) + '\n')
        else:
            for ind, row in enumerate(charge_data):
                output.write(str(ind+1) + ' ' + str(row) + '\n')
        #for ind, row in enumerate(charge_data):
        #    count += 1
        #    if scaling == "AMBER" : 
        #        output.write(str(count) + ' ' + str(row*Scale) + '\n')
        #    else: 
        #        output.write(str(count) + ' ' + str(row) + '\n')
    if verbose == True :
        print('The charges have been saved to the following two files:')
        print(savefile, '  with headers')
        print(savefile2, '  without headers')

def replace_charge_to_prmtop(No_Atoms, newfiledir, old_prmtop, charge_file, new_name):
    """
    replace calculated charges to amber force field files  
    """
    
    hartree2kcalmol = 627.51 
    Angs2Bohr = 0.5291772
    # in Amber.parm E(in Kcal/mol) = ( q1 * q2 ) /r(in Angs)
    # q in electron charge unit
    Scale=np.sqrt(hartree2kcalmol*Angs2Bohr) # scaling for amber formatted charge

    #rounding up 
    nrow = math.ceil(No_Atoms/5) # in amber charge section are divided into 5 columns

    ######### Defining input output file 
    
    new_prmtop=f'{newfiledir}/{new_name}.prmtop'
    # print("Generating state-depedent force field:\n",new_prmtop)
    chargefile = charge_file
    
    ochg = []
    inprm = open(old_prmtop,'r')
    
    lines = inprm.readlines()
    m=0
    for line in lines:    
        res1 = re.findall(r'%FLAG CHARGE',line)
        if res1 !=[]:
            n=lines.index(line)
            m = n 
    #print(m)

    for i in range(m+2,m+2+nrow):
        tt=lines[i].split()
        for i in range(len(tt)):
            ochg.append(float(tt[i]))
    
    #print('.prmtop to be replaced =', old_prmtop)
    #################################################
    
    #for i in charge_list:
        #print('replacing charge with charge on ', i , 'state')
        
    inchg=np.loadtxt(charge_file,unpack=True,usecols=[1])
        
    new_prmtop = new_prmtop
    #print('new .prmtop file =', new_prmtop)
        
    outprm=open(new_prmtop,'w')
    #========replace charge==========            
    for i in range(No_Atoms):
        ochg[i]=inchg[i]*Scale

    n0=0
    for i in range(len(lines)):
        if i>m+1 and i<m+2+nrow:
            for n in range(5):
                k=n0*5+n        
                outprm.write('{:16.8E}'.format(ochg[k]))
            n0 +=1
        #        print(n0)
            outprm.write('\n')
        else:
            outprm.write('{}'.format(lines[i]))
    print(arrow_str,"Generated state-dependent force field file: \n", new_prmtop, '\n')

def get_coupling(qchem_output_dir,name, state1, state2, Nstate): 
    Nstate = int(Nstate)
    """
    acquiring FCD coupling for state of interest 
    
    """
    
    output_file= f'{qchem_output_dir}/{name}.out'
    outCoup = f'{qchem_output_dir}/{name}_coupling.out'
    
    #Reading line by line from the output file
    lines = []
    with open(output_file) as f:
        lines = f.readlines()
        
        
    contents = []
    count = 0
    keywords_to_find = '            FCD Couplings Between Singlet Excited States\n'
    
    # total pairwise combination (N) = n_factorial/(n_min_r_factorial * 2); n = Nstate
    N = int(math.factorial(Nstate) / (2 * math.factorial(Nstate-2))) 
    # print(N)
    skip_header = 3
    include_footer = 7
    #N = 310 # resolve this total number of lines
    for line in lines:
        count += 1
        if line == keywords_to_find :
            start = count  + skip_header # to include the title too
            end = count + N  + skip_header
    
    #print('start = ', start)
    #print('end = ',end)
    contents=lines[start:end]
    raw_start =start - skip_header
    raw_end  = end -skip_header + include_footer
    raw_contents = lines[raw_start:raw_end]
#     ############################################################
    savefile = f'{qchem_output_dir}/FCD_coupling_raw.txt'
    #savefile2 = f'{filedir}/FCD_coupling_{CASE}.txt'
    
    with open(savefile, "w") as output:
        for row in contents:
            output.write(str(row) + '\n')
    test_table= pd.read_table(savefile, header=None,delim_whitespace=True )
    a = state1
    b = state2
    coupling = [] # FCD coupling squared
    for i in range(len(test_table)):
        #column[0,1,2,3,4,5,6,7] = state A, state B , x, y , z, FCD Coupling(ev),  
        if test_table[0][i] == a and test_table[1][i] == b :
            c = test_table[5][i]
            print("Diabatic coupling between states ", a, " and " , b , " = ", c, ' eV')
            # print("coupling^2 coefficient between state ", a, " and " , b , " : ", c*c)
            coupling.append(c)
        elif test_table[0][i] == b and test_table[1][i] == a :
            c = test_table[5][i]
            print("Diabatic coupling between states ", a, " and " , b , " = ", c, ' eV')
            # print("coupling^2 coefficient between state ", a, " and " , b , " : ", c*c)
            coupling.append(c)
    return coupling


def construct_prmtop_inter(template_dir,running_dir,testing_str,case,Drange,Arange,Brange,
                           CT_donor_fragments,CT_acceptor_fragment,
                           solute_dir,solvent_dir,
                           solute_resname,solvent_resname,solvent_num,
                           box_size,center_of_mass,force_field_type,N_cis,charge_file, dependency=''):
    """
    Create prmtop template with net zero charge   
    1. prepare the pdb file by their own in the target running folder first
    
    template_dir: slrum script template for the project 
    running_dir:  where to run this quantum chemistry jobs and prepare the force field 
    testing_str: a string containing several space separated values utilized as part of the name of quantum chemistry jobs
    case: name to specify the purpose or subject of calculation like `triad`
    solvent_dir: solvent related quantum chemistry calculation results and structure files 
    solvent_name: name abbriviation of the solvent, usually 3 captial letters
    solute_name: name abbriviation of the solute, usually short and understandable
    
    solvent_num: number of solvent, int
    box_size: length of cubic simulation box, unit in A
    center_of_mass: center of mass of the solute with 3-D vector coordinate (center_of_mass,center_of_mass,center_of_mass), unit in A
    force_field_type: name of the force field utilized to describe the system
    charge_file: zero-charge file or state-dependent charge file utilized to constrcut the FF
    
    N_cis: number of singlets are calculated, int 
    case: series name for the calculation 
    
    This function prepare and submit prmtop construction job based on a construct_prmtop_template.sh 
    
    other parameter to control the construc prmptop such as solvent type, no of solvent molecule 
    box size , etc will be included in next update

    check the template for platform and cluster specific parameters
    in the next release the cluster control can be performed from external function ctr.job_submit_setup()
    """

    construct_prmtop_file = shell_copy_template(template_dir,running_dir, filename='construct_prmtop_inter')
    os.chdir(running_dir)
    QC_DIR = f'{running_dir}'
    if len(solute_dir.split(",")) == 1:
        A_dir = solute_dir
        D_dir = solute_dir
    elif len(solute_dir.split(",")) == 2:
        D_dir = solute_dir.split(",")[0]
        A_dir = solute_dir.split(",")[1]
    if len(solvent_dir.split(",")) == 1:
        a_dir = solvent_dir
        d_dir = solvent_dir
    elif len(solvent_dir.split(",")) == 2:
        d_dir = solvent_dir.split(",")[0]
        a_dir = solvent_dir.split(",")[1] 
    # d_dir = solvent_dir.split(",")[0]
    # a_dir = solvent_dir.split(",")[1]
    
    if len(solute_resname.split(",")) == 2:
        D_res = solute_resname.split(",")[0] 
        A_res = solute_resname.split(",")[1]
    if len(solvent_resname.split(",")) == 2:
        d_res = solvent_resname.split(",")[0]
        a_res = solvent_resname.split(",")[1]

    solvent_name = d_res + a_res
    solute_name = solute_name1 = solute_resname
    sed_inplace(construct_prmtop_file,'WORKDIR=','WORKDIR='+QC_DIR)
    sed_inplace(construct_prmtop_file,'SRCDIR=','SRCDIR='+template_dir)
    sed_inplace(construct_prmtop_file,'RUNDIR=','RUNDIR='+QC_DIR)
    sed_inplace(construct_prmtop_file,'system=','system='+case+'_')
    sed_inplace(construct_prmtop_file,'GIVEN_STRUC=','GIVEN_STRUC='+ testing_str)
    
    sed_inplace(construct_prmtop_file,"D_resname=","D_resname="+D_res)
    sed_inplace(construct_prmtop_file,"A_resname=","A_resname="+A_res)
    sed_inplace(construct_prmtop_file,"d_resname=","d_resname="+d_res)
    sed_inplace(construct_prmtop_file,"a_resname=","a_resname="+a_res)
    
    sed_inplace(construct_prmtop_file,"D_dir=","D_dir="+D_dir)
    sed_inplace(construct_prmtop_file,"A_dir=","A_dir="+A_dir)
    sed_inplace(construct_prmtop_file,"d_dir=","d_dir="+d_dir)
    sed_inplace(construct_prmtop_file,"a_dir=","a_dir="+a_dir)
    
    sed_inplace(construct_prmtop_file,'forcefield=gaff','forcefield='+force_field_type)
    sed_inplace(construct_prmtop_file,'box_side=60','box_side='+str(box_size))
    
    sed_inplace(construct_prmtop_file,'solvent_mol=2700','solvent_mol='+str(solvent_num))
    sed_inplace(construct_prmtop_file,'solvent_dir=','solvent_dir='+solvent_dir)
    sed_inplace(construct_prmtop_file,'solvent=','solvent='+solvent_name)
    
    sed_inplace(construct_prmtop_file,'com=30','com='+str(center_of_mass))
    sed_inplace(construct_prmtop_file,'chargefile=','chargefile='+charge_file)
    sed_inplace(construct_prmtop_file,'code=','code='+solute_name)
    sed_inplace(construct_prmtop_file,'solute=','solute='+solute_name1)

    sed_inplace(construct_prmtop_file,'STATE_NUM=','STATE_NUM='+str(N_cis))
    
    os.chdir(running_dir)
    # if dependency == '':
    #     data = os.popen('sbatch construct_prmtop_flexible.sh').read()
    # else :
    #     data = os.popen('sbatch --dependency=afterok:'+dependency+ ' construct_prmtop_flexible.sh').read()
    data = os.popen('sh construct_prmtop_inter.sh').read()
    #os.chdir(basedir)
    return data 

def construct_prmtop_intra(template_dir,running_dir,testing_str,case,Drange,Arange,Brange,
                           CT_donor_fragments,CT_acceptor_fragment,
                           solute_dir,solvent_dir,
                           solute_resname,solvent_resname,solvent_num,
                           box_size,center_of_mass,force_field_type,N_cis,charge_file, dependency=''):
    """
    Create prmtop template with net zero charge   
    1. prepare the pdb file by their own in the target running folder first
    
    template_dir: slrum script template for the project 
    running_dir:  where to run this quantum chemistry jobs and prepare the force field 
    testing_str: a string containing several space separated values utilized as part of the name of quantum chemistry jobs
    case: name to specify the purpose or subject of calculation like `triad`
    solvent_dir: solvent related quantum chemistry calculation results and structure files 
    solvent_name: name abbriviation of the solvent, usually 3 captial letters
    solute_name: name abbriviation of the solute, usually short and understandable
    
    solvent_num: number of solvent, int
    box_size: length of cubic simulation box, unit in A
    center_of_mass: center of mass of the solute with 3-D vector coordinate (center_of_mass,center_of_mass,center_of_mass), unit in A
    force_field_type: name of the force field utilized to describe the system
    charge_file: zero-charge file or state-dependent charge file utilized to constrcut the FF
    
    N_cis: number of singlets are calculated, int 
    case: series name for the calculation 
    
    This function prepare and submit prmtop construction job based on a construct_prmtop_template.sh 
    
    other parameter to control the construc prmptop such as solvent type, no of solvent molecule 
    box size , etc will be included in next update

    check the template for platform and cluster specific parameters
    in the next release the cluster control can be performed from external function ctr.job_submit_setup()
    """

    construct_prmtop_file = shell_copy_template(template_dir,running_dir, filename='construct_prmtop_intra')
    os.chdir(running_dir)
    QC_DIR = f'{running_dir}'
    if len(solute_dir.split(",")) == 1:
        A_dir = solute_dir
        D_dir = solute_dir
    elif len(solute_dir.split(",")) == 2:
        D_dir = solute_dir.split(",")[0]
        A_dir = solute_dir.split(",")[1]
    if len(solvent_dir.split(",")) == 1:
        a_dir = solvent_dir
        d_dir = solvent_dir
    elif len(solvent_dir.split(",")) == 2:
        d_dir = solvent_dir.split(",")[0]
        a_dir = solvent_dir.split(",")[1] 
    # d_dir = solvent_dir.split(",")[0]
    # a_dir = solvent_dir.split(",")[1]
    
    if len(solute_resname.split(",")) == 2:
        D_res = solute_resname.split(",")[0] 
        A_res = solute_resname.split(",")[1]
    elif len(solute_resname.split(",")) == 1:
        D_res  = solute_resname
        A_res  = ''
    
    if len(solvent_resname.split(",")) == 2:
        d_res = solvent_resname.split(",")[0]
        a_res = solvent_resname.split(",")[1]
    elif len(solvent_resname.split(",")) == 1:
        d_res = solvent_resname
        a_res = ''

    solvent_name = d_res + a_res
    solute_name = solute_name1 = solute_resname
    sed_inplace(construct_prmtop_file,'WORKDIR=','WORKDIR='+QC_DIR)
    sed_inplace(construct_prmtop_file,'SRCDIR=','SRCDIR='+template_dir)
    sed_inplace(construct_prmtop_file,'RUNDIR=','RUNDIR='+QC_DIR)
    sed_inplace(construct_prmtop_file,'system=','system='+case+'_')
    sed_inplace(construct_prmtop_file,'GIVEN_STRUC=','GIVEN_STRUC='+ testing_str)
    
    sed_inplace(construct_prmtop_file,"D_resname=","D_resname="+D_res)
    sed_inplace(construct_prmtop_file,"A_resname=","A_resname="+A_res)
    sed_inplace(construct_prmtop_file,"d_resname=","d_resname="+d_res)
    sed_inplace(construct_prmtop_file,"a_resname=","a_resname="+a_res)
    
    sed_inplace(construct_prmtop_file,"D_dir=","D_dir="+D_dir)
    sed_inplace(construct_prmtop_file,"A_dir=","A_dir="+A_dir)
    sed_inplace(construct_prmtop_file,"d_dir=","d_dir="+d_dir)
    sed_inplace(construct_prmtop_file,"a_dir=","a_dir="+a_dir)
    
    sed_inplace(construct_prmtop_file,'forcefield=gaff','forcefield='+force_field_type)
    sed_inplace(construct_prmtop_file,'box_side=60','box_side='+str(box_size))
    
    sed_inplace(construct_prmtop_file,'solvent_mol=2700','solvent_mol='+str(solvent_num))
    sed_inplace(construct_prmtop_file,'solvent_dir=','solvent_dir='+solvent_dir)
    sed_inplace(construct_prmtop_file,'solvent=','solvent='+solvent_name)
    
    sed_inplace(construct_prmtop_file,'com=30','com='+str(center_of_mass))
    sed_inplace(construct_prmtop_file,'chargefile=','chargefile='+charge_file)
    sed_inplace(construct_prmtop_file,'code=','code='+solute_name)
    sed_inplace(construct_prmtop_file,'solute=','solute='+solute_name1)

    sed_inplace(construct_prmtop_file,'STATE_NUM=','STATE_NUM='+str(N_cis))
    
    os.chdir(running_dir)
    # if dependency == '':
    #     data = os.popen('sbatch construct_prmtop_flexible.sh').read()
    # else :
    #     data = os.popen('sbatch --dependency=afterok:'+dependency+ ' construct_prmtop_flexible.sh').read()
    data = os.popen('sh construct_prmtop_intra.sh').read()
    #os.chdir(basedir)
    return data #    return

def get_charge_on_Donor(qchem_output_dir, name, Nstate):
    output_file= f'{qchem_output_dir}/{name}.out'
    str_to_find = '       Fragment Charges of Singlet Excited State with Nuclear Charges\n'
    #Reading line by line from the output file
    with open(output_file) as f:
        lines = f.readlines()
    start,end = 0,0
    for line in lines:
        if line == str_to_find :
            start = lines.index(line)+4
            end   = lines.index(line)+Nstate+4
    data = lines[start:end]
    return data


import pickle 

def main_ConstructFF():
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
    job_info = init_job_info(case_id_list)
    startime = time.time 
    print("ConstructFF begins at "+str(startime))
    # print("Input parameters ")
    # print(simulation_parameter)
    for ind, ctr_ele in enumerate(job_control): 
        print(ind, job_control)
        run_step3(ind, job_control, dict_of_simulation,job_info)
        
        

    print(job_info)
    # write list to binary file
    def write_list(a_list):
        # store list in binary file so 'wb' mode
        with open(work_dir+'/listfile', 'wb') as fp:
            pickle.dump(a_list, fp)
            print('Done writing list into a binary file')

    # Read list to memory
    def read_list():
        # for reading also binary mode is important
        with open(work_dir+'/listfile', 'rb') as fp:
            n_list = pickle.load(fp)
            return n_list


    write_list(job_info)
    job_info1 = read_list()
    print('info1 \n', job_info1)
    # run_submitQC(simulation_parameter)
    end_time = time.time  
    print("ConstructFF ends at "+str(end_time))

if __name__ == "__main__":
    main_ConstructFF()
