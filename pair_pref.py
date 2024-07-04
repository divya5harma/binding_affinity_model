

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 16:22:33 2024

@author: divya
"""

############# Binding propensity ##########
import pandas as pd
import os
from Bio import PDB
import numpy as np
pdb_path = '/home/divya/Documents/ab_analysis/pdbs_all_final/'
cognate_path = '/home/divya/Documents/ab_analysis/cognate_recps/'
# all_files = os.listdir(pdb_path)
# # Filter PDB files
# pdb_files = [file for file in all_files if file.endswith('.pdb.csv')]
# pdb_files = pdb_files[0:5]
# print(pdb_files)
# pdb_files_seq = [file for file in all_files if file.endswith('.pdb')]
# pdb_files_seq = pdb_files_seq[0:5]
# print(pdb_files_seq)
# df1 = pd.read_csv(pdb_path+'all_interaction_6C6Z.pdb.csv',header=None)
# epis = df1.iloc[0]
# print(epis)

# dfz = pd.read_csv('pdb_all_files.csv')
# dfz = pd.read_excel('virus_pdbs.xlsx',sheet_name='protozoa')
dfz = pd.read_csv('cognate_groups2.csv')

## get epitope and paratope sequence
for x in range(len(dfz)):
    cog_pdb = dfz['cognate_pdb'][x].lower()
    f_ab = dfz['ab_pdb'][x]
    f_aa = f_ab.split(',')
    print(f_aa)
    p = []
    e = []
    pdb = []
    pairs = []
    ab_prot = []
    epi_prot = []
    for i in range(len(f_aa)):
        pairs.append([])
        f = f_aa[i].lower()
        res_para_seq = ''
        res_epi_seq = ''
        file = open(pdb_path+'all_interaction_'+f+'.pdb.csv','r')
        file_detail = file.readlines()
        for x,line in enumerate(file_detail):
            l = line.split(',')
            l_clean = [value.strip('"\n') for value in l if value.strip('"') and value != '\n']
            if l_clean:
                ## epitope seq
                if l_clean[0] == '*** paratope contact list ***':
                    # print('E: ',l_clean)
                    contacts = file_detail[x+1:]
                    contacts = [value.strip('"\n') for value in contacts if value.strip('"') and value != '\n']
                    # print(contacts)
                    for j in contacts:
                        res = j.split(',')
                        p1 = res[0][-1]
                        # print(p1)
                        if (len(res))>1:
                            for y in res[1:]:
                                # print(res)
                                if y != '':
                                    p = p1 + y[-1]
                                    # print(p)
                                    pairs[i].append(p)
                        else:
                            pairs[i].append(p1)
                                     
                
                # Filter out empty strings and newline characters, and remove the double quotes
                # res_clean = [value.strip('"') for value in res if value.strip('"') and value != '\n']
                # res = [s.strip('\n') for s in res]
               
            
           
        pdb.append(f)
   

# ######## get sequence of protein from the pdb file


        parser = PDB.PDBParser(QUIET=True)
    
        # Load the PDB file
        structure = parser.get_structure('protein', pdb_path+f+'.pdb')
    
        # Initialize a dictionary to store sequences for each chain
        chain_sequences = {}
    
        # Iterate over all models in the structure
        for model in structure:
            # Iterate over all chains in the model
            for chain in model:
                # Initialize an empty string to store the sequence
                sequence = ''
    
                # Iterate over all residues in the chain
                for residue in chain:
                    # Check if the residue is an amino acid
                    if PDB.is_aa(residue):
                        # Append the one-letter code of the amino acid to the sequence
                        sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
    
                # Store the sequence in the dictionary with the chain ID as the key
                chain_sequences[chain.id] = sequence
    
        ab_seq = chain_sequences['H'] + chain_sequences['L']
        epi_prot_seq = chain_sequences['A']
        ab_prot.append(ab_seq)
        epi_prot.append(epi_prot_seq)
                    
    ####### get the pair preference with standard deviation
    dff = pd.DataFrame({'pdb':pdb,'ab_prot':ab_prot,'ag_prot':epi_prot})
    dict2_ag = {}
    for i in range(len(dff)):
        amino_acids_ab = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
        
        amino_acids_ag = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
        mains = {'AA': 0, 'AC': 0, 'AD': 0, 'AE': 0, 'AF': 0, 'AG': 0, 'AH': 0, 'AI': 0, 'AK': 0, 'AL': 0, 'AM': 0, 'AN': 0, 'AP': 0, 'AQ': 0, 'AR': 0, 'AS': 0, 'AT': 0, 'AV': 0, 'AW': 0, 'AY': 0, 'CA': 0, 'CC': 0, 'CD': 0, 'CE': 0, 'CF': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CK': 0, 'CL': 0, 'CM': 0, 'CN': 0, 'CP': 0, 'CQ': 0, 'CR': 0, 'CS': 0, 'CT': 0, 'CV': 0, 'CW': 0, 'CY': 0, 'DA': 0, 'DC': 0, 'DD': 0, 'DE': 0, 'DF': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DK': 0, 'DL': 0, 'DM': 0, 'DN': 0, 'DP': 0, 'DQ': 0, 'DR': 0, 'DS': 0, 'DT': 0, 'DV': 0, 'DW': 0, 'DY': 0, 'EA': 0, 'EC': 0, 'ED': 0, 'EE': 0, 'EF': 0, 'EG': 0, 'EH': 0, 'EI': 0, 'EK': 0, 'EL': 0, 'EM': 0, 'EN': 0, 'EP': 0, 'EQ': 0, 'ER': 0, 'ES': 0, 'ET': 0, 'EV': 0, 'EW': 0, 'EY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'GA': 0, 'GC': 0, 'GD': 0, 'GE': 0, 'GF': 0, 'GG': 0, 'GH': 0, 'GI': 0, 'GK': 0, 'GL': 0, 'GM': 0, 'GN': 0, 'GP': 0, 'GQ': 0, 'GR': 0, 'GS': 0, 'GT': 0, 'GV': 0, 'GW': 0, 'GY': 0, 'HA': 0, 'HC': 0, 'HD': 0, 'HE': 0, 'HF': 0, 'HG': 0, 'HH': 0, 'HI': 0, 'HK': 0, 'HL': 0, 'HM': 0, 'HN': 0, 'HP': 0, 'HQ': 0, 'HR': 0, 'HS': 0, 'HT': 0, 'HV': 0, 'HW': 0, 'HY': 0, 'IA': 0, 'IC': 0, 'ID': 0, 'IE': 0, 'IF': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IK': 0, 'IL': 0, 'IM': 0, 'IN': 0, 'IP': 0, 'IQ': 0, 'IR': 0, 'IS': 0, 'IT': 0, 'IV': 0, 'IW': 0, 'IY': 0, 'KA': 0, 'KC': 0, 'KD': 0, 'KE': 0, 'KF': 0, 'KG': 0, 'KH': 0, 'KI': 0, 'KK': 0, 'KL': 0, 'KM': 0, 'KN': 0, 'KP': 0, 'KQ': 0, 'KR': 0, 'KS': 0, 'KT': 0, 'KV': 0, 'KW': 0, 'KY': 0, 'LA': 0, 'LC': 0, 'LD': 0, 'LE': 0, 'LF': 0, 'LG': 0, 'LH': 0, 'LI': 0, 'LK': 0, 'LL': 0, 'LM': 0, 'LN': 0, 'LP': 0, 'LQ': 0, 'LR': 0, 'LS': 0, 'LT': 0, 'LV': 0, 'LW': 0, 'LY': 0, 'MA': 0, 'MC': 0, 'MD': 0, 'ME': 0, 'MF': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'MK': 0, 'ML': 0, 'MM': 0, 'MN': 0, 'MP': 0, 'MQ': 0, 'MR': 0, 'MS': 0, 'MT': 0, 'MV': 0, 'MW': 0, 'MY': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'NY': 0, 'PA': 0, 'PC': 0, 'PD': 0, 'PE': 0, 'PF': 0, 'PG': 0, 'PH': 0, 'PI': 0, 'PK': 0, 'PL': 0, 'PM': 0, 'PN': 0, 'PP': 0, 'PQ': 0, 'PR': 0, 'PS': 0, 'PT': 0, 'PV': 0, 'PW': 0, 'PY': 0, 'QA': 0, 'QC': 0, 'QD': 0, 'QE': 0, 'QF': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QK': 0, 'QL': 0, 'QM': 0, 'QN': 0, 'QP': 0, 'QQ': 0, 'QR': 0, 'QS': 0, 'QT': 0, 'QV': 0, 'QW': 0, 'QY': 0, 'RA': 0, 'RC': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RY': 0, 'SA': 0, 'SC': 0, 'SD': 0, 'SE': 0, 'SF': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SK': 0, 'SL': 0, 'SM': 0, 'SN': 0, 'SP': 0, 'SQ': 0, 'SR': 0, 'SS': 0, 'ST': 0, 'SV': 0, 'SW': 0, 'SY': 0, 'TA': 0, 'TC': 0, 'TD': 0, 'TE': 0, 'TF': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TK': 0, 'TL': 0, 'TM': 0, 'TN': 0, 'TP': 0, 'TQ': 0, 'TR': 0, 'TS': 0, 'TT': 0, 'TV': 0, 'TW': 0, 'TY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'VH': 0, 'VI': 0, 'VK': 0, 'VL': 0, 'VM': 0, 'VN': 0, 'VP': 0, 'VQ': 0, 'VR': 0, 'VS': 0, 'VT': 0, 'VV': 0, 'VW': 0, 'VY': 0, 'WA': 0, 'WC': 0, 'WD': 0, 'WE': 0, 'WF': 0, 'WG': 0, 'WH': 0, 'WI': 0, 'WK': 0, 'WL': 0, 'WM': 0, 'WN': 0, 'WP': 0, 'WQ': 0, 'WR': 0, 'WS': 0, 'WT': 0, 'WV': 0, 'WW': 0, 'WY': 0, 'YA': 0, 'YC': 0, 'YD': 0, 'YE': 0, 'YF': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YK': 0, 'YL': 0, 'YM': 0, 'YN': 0, 'YP': 0, 'YQ': 0, 'YR': 0, 'YS': 0, 'YT': 0, 'YV': 0, 'YW': 0, 'YY': 0}
        pair_pref = {'AA': 0, 'AC': 0, 'AD': 0, 'AE': 0, 'AF': 0, 'AG': 0, 'AH': 0, 'AI': 0, 'AK': 0, 'AL': 0, 'AM': 0, 'AN': 0, 'AP': 0, 'AQ': 0, 'AR': 0, 'AS': 0, 'AT': 0, 'AV': 0, 'AW': 0, 'AY': 0, 'CA': 0, 'CC': 0, 'CD': 0, 'CE': 0, 'CF': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CK': 0, 'CL': 0, 'CM': 0, 'CN': 0, 'CP': 0, 'CQ': 0, 'CR': 0, 'CS': 0, 'CT': 0, 'CV': 0, 'CW': 0, 'CY': 0, 'DA': 0, 'DC': 0, 'DD': 0, 'DE': 0, 'DF': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DK': 0, 'DL': 0, 'DM': 0, 'DN': 0, 'DP': 0, 'DQ': 0, 'DR': 0, 'DS': 0, 'DT': 0, 'DV': 0, 'DW': 0, 'DY': 0, 'EA': 0, 'EC': 0, 'ED': 0, 'EE': 0, 'EF': 0, 'EG': 0, 'EH': 0, 'EI': 0, 'EK': 0, 'EL': 0, 'EM': 0, 'EN': 0, 'EP': 0, 'EQ': 0, 'ER': 0, 'ES': 0, 'ET': 0, 'EV': 0, 'EW': 0, 'EY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'GA': 0, 'GC': 0, 'GD': 0, 'GE': 0, 'GF': 0, 'GG': 0, 'GH': 0, 'GI': 0, 'GK': 0, 'GL': 0, 'GM': 0, 'GN': 0, 'GP': 0, 'GQ': 0, 'GR': 0, 'GS': 0, 'GT': 0, 'GV': 0, 'GW': 0, 'GY': 0, 'HA': 0, 'HC': 0, 'HD': 0, 'HE': 0, 'HF': 0, 'HG': 0, 'HH': 0, 'HI': 0, 'HK': 0, 'HL': 0, 'HM': 0, 'HN': 0, 'HP': 0, 'HQ': 0, 'HR': 0, 'HS': 0, 'HT': 0, 'HV': 0, 'HW': 0, 'HY': 0, 'IA': 0, 'IC': 0, 'ID': 0, 'IE': 0, 'IF': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IK': 0, 'IL': 0, 'IM': 0, 'IN': 0, 'IP': 0, 'IQ': 0, 'IR': 0, 'IS': 0, 'IT': 0, 'IV': 0, 'IW': 0, 'IY': 0, 'KA': 0, 'KC': 0, 'KD': 0, 'KE': 0, 'KF': 0, 'KG': 0, 'KH': 0, 'KI': 0, 'KK': 0, 'KL': 0, 'KM': 0, 'KN': 0, 'KP': 0, 'KQ': 0, 'KR': 0, 'KS': 0, 'KT': 0, 'KV': 0, 'KW': 0, 'KY': 0, 'LA': 0, 'LC': 0, 'LD': 0, 'LE': 0, 'LF': 0, 'LG': 0, 'LH': 0, 'LI': 0, 'LK': 0, 'LL': 0, 'LM': 0, 'LN': 0, 'LP': 0, 'LQ': 0, 'LR': 0, 'LS': 0, 'LT': 0, 'LV': 0, 'LW': 0, 'LY': 0, 'MA': 0, 'MC': 0, 'MD': 0, 'ME': 0, 'MF': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'MK': 0, 'ML': 0, 'MM': 0, 'MN': 0, 'MP': 0, 'MQ': 0, 'MR': 0, 'MS': 0, 'MT': 0, 'MV': 0, 'MW': 0, 'MY': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'NY': 0, 'PA': 0, 'PC': 0, 'PD': 0, 'PE': 0, 'PF': 0, 'PG': 0, 'PH': 0, 'PI': 0, 'PK': 0, 'PL': 0, 'PM': 0, 'PN': 0, 'PP': 0, 'PQ': 0, 'PR': 0, 'PS': 0, 'PT': 0, 'PV': 0, 'PW': 0, 'PY': 0, 'QA': 0, 'QC': 0, 'QD': 0, 'QE': 0, 'QF': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QK': 0, 'QL': 0, 'QM': 0, 'QN': 0, 'QP': 0, 'QQ': 0, 'QR': 0, 'QS': 0, 'QT': 0, 'QV': 0, 'QW': 0, 'QY': 0, 'RA': 0, 'RC': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RY': 0, 'SA': 0, 'SC': 0, 'SD': 0, 'SE': 0, 'SF': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SK': 0, 'SL': 0, 'SM': 0, 'SN': 0, 'SP': 0, 'SQ': 0, 'SR': 0, 'SS': 0, 'ST': 0, 'SV': 0, 'SW': 0, 'SY': 0, 'TA': 0, 'TC': 0, 'TD': 0, 'TE': 0, 'TF': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TK': 0, 'TL': 0, 'TM': 0, 'TN': 0, 'TP': 0, 'TQ': 0, 'TR': 0, 'TS': 0, 'TT': 0, 'TV': 0, 'TW': 0, 'TY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'VH': 0, 'VI': 0, 'VK': 0, 'VL': 0, 'VM': 0, 'VN': 0, 'VP': 0, 'VQ': 0, 'VR': 0, 'VS': 0, 'VT': 0, 'VV': 0, 'VW': 0, 'VY': 0, 'WA': 0, 'WC': 0, 'WD': 0, 'WE': 0, 'WF': 0, 'WG': 0, 'WH': 0, 'WI': 0, 'WK': 0, 'WL': 0, 'WM': 0, 'WN': 0, 'WP': 0, 'WQ': 0, 'WR': 0, 'WS': 0, 'WT': 0, 'WV': 0, 'WW': 0, 'WY': 0, 'YA': 0, 'YC': 0, 'YD': 0, 'YE': 0, 'YF': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YK': 0, 'YL': 0, 'YM': 0, 'YN': 0, 'YP': 0, 'YQ': 0, 'YR': 0, 'YS': 0, 'YT': 0, 'YV': 0, 'YW': 0, 'YY': 0}
       
        pdbs = dff['pdb'][i]
        ag_protein = dff['ag_prot'][i]
        ab_protein = dff['ab_prot'][i]
        
        for res in ag_protein:
            if res in amino_acids_ag.keys():
                amino_acids_ag[res] += 1
        
        for res_ab in ab_protein:
            if res_ab in amino_acids_ab.keys():
                amino_acids_ab[res_ab] += 1
        
        for pair in pairs[i]:
            if pair in mains.keys():
                mains[pair] += 1
    
        # print(mains)
        for key,value in mains.items():
            k1 = key[0]
            k2 = key[1]
            nij = value
            ni = amino_acids_ab[k1]
            nj = amino_acids_ag[k2]
            if ni != 0 and nj != 0:
                pp = (nij / (ni*nj))*100
                pair_pref[key] = pp
            else:
                pp = 0
                pair_pref[key] = pp
        
        for key, value in pair_pref.items():
            if key in dict2_ag:
                dict2_ag[key].append(value)
            else:
                dict2_ag[key] = [value]
        
  
    averages_ag = {key: np.mean(values) for key, values in dict2_ag.items()}
    std_deviations_ag = {key: np.std(values) for key, values in dict2_ag.items()}   
    # convert dictionary to df
    # Convert dictionaries to DataFrame
    df1 = pd.DataFrame({'aa': list(averages_ag.keys()), 
                        'value': list(averages_ag.values()), 
                        'sd': list(std_deviations_ag.values())})
    
    
    df1_sorted = df1.sort_values(by='value', ascending=False)
    # write dataframes to excel file in different sheets
    # create excel writer
    os.chdir('/home/divya/Documents/ab_analysis/final_files/'+cog_pdb+'/')
    writer = pd.ExcelWriter('pair_preference_'+cog_pdb+'ab_final.xlsx')
    # write dataframe to excel sheet named 'marks'
    df1_sorted.to_excel(writer, 'pair_pref',index=False)
    # save the excel file
    writer.save()
                 
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            