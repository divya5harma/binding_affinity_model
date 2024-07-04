#Propensity for antibody

############# Binding propensity ##########
import pandas as pd
import os
from Bio import PDB
import numpy as np
os.chdir('/home/divya/Documents/ab_analysis/reanalysis_abs/')
pdb_path = '/home/divya/Documents/ab_analysis/reanalysis_abs/interactions_all_abs/'
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
nopdb = []

## get epitope and paratope sequence
for x in range(len(dfz)):
    cog_pdb = dfz['cognate_pdb'][x].lower()
    f_ab = dfz['ab_pdb'][x]
    f_aa = f_ab.split(',')
    print(f_aa)
    p = []
    e = []
    pdb = []
    ab_prot = []
    epi_prot = []
    for i in f_aa:
        f = i.lower()
        res_para_seq = ''
        res_epi_seq = ''
        file = open(pdb_path+'all_interaction_'+f+'.pdb.csv','r')
        file_detail = file.readlines()
        for line in file_detail:
            l = line.split(',')
            l_clean = [value.strip('"') for value in l if value.strip('"') and value != '\n']
            if l_clean:
                ## epitope seq
                if l_clean[0] == 'Epitope residues':
                    # print('E: ',l_clean)
                    res = l_clean[1:]
                    # Filter out empty strings and newline characters, and remove the double quotes
                    # res_clean = [value.strip('"') for value in res if value.strip('"') and value != '\n']
                    res = [s.strip('\n') for s in res]
                    for i in res:
                        res_epi_seq += i[-1]
                
                elif l_clean[0] == 'Paratope residues':
                    # print('P: ',l_clean)
                    resp = l_clean[1:]
                    resp = [s.strip('"\n') for s in resp]
                    # Filter out empty strings and newline characters, and remove the double quotes
                    # resp_clean = [value.strip('"') for value in resp if value.strip('"') and value != '\n']
                    for i in resp:
                        res_para_seq += i[-1]  
        pdb.append(f)
        e.append(res_epi_seq)
        p.append(res_para_seq)


    ######## get sequence of protein from the pdb file
    
    
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

            
####### get the propensity with standard deviation
    dff = pd.DataFrame({'pdb':pdb,'para_interface':p,'epi_interface':e,'ab_prot':ab_prot,'ag_prot':epi_prot})
   
    dict2_ag = {}
    dict2_ab = {}
    for i in range(len(dff)):
        try:
            amino_acids_interface = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
            amino_acids = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
            len_interface = 0
            len_protein = 0
            amino_acids_ab_interface = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
            amino_acids_ab = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
            len_interface_ab = 0
            len_protein_ab = 0
            
            pdbs = dff['pdb'][i]
            epi_interface = dff['epi_interface'][i]
            para_interface = dff['para_interface'][i]
            ag_protein = dff['ag_prot'][i]
            ab_protein = dff['ab_prot'][i]
            
            ## propensity of antigen
            len_interface = len_interface + len(epi_interface)
            len_protein = len_protein + len(ag_protein)
            for aa_res in epi_interface:
                if aa_res in amino_acids_interface.keys():
                    amino_acids_interface[aa_res] += 1
            for res in ag_protein:
                if res in amino_acids.keys():
                    amino_acids[res] += 1
            
            ## propensity of antibody
            len_interface_ab = len_interface + len(para_interface)
            len_protein_ab = len_protein + len(ab_protein)
            for aa_res_ab in para_interface:
                if aa_res_ab in amino_acids_ab_interface.keys():
                    amino_acids_ab_interface[aa_res_ab] += 1
            for res_ab in ab_protein:
                if res_ab in amino_acids_ab.keys():
                    amino_acids_ab[res_ab] += 1
            
            # Divide values of corresponding keys
            num_ag = {key: amino_acids_interface[key] / amino_acids[key] if amino_acids[key] != 0 else 0 for key in amino_acids_interface}
            denom_ag = len_interface/len_protein
            aa_propensity_ag = {key: value / denom_ag for key, value in num_ag.items()}
            
            # Divide values of corresponding keys
            num_ab = {key: amino_acids_ab_interface[key] / amino_acids_ab[key] if amino_acids_ab[key] != 0 else 0 for key in amino_acids_ab_interface}
            denom_ab = len_interface_ab/len_protein_ab
            aa_propensity_ab = {key: value / denom_ab for key, value in num_ab.items()}
        
            for key, value in aa_propensity_ag.items():
                if key in dict2_ag:
                    dict2_ag[key].append(value)
                else:
                    dict2_ag[key] = [value]
            for key, value in aa_propensity_ab.items():
                if key in dict2_ab:
                    dict2_ab[key].append(value)
                else:
                    dict2_ab[key] = [value]
        except:
            nopdb.append(pdbs)
        # Calculate averages and standard deviations
    averages_ag = {key: np.mean(values) for key, values in dict2_ag.items()}
    std_deviations_ag = {key: np.std(values) for key, values in dict2_ag.items()}
    averages_ab = {key: np.mean(values) for key, values in dict2_ab.items()}
    std_deviations_ab = {key: np.std(values) for key, values in dict2_ab.items()}
        # convert dictionary to df
        
    # Convert dictionaries to DataFrame
    df1 = pd.DataFrame({'aa': list(averages_ag.keys()), 
                        'value': list(averages_ag.values()), 
                        'sd': list(std_deviations_ag.values())})
        
   
    df2 = pd.DataFrame({'aa': list(averages_ab.keys()), 
                        'value': list(averages_ab.values()), 
                        'sd': list(std_deviations_ab.values())})
         
    
    # create excel writer
    os.chdir('/home/divya/Documents/ab_analysis/reanalysis_abs/'+cog_pdb+'/')
    writer = pd.ExcelWriter('propensity_'+cog_pdb+'_ab.xlsx')
    # write dataframe to excel sheet named 'marks'
    df1.to_excel(writer, 'antigen')
    df2.to_excel(writer,'antibody')
    # save the excel file
    writer.save()         