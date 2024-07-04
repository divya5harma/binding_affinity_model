#!/usr/bin/env python

import math
import pandas as pd

def hydro(file, chname1):
    print(file)
    atom = []
    hy_index = {
        "ala": 0.87,
        "asp": 0.66,
        "cys": 1.52,
        "glu": 0.67,
        "phe": 2.87,
        "gly": 0.10,
        "his": 0.87,
        "ile": 3.15,
        "lys": 1.64,
        "leu": 2.17,
        "met": 1.67,
        "asn": 0.09,
        "pro": 2.77,
        "gln": 0.00,
        "arg": 0.85,
        "ser": 0.07,
        "thr": 0.07,
        "val": 1.87,
        "trp": 3.77,
        "tyr": 2.67
    }
    count = []
    residue = []
    res_no = []
    chain_e = []
    shydro = []
    # file3 = open("sh.csv","w")
    # file3.write("res1"+","+"resno"+","+"chain"+","+"sh"+","+"\n") 
    with open(file, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") and "CA" in line:
                ch = line[21]
                if chname1 == "na":
                    atom.append(line)
                elif ch == chname1:
                    atom.append(line)
            elif line.startswith("ENDMDL"):
                break

        
        for line in atom:
            c = hydro = 0
            ares1 = line[17:20]  # residue name
            # print(ares1)
            ano1 = line[22:27]  # atom number
            achain = line[21]
            x1, y1, z1 = float(line[30:38]), float(line[38:46]), float(line[46:54])
            
            for j in atom:
                ares2 = j[17:20]  # residue name
                ano2 = j[22:27]  # atom number
                x2, y2, z2 = float(j[30:38]), float(j[38:46]), float(j[46:54])
                dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

                if dist > 0 and dist <= 8:
                    st = ares2.lower()
                    val = hy_index.get(st, 0)
                    hydro += val 
            shydro.append(hydro)
            residue.append(ares1)
            res_no.append(ano1)
            chain_e.append(achain)

    df1 = pd.DataFrame({'res':residue,'resno':res_no,'chain':chain_e,'sh':shydro})  
    # print(df1)
    # df1['resno'] = df1['resno'].astype(str)
            # file3.write(str(ares1)+","+str(ano1)+","+str(achain)+","+str(hydro)+"\n")
    # print(df1['resno'].to_list())
    avg_sh = []
    pdb = file.split('.pdb')[0].split('pdbs_all_final/')[1]
    df2 = pd.read_csv('res_numbers.csv')
    for i in range(len(df2)):
        epi = df2['l_resno'][i]
        if df2['pdb'][i] == pdb:
            if isinstance(epi, float):
                avg_sh.append(0)
            # print(epi)
            elif isinstance(epi, str):
                ep = epi.split(',')
                e = [item.strip() for item in ep]
                print(e)
                for j in range(len(df1)):
                    resno = df1['resno'][j].strip()
                    # resno = int(resno)
                    sh = df1['sh'][j]
                    if resno in e:
                        # print(resno)
                        avg_sh.append(sh)
    print(avg_sh)            
    mean_sh = sum(avg_sh) / len(avg_sh)
    # tot_sh = sum(avg_sh)
    # file3.close()
    return mean_sh
dfz = pd.read_csv('res_numbers.csv')
result = dfz['pdb'].tolist()
name = []
final_sh = []
for i in result:
    x = hydro('/home/divya/Documents/ab_analysis/pdbs_all_final/'+i+'.pdb','L')
    name.append(i)
    final_sh.append(x)
    
dff = pd.DataFrame(
    {'name': name,
     'sh': final_sh
    }) 
print(dff.head())
dff.to_csv('final_sh_para_l.csv')
# dff.to_clipboard(excel=True,index=False)    
    
