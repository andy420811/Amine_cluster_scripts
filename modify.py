from constants import *

def get_start_end(raw_data):
    raw_data[1].split()
    for line in raw_data:
        if 'SS' in line:
            tmp = line.split('"')[1].split()

    start = []
    end = []
    for i in tmp:
        start.append(i[2:].split(":")[0])
        end.append(i[2:].split(":")[1])
    
    return start,end
        
def add_ESP_charges(raw_data , solute_ESP_charges):
    raw_data[1] = ''
    for i,line in enumerate(raw_data[2:]):
        if i < len(solute_ESP_charges):
            raw_data[i+2] = raw_data[i+2][0] + '-' + str(solute_ESP_charges[i]) + '\t' + raw_data[i+2][1:]
        else:
            if raw_data[i+2][0] == 'O':
                raw_data[i+2] = raw_data[i+2][0] + '-' + str(H2O_ESP_charges[0]) + '\t' + raw_data[i+2][1:]
            if raw_data[i+2][0] == 'H':
                raw_data[i+2] = raw_data[i+2][0] + '-' + str(H2O_ESP_charges[1]) + '\t' + raw_data[i+2][1:]

def add_Atom_types(raw_data , solute_Atom_types):
    for i,line in enumerate(raw_data[2:]):
        if i < len(solute_Atom_types):
            raw_data[i+2] = raw_data[i+2][0] + '-' + str(solute_Atom_types[i]) + raw_data[i+2][1:]
        else:
            if raw_data[i+2][0] == 'O':
                raw_data[i+2] = raw_data[i+2][0] + '-' + str(H2O_Atom_types[0]) + raw_data[i+2][1:]
            if raw_data[i+2][0] == 'H':
                raw_data[i+2] = raw_data[i+2][0] + '-' + str(H2O_Atom_types[1]) + raw_data[i+2][1:]


def add_oniom_level(raw_data , start , end):
    s = int(start[0]) - 1
    e = int(end[0]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = line + ' H'
    s = int(start[1]) - 1
    e = int(end[1]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = line + ' M'
    s = int(start[2]) - 1
    e = int(end[2]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = line + ' L'
    return raw_data
        
def add_fixed_atoms(raw_data , start , end , fixed_atoms = (0 , -1 , -1)):
    s = int(start[0]) - 1
    e = int(end[0]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = '\t'.join([line.split()[0]] + [' ' + str(fixed_atoms[0]) + ' '] + line.split()[1:])
    s = int(start[1]) - 1
    e = int(end[1]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = '\t'.join([line.split()[0]] + [' ' + str(fixed_atoms[1]) + ' '] + line.split()[1:])
    s = int(start[2]) - 1
    e = int(end[2]) - 1
    for i,line in enumerate(raw_data[s:e] , s):
        raw_data[i] = '\t'.join([line.split()[0]] + [' ' + str(fixed_atoms[2]) + ' '] + line.split()[1:])
    return raw_data
    
        
def add_header(raw_data , command , title , Chag_Mult):
    tmp = []
    tmp += command
    tmp.append('')
    tmp.append(title)
    tmp.append('')
    tmp.append(Chag_Mult)
    raw_data = tmp + raw_data[2:]
    
    return raw_data

def add_vdw_para(raw_data):
    raw_data.append('')
    raw_data.append('VDW HW 0.6000008633 0.01571502392')
    raw_data.append('VDW NT 1.8240008285 0.17016267942')
    raw_data.append('VDW HO 0.6000008633 0.01571502392')
    
    return raw_data

input_path = "GAUSSIAN/oniom/MEA/test2/MEA_out.xyz"
# input_path = "./clusters/result/structure/MEA/default/run4/s250/solute_0.xyz"
command = [
    '%mem=190GB' , 
    # '#p amber=print oniom(b3lyp/6-31g(d,p):pm7:amber=(Print,SoftFirst)) opt',
    '#p oniom(b3lyp/6-31g(d,p):pm7:amber=(Print,SoftFirst)) opt',
    ]
title = 'MEA oniom'
Chag_Mult = '0 1 0 1 0 1 0 1 0 1 0 1'


data = []
with open(input_path) as f :
    raw_data = f.read().splitlines()

start , end = get_start_end(raw_data)

esp_data = add_ESP_charges(raw_data , MEA_ESP_charges)

types_data = add_Atom_types(raw_data , MEA_Atom_types)

oniom_data = add_oniom_level(raw_data , start, end)

fixed_atoms_data = add_fixed_atoms(raw_data , start, end)

data = add_header(raw_data , command ,title , Chag_Mult)

data = add_vdw_para(data)

with open("./clusters/result/gaussian input/MEA_3.com" , 'w') as f:
    data += ['\n']
    f.write('\n'.join(data))
