import os
import ase
import numpy as np
import re
from ase.io import *
from ase.visualize import view
from itertools import repeat
import plot
import pickle
from scipy.signal import argrelextrema


def save_molecules(cluster):
    if not os.path.exists('./clusters'):
        os.mkdir('./clusters')
    with open('./clusters/tmp.pickle', 'wb') as f:
        pickle.dump(cluster , f)
def load_molecules():
    if os.path.exists('./clusters/tmp.pickle'):
        with open('./clusters/tmp.pickle', 'rb') as f:
            c = pickle.load(f)
        return c
    return None

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
class Cluster(ase.Atoms):
    
    def __init__(self , Atoms = None , G_solv = None , Setup = None):
        if(Atoms != None):
            ase.Atoms.__init__(self , Atoms)
            self.SolvShell = []
            self.E_cluster = list(self.info.keys())[0]
            self.size = len(self.numbers)
            self.mol_name = []
            self.mol_num = {}
            self.solvent_R = []
            self.molecules = self.get_molecules()
            self.spec_solu_solv()
            self.resolution = 20
            
        if(Setup != None):
            self.setup = Setup
                        
        if(G_solv != None):
            self.G_solv = G_solv + (-8.3145)/1000.0*300.0/4.18*np.log(55.34/self.setup["n_solv"])
            
    def spec_solu_solv(self):
        self.molecules = {'solute' if len(v) < 5 else 'solvent' : v for k,v in self.molecules.items()}
        self.mol_num = {'solute' if v < 5 else 'solvent' : v for k,v in self.mol_num.items()}
        self.tot_mol = (self.mol_num['solute']+self.mol_num['solvent'])
    @staticmethod
    def makedir(path):
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

    @staticmethod
    def view_ListofMol(mol , p = False):
        tmp = ase.Atoms()
        for i in mol:
           tmp.extend(i)
        view(tmp)
        if p:
            os.system('pause')
    @staticmethod
    def count_NumofAtom(mol):
        tmp = ase.Atoms()
        for i in mol:
           tmp.extend(i)
        return len(tmp.numbers)

        
            
    @staticmethod
    def merge_ListofMol(mol , info = None):
        tmp = ase.Atoms()
        if info != None:
            tmp.info = info
        for i in mol:
           tmp.extend(i)
        return tmp            
            
    @staticmethod
    def max_value(inputlist):
        return max([max(sublist) for sublist in inputlist])
    
    def set_resolution(self , resolution):
        self.resolution = resolution
    
    def get_rNc(self ,rs , resolution = None):
        if resolution == None:
            resolution = self.resolution

        test = np.linspace(np.min(rs)-1.0,np.max(rs)+1.0,resolution)
        dr = test[1] - test[0]
        print(dr)
        n = []
        for t in test:
            count = 0
            for r in rs:
                if r>=t and r<t+dr:
                    count+=1
            n.append((t+dr,count))
        return n
        
    def get_rNg(self,rs,resolution = None):
        if resolution == None:
            resolution = self.resolution
        test = np.linspace(np.min(rs)-1.0,np.max(rs)+1.0,resolution)
        dr = test[1] - test[0]
        rho_bulk = (self.tot_mol-1)/(4.0/3.0*np.pi*(np.max(rs)**3.0))
        print(dr)
        rNg=[]
        for t in test:
            count = 0
            for r in rs:
                if r>=t and r<t+dr:
                    count+=1
            rNg.append((t+dr,count/(4.0*np.pi*dr*(t**2.0))/rho_bulk))
            #rNg.append((t+dr,(4.0*np.pi*dr*(t**2.0))/rho_bulk))
        return rNg
    
    # def get_rNg(self , rs):
    #     nt = self.tot_mol
    #     return [(r,(g+1)/(r**2)/nt*(rs[-1]**2)) for g,r in enumerate(rs)]
    
    
    @staticmethod
    def get_rNdiffc(rNc):
        rNdiffc = [(0,0)]
        pre_c = 0
        for i , (r,c) in enumerate(rNc):
            if i == 0:
                pre_c = c
                continue
            diffc = c - pre_c
            pre_c = c
            rNdiffc.append((r,diffc))
        return rNdiffc
    
    def plot(self , d1 , d2 , ax , fig , pic):
        lines = plot.plot_errorbar2([d1,d2] , *ax)
        fig.set_figwidth(15)
        fig.set_figheight(9)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plot.plt.suptitle(f"Radial distribution of cluster" , fontsize=25)
        if pic :
            plot.plt.show()
        return lines

    def radial_distribution(self , sbin , pic = False , picPATH = None):
        mol_solu = self.molecules['solute'].copy()[0]
        mol_solv = self.molecules['solvent'].copy()
        p_center = mol_solu.get_center_of_mass()
        rs = []
        rNmol = []
        for i,mol in enumerate(mol_solv):
            p = mol.get_center_of_mass()
            r = np.linalg.norm(p-p_center)
            rs.append(r)
            rNmol.append((r,i))
        rs = sorted(rs)
        rNmol = sorted(rNmol)
        
        rNc = self.get_rNc(rs , sbin)
        rNg = self.get_rNg(rs , sbin)
        # rNd = [(r,1) for g,r in enumerate(rs)]
        # rNdiffc = self.get_rNdiffc(rNc)
        
        # plot.plt.scatter(*zip(*rNd))
        if pic or picPATH != None:
            fig , ax = plot.plt.subplots(2,1)
            self.plot(rNg , rNc , ax , fig , pic)
            self.makedir(picPATH)
            plot.plt.savefig(picPATH)
            plot.plt.close()
        return (rNmol , rNg , rNc)
        
    def sort_molsByr(self , rNmol):
        self.solvent_R = [i[0] for i in  rNmol]
        self.molecules['solvent'] = [self.molecules['solvent'][i[1]] for i in  rNmol]
        
    def find_Middleof3maxR(self , rNg , rNc , Min_1stShell = 10):
        """find 2 radius to separte the cluster into 3 parts
        Radius = middle points between the 2 highest peaks
        
        the highest peak constraints:
        1. local maximum (if there are only 2 values,the third value will choose the final point)
        2. highest value in local maximum
        2. the first Rm should include Min_1stShell's molecules
        
        Args:
            rNg (list of tuples): _description_
            rNc (list of tuples): _description_
            Min_1stShell (int, optional): The minimum molecules include in 1st solvation shell. Defaults to 10.

        Returns:
            Tuple: (R1 , R2)
        """        
        
        tmp3 = 0.0
        rNt = []
        rmax = max([r for r,c in rNc])
        for r,c in rNc:
            tmp3 += c
            rNt.append(tmp3)
        tmp2 = np.array([g for r,g in rNg])
        tmp2 = argrelextrema(tmp2 , np.greater)
        tmp = [rNg[i] for i in tmp2[0] if i > 1 and rNt[i+1] > Min_1stShell]
        tmp = sorted([(v , i , r) for i, (r , v) in enumerate(tmp)])
        tmp = sorted([(r,i) for v,i,r in tmp[-1:-6:-1]])
        if len(tmp) == 2:
            Rm = [(tmp[0][0]+tmp[1][0])/2 , (tmp[1][0]+rmax)/2]
        else :
            Rm = [(tmp[0][0]+tmp[1][0])/2 , (tmp[1][0]+tmp[2][0])/2]
            
        return Rm
            
    def writefile(self , Rm , path , type , p=False , picPATH = None):
        start = 0
        end = 0
        info = {'SS' : ''}
        Erow = 0
        Srow = 3
        start_end = []
        for j,R_set in enumerate(Rm):
            if j == 0:
                Erow = Srow + self.count_NumofAtom(self.molecules['solute'])
            for i,r in enumerate(self.solvent_R , start = start):
                if r <= R_set:
                    end = i
                else :
                    break
            self.SolvShell.append(self.molecules['solvent'][start:end + 1:1])
            
            Erow += self.count_NumofAtom(self.SolvShell[j])
            info['SS'] += f' {j}/{Srow}:{Erow} '            
            start_end.append((Srow,Erow))
            Srow = Erow
            start = end + 1
        j += 1
        self.SolvShell.append(self.molecules['solvent'][start::])
        Erow += self.count_NumofAtom(self.SolvShell[j])
        info['SS'] += f' {j}/{Srow}:{Erow} '
        start_end.append((Srow,Erow))

        final_struc = self.merge_ListofMol(self.molecules['solute'] + self.SolvShell[0]+self.SolvShell[1]+self.SolvShell[2] , info = info)
        self.makedir(path)
        
        # write(path , final_struc , 'proteindatabank' )
        write(path , final_struc , type )

        if p:
            self.view_ListofMol(self.SolvShell[0]+self.molecules['solute'])
            self.view_ListofMol(self.SolvShell[0]+self.molecules['solute']+self.SolvShell[1])
            self.view_ListofMol(self.SolvShell[0]+self.molecules['solute']+self.SolvShell[1]+self.SolvShell[2] , True)
        if picPATH != None:
            self.makedir(picPATH)
            write(picPATH + f'SS1.png' , self.merge_ListofMol(self.molecules['solute'] + self.SolvShell[0]))
            write(picPATH + f'SS2.png' , self.merge_ListofMol(self.molecules['solute'] + self.SolvShell[0]+self.SolvShell[1]))
            write(picPATH + f'SS3.png' , self.merge_ListofMol(self.molecules['solute'] + self.SolvShell[0]+self.SolvShell[1]+self.SolvShell[2]))
            
        return start_end
    
    def gen_length_table(self):
        length = [[[]for i in repeat(None, 40)] for i in repeat(None, 40)]
        length[0][0] =  0.75    ;                                           
        length[0][5] =  1.09    ;   length[5][5] = [1.54 , 1.34 , 1.20]     
        length[0][6] =  1.01    ;   length[5][6] = [1.43 , 1.38 , 1.16]     ;   length[6][6] = [1.45 , 1.25 , 1.10]     ;   
        length[0][7] =  0.96    ;   length[5][7] = [1.43 , 1.20 , 1.13]     ;   length[6][7] = [1.36 , 1.22]            ;   length[7][7] = [1.48 , 1.21]
        length[0][8] =  0.92    ;   length[5][8] =  1.35                    ;   length[6][8] = 1.36                     ;   length[7][8] = 1.42
        length[0][16] = 1.27    ;   length[5][16] = 1.77                    ;   length[6][16] = 1.75                    ;   length[7][16] = 1.72
        length[0][34] = 1.27    ;   length[5][34] = 1.94                    ;   length[6][34] = 1.89                    ;   length[7][34] = 2.06
        
        self.length = length
        
    def get_length_range(self , i , j) :
        i = min([i,j])
        j = max([i,j])
        if not hasattr(self , 'length'):
            self.gen_length_table()
        tmp1 = self.length[i-1][j-1]
        
        if isinstance(tmp1 , list):
            return (max(tmp1) + 0.5, min(tmp1) - 0.5)
        else:
            return (tmp1 + 0.5 , tmp1 - 0.5)
        
    def get_molecules(self):
        molecules = {}
        self.distm = self.get_all_distances()
        distm = self.distm.copy()
        chkm = []
        for i in range(self.size):
            if i not in chkm:
                mol,chk = self.find_molecule(i , distm , [] , [] )
                chkm.extend(chk)
                s = str()
                for n in mol:
                    s += n[1]
                mol = ase.Atoms(symbols = s , positions = [self.positions[n[0]] for n in mol])
                self.mol_name.append(s)
                
                if molecules.get(s) == None:
                    molecules[s] = [mol]
                else:
                    molecules[s].append(mol)
                    
                if self.mol_num.get(s) == None:
                    self.mol_num[s] = 1
                else :
                    self.mol_num[s] += 1
                    
        if len(chkm) != self.size:
            e = self.size - len(chkm)
            print(bcolors.WARNING + f'size mismatch for molecules : {e} atom(s) missing' + bcolors.ENDC)
        return molecules
        
    def find_molecule( self , i , distm ,  molecule = [] , chkm = [] ) :
        for j,d in enumerate(distm[i]):
            if j in chkm :
                continue
            if j == i :
                chkm.append(j)
                molecule.append((j,self.symbols[j]))
                
            max_l , min_l = self.get_length_range(self.numbers[i] ,self.numbers[j])
            if d < max_l and d > min_l:
                chkm.append(j)
                if (j,self.symbols[j]) not in molecule:
                    molecule.append((j,self.symbols[j]))
                if self.symbols[j] != 'H':
                    molecule,chkm = self.find_molecule( j , distm , molecule = molecule  , chkm = chkm )
                    # molecule.extend(tmp1)
                    # chkm.extend(tmp2)     
        molecule.sort(key=lambda tup: tup[0])
        chkm.sort()
        return (molecule,chkm)
            
    def get_E_cluster(self) :
        return self.E_cluster
    
    def set_Gsolv(self , G_solv):
        self.G_solv = G_solv

    def read_G_solv_file(self , path) :
        clusters = []
        max_run = 0
        with open( path , "r" ) as f:
            mol_dataset = f.readlines()
        for data in mol_dataset :
            data = data.split("|")
            
            setup = self.get_setup(data)
#            G_solv = float(re.findall( r"[-+]?(?:\d*\.\d+|\d+)" , data[-2])[0])
            G_solv = float(re.findall( r"[-+]?(?:\d*\.\d+|\d+)" , data[-2])[0])
            
            if setup["run"] > max_run :
                max_run = setup["run"]
            NUM_RUN = max_run
            
            clusters.append(Cluster(G_solv=G_solv,Setup = setup))
        return clusters
            
    def get_setup(self,data):
        setup = dict()
        
        path = data[0].split("/")
        
        strct_name = path[0]
        setting = path[1]
        run = int(re.findall( r"\d+" , path[2] )[0])
        n_solv = int(re.findall( r"\d+" , path[3] )[0] )
        for i in ( "strct_name" , "setting" , "run" , "n_solv" ):
            setup[i] = locals()[i]
        return setup

    def read_trajectory(self , path , num = ':'):
        tmp = read(path,num)
        
        if isinstance(tmp , list):
            return [Cluster(cluster) for cluster in tmp]
        else:
            c = load_molecules()
            if c != None:
                return c
            return Cluster(tmp)
    
    def read_cluster(self , path):
        return Cluster(read(path))

def pdb_resitype_modified(path , start_end , solute_sz):
    fmt = '{:>4}  {:>5} {:<4} {:>3}     {:>1}    {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}'
    with open(path) as f:
        data = f.read().splitlines()
    start,end = start_end[0]
    for i,line in enumerate(data[start-2:end-2],start):
        data[i-2] = fmt.format(*(line.split()[:3] + ['FIR'] + line.split()[4:]))
        
    for i,line in enumerate(data[1:solute_sz+1],1):
        data[i] = fmt.format(*(line.split()[:3] + ['SOU'] + line.split()[4:]))
        
    start,end = start_end[1]
    for i,line in enumerate(data[start-2:end-2],start):
        data[i-2] = fmt.format(*(line.split()[:3] + ['SEC'] + line.split()[4:]))
        
    start,end = start_end[2]
    for i,line in enumerate(data[start-2:end-2],start):
        data[i-2] = fmt.format(*(line.split()[:3] + ['THR'] + line.split()[4:]))
    data.append('')
    
    with open(path , 'w') as f:
        f.write('\n'.join(data))

if __name__ == '__main__':
    
    # c = Cluster.read_trajectory(Cluster(),"plot_py/final_ensemble.xyz" , num = 9)
    # qcg_final_ensemble_inPath = "./qcg/clus_sz_ck/MEA/default/run{1}/s{0}/ensemble/final_ensemble.xyz"
    # rNg_plot_outPath = './clusters/result/solvation shell/MEA/default/run{1}_s{0}_{2}_rNg.png'
    # struture_outPath = './clusters/result/structure/MEA/default/run{1}/s{0}/solute_{2}.xyz'
    # SS_pic_outPath = './clusters/result/solvation shell/MEA/default/run{1}_s{0}_{2}_'
    
    
    # qcg_final_ensemble_inPath = "./qcg/clus_sz_ck/p-MEA/default/run{1}/s{0}p/ensemble/final_ensemble.xyz"
    # rNg_plot_outPath = './clusters/result/solvation shell/p-MEA/default/run{1}_s{0}_{2}_rNg.png'
    # struture_outPath = './clusters/result/structure/p-MEA/default/run{1}/s{0}/solute_{2}.pdb'
    # SS_pic_outPath = './clusters/result/solvation shell/p-MEA/default/run{1}_s{0}_{2}_'
    
    qcg_final_ensemble_inPath = "./qcg/clus_sz_ck/MEA/default/run{1}/s{0}/ensemble/final_ensemble.xyz"
    rNg_plot_outPath = './clusters/result/solvation shell/MEA/default/run{1}_s{0}_{2}_rNg.png'
    struture_outPath = './clusters/result/structure/MEA/default/run{1}/s{0}/solute_{2}.pdb'
    SS_pic_outPath = './clusters/result/solvation shell/MEA/default/run{1}_s{0}_{2}_'
    
    
    for k in range(0,20):
        for j in [1,2,3,4,5]:
            for i in [100,125,150,175,200,250]:
            # for i in [200,250]:
                # c = Cluster.read_trajectory(Cluster(),f"./qcg/clus_sz_ck/p-MEA/default/run{j}/s{i}p/ensemble/final_ensemble.xyz" , num = k)
                c = Cluster.read_trajectory(Cluster(),qcg_final_ensemble_inPath.format(i,j) , num = k)
                # view(read("plot_py/final_ensemble.xyz",':')[0])
                # os.system("pause")
                print(bcolors.HEADER + str(i) + bcolors.ENDC)
                c.set_resolution(20)
                rNmol , rNg , rNc = c.radial_distribution(20 , False , rNg_plot_outPath.format(i,j,k))  
                rNmol , rNg , rNc = c.radial_distribution(20 , False )  
                c.sort_molsByr(rNmol)
                
                rNt = []
                Rm = c.find_Middleof3maxR(rNg , rNc , 13)
                print(Rm)
                #start_end = c.writefile(Rm , struture_outPath.format(i,j,k) , 'extxyz' )
                start_end = c.writefile(Rm , struture_outPath.format(i,j,k) , 'proteindatabank' )
                pdb_resitype_modified(struture_outPath.format(i,j,k) , start_end , 12)