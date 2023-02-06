from matplotlib import pyplot as plt
from matplotlib import font_manager,figure
import sys
from statistics import stdev,mean
import numpy as np
import re
import ase
from ase.io import read,write
from scipy.optimize import curve_fit


#font_manager.findSystemFonts(fontpaths=None, fontext="ttf")
font_manager.findSystemFonts(fontpaths=None)
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

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

NUM_RUN = 5
class crestout():
    def __init__(self , path):
        self.read_crest(path)

    def read_crest(self , path):
        with open(path , "r") as f:
            lines = f.readlines()
        tmp = []
        for i ,line in enumerate(lines):
            if "Command line input:" in line:
                tmp.append(i)
            if "input parameters" in line:
                tmp.append(i)
            if "Total Energy of solute:" in line:
                tmp.append(i)
            if "Total energy of solvent:" in line:
                tmp.append(i)
                
class Cluster(ase.Atoms):
    def __init__(self , Atoms = None , G_solv = None , Setup = None):
        if(Atoms != None):
            self.arrays = Atoms.arrays
            self.numbers = Atoms.numbers
            self.positions = Atoms.positions
            self.symbols = Atoms.symbols
            self.info = Atoms.info
            self.ase_objtype = Atoms.ase_objtype
            self.calc = Atoms.calc
            self.E_cluster = list(self.info.keys())[0]

        if(Setup != None):
            self.setup = Setup
                        
        if(G_solv != None):
            self.G_solv = G_solv + (-8.3145)/1000.0*300.0/4.18*np.log(55.34/self.setup["n_solv"])
            
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


def read_clusters( ensemble_trajection_path) :
    return [Cluster(cluster) for cluster in read(ensemble_trajection_path,':')]


def plot_errorbar2(input , ax_line = None , ax_bar = None):
    lines = []
    barcolor = ['#DDA15E' , '#69140E' , '#9b591c' , '#69140E' , '#69140E']
    
    data1 = input[0]
    data2 = input[1]
    # if no input ax , it will be a error line drawn in default
    if ax_line and ax_bar == None:
        ax_line =  plt.subplot()
    
    if ax_line != None :
        lines.append(ax_line.scatter(*zip(*data1) ,  color = 'black' , alpha = 0.2 , s = 15))
        lines.append(*ax_line.plot(*zip(*data1),lw=1.5  , color = '#283618' , alpha = 1 ))
        ax_line.set_xlabel('r (Å)',fontsize=20)
        ax_line.set_ylabel('g(r)',fontsize=20)
        ax_line.tick_params(axis='both', which='major', labelsize=15)
        ax_line.set_xlim((0 , data1[-1][0] + 1))
        ax_line.legend(loc='best',frameon=False)
        
    if ax_bar != None:
        
        ax_bar.bar(*zip(*data2),color=barcolor[0] , alpha=0.5 , label = f'Increase molecule' ,width= 0.1)        
        plt.axhline(y=0,linewidth=1, color='black')
        ax_bar.set_ylabel('molecule (n)' , fontsize=20)
        ax_bar.tick_params('y')
        ax_bar.set_xlabel('r (Å)',fontsize=20)
        ax_bar.tick_params(axis='both', which='major', labelsize=15)
        ax_bar.set_xlim((0 , data1[-1][0] + 1))
        ax_bar.legend(loc='best',frameon=False)
        
    return lines

def plot_errorbar(input , ax_line = None , ax_bar = None , n = [-sys.maxsize-1 , sys.maxsize]):
    data = []
    lines = []
    y_stdev = []
    y_avg = []
    runs = []
    x = []
    barcolor = ['#DDA15E' , '#69140E' , '#9b591c' , '#69140E' , '#69140E']
    
    # if no input ax , it will be a error line drawn in default
    if ax_line and ax_bar == None:
        ax_line =  plt.subplot()
    if isinstance(input , dict):
        for keys,value in input.items() :
            if keys <= float(n[1]) and keys >= float(n[0]):
                runs.append(len(value))
                data.extend([(keys,v) for v in value])
                if len(value) != NUM_RUN:
                    print(bcolors.WARNING + f"WRONG number at {keys} : {len(value)-NUM_RUN}" + bcolors.ENDC)
                if len(value) == 1:
                    y_stdev.append(0)
                else: 
                    y_stdev.append(stdev(value))
                y_avg.append(mean(value))
                x.append(float(keys))
            
    
    if ax_line != None :
        lines.append(ax_line.scatter(*zip(*data) ,  color = 'black' , alpha = 0.2 , s = 15))
        if y_stdev != None:
            lines.append(ax_line.errorbar(x,y_avg , yerr = y_stdev , label = "stdev" , fmt = '.' , capsize = 5 , capthick = 1 , ecolor = '#606C38' ))
        lines.append(*ax_line.plot(x,y_avg,lw=1.5  , color = '#283618' , alpha = 1 ))
        ax_line.set_xlabel('number of solvent molecules',fontsize=20)
        ax_line.set_ylabel('G_solv (kcal/mol)',fontsize=20)
        ax_line.legend(loc='best',frameon=False)
    if ax_bar != None:
        
        runs = get_diffind(runs)
        x_plot = [x[runs[i][0]:runs[i+1][0]] for i in range(len(runs)-1)]
        y_stdev_plot = [y_stdev[runs[i][0]:runs[i+1][0]] for i in range(len(runs)-1)]
        
        for i , (j , k) in enumerate(zip(x_plot , y_stdev_plot)):
            ax_bar.bar(j,k,color=barcolor[runs[i][2]] , alpha=0.5 , label = f'stdev in {runs[i][1]} runs' ,width= 4)        
            
            
        ax_bar.set_ylabel('stdev (kcal/mol)' , fontsize=20)
        ax_bar.set_ylim(0,max(y_stdev)*1.5)
        ax_bar.tick_params('y')
        ax_bar.set_xlabel('number of solvent molecules',fontsize=20)
        ax_bar.legend(loc='best',frameon=False)
        
    return (x , y_avg , y_stdev , lines)

def get_diffind(i):
    index  = []
    map = []
    colormap = []
    j = 0
    value = 0
    prevalue = 0
    for j,value in enumerate(i):
        if prevalue != value:
            prevalue = value
            index.append(j)
            map.append(value)
            if value in map:
                colormap.append(map.index(value))
    index.append(j+1)
    map.append(value)
    colormap.append(map.index(value))
    return list(zip(index,map,colormap))
        
def read_cluster_solv(strct_name , setting , G_solv_path = "plot_py/G_solv.txt"):
    
    # G_solvs : a dictionary for a key number of solvent molecules to the G_solv list include n runs
    G_solvs = {}
    n_runs = {}
    for cluster in Cluster().read_G_solv_file( G_solv_path ):
        #if cluster.setup["setting"] == "default":
        
        # limit of the selecting data 
        if cluster.setup["strct_name"] == strct_name and cluster.setup["setting"] == setting :
            try:
                G_solvs[cluster.setup["n_solv"]].append(cluster.G_solv)
                n_runs["n_solv"] = len(G_solvs[cluster.setup["n_solv"]])
            except:
                G_solvs[cluster.setup["n_solv"]] = [cluster.G_solv]
            
    G_solvs = dict(sorted(G_solvs.items()))
    
    return (n_runs  , G_solvs)

def plot_fitting( x , y_avg , y_stdev , ax , lines) :
    
    def linearFunc(x,intercept,slope,n):
        y = intercept + slope * (np.array(x)**(n))
        return y
    
    # a_fit,cov=curve_fit(linearFunc,x,y_avg,sigma=y_stdev,absolute_sigma=True,maxfev = 5000 , bounds= ([-np.inf , -np.inf , -1.0 ],[ np.inf,np.inf, - 0.01]))
    if min(y_stdev) == 0:
        a_fit,cov=curve_fit(linearFunc,x,y_avg,maxfev = 30000 )
        print(bcolors.WARNING + "The fitting line do not consider the stddev" + bcolors.ENDC)
    else :
        #a_fit,cov=curve_fit(linearFunc,x,y_avg,sigma=y_stdev,absolute_sigma=True,maxfev = 30000 , p0= ( [ 0 , 1.0 , 0.1 ] ))
        a_fit,cov=curve_fit(linearFunc,x,y_avg,sigma=y_stdev,absolute_sigma=True,maxfev = 30000)
    inter = a_fit[0]
    slope = a_fit[1]
    n = a_fit[2]
    
    d_inter = np.sqrt(cov[0][0])
    d_slope = np.sqrt(cov[1][1])
    d_n = np.sqrt(cov[2][2])
    
    # Compute a best fit line from the fit intercept and slope.
    yfit = inter + slope*(np.array(x)**(n))
    # Create a graph of the fit to the data. We just use the ordinary plot
    # command for this.
    lines.append(*ax.plot(x,yfit,label = f'Fit : y = {slope:.2f} * (x ^ {n:.2f}) + {inter:.2f}\nUncertainty : slope : {d_slope:.2f}, intercept : {d_inter:.2f}, n : {d_n:.2f}' , color = "#BC6C25"))
    ax.set_xlabel('number of solvent molecules',fontsize=20)
    ax.set_ylabel('G$_{solv}$ (kcal/mol)',fontsize=20)
    ax.legend(loc='best',frameon=False)

    print(f'The slope = {slope}, with uncertainty {d_slope}')
    print(f'The intercept = {inter}, with uncertainty {d_inter}')
    print(f'The n = {n}, with uncertainty {d_n}')
    
    chisqr = sum((y_avg-linearFunc(x,inter,slope,n))**2/np.array(y_stdev)**2)
    dof = len(y_avg) - 2
    chisqr_red = chisqr/dof
    
    print(f'Reduced chi^2 = {chisqr_red}')

def plot_cluster(strct_name , setting , config = (2,1) , mode = (True , True ,True), n = [-sys.maxsize-1 , sys.maxsize] , path = "" , preData = []) :
    #mode[0] error line , mode[1] error bar , mode[2] fitting line
    data = preData
    n_runs , G_solvs = read_cluster_solv(strct_name , setting , path)
    data.append(G_solvs)

    fig , ax = plt.subplots(*config)
    if len(ax) == 1:
        x , y_avg , y_stdev ,lines = plot_errorbar(G_solvs , ax , n)
        if mode[2] :
            plot_fitting( x , y_avg , y_stdev , ax , lines)
    else:
        x , y_avg , y_stdev , lines = plot_errorbar(G_solvs , *ax , n)
        if mode[2] :
            plot_fitting( x , y_avg , y_stdev , ax[0] , lines)
            
            
    fig.set_figwidth(15)
    fig.set_figheight(9)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    return (ax , data , lines)
    

