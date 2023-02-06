import plot
from plot import plt as plt
import numpy as np

exp1 = -265.9 # proton solvation energy (https://pubs.acs.org/doi/10.1021/jp063552y) Kelly
exp2 = -262.23 # proton solvation energy (https://aip.scitation.org/doi/10.1063/1.477096) G. J. Tawa

def extract_molecule():

    data = []
    k = 0
    for k , (i , j) in enumerate([("MEA" , 'default') , ("MEA" , 'w1.0') , ("p-MEA" , 'default') , ("p-MEA" , 'w1.0') , ("Proton" , "default")]) :
    
        print(plot.bcolors.HEADER + f"cur : {i} with {j}" + plot.bcolors.ENDC)
        ax , data = plot.plot_cluster(i,j , (2,1) , (True ,True ,True) , n = [0 , 500] , path = "plot_py/G_solv.txt" , preData=data)
        plt.figure(k+1)
        plt.suptitle(f"{i} with wall scale : {j}" , fontsize=25)
        # if i == "MEA":
        #     ax[0].set_ylim([-17.5 , -6])
        # elif i == "p-MEA":
        #     ax[0].set_ylim([-170 , -100])
        
        plt.savefig(f"figtmp/{i}_{j}.png")
        
def replace_none(lines , set_label):
    for i,l in enumerate(set_label):
        if l == None:
            set_label[i] = lines[i].get_label()
            
    return set_label
    

def main():
    
    data = []
    k = 0
    for k , (i , j) in enumerate([("MEA" , 'default') , ("MEA" , 'w1.0') , ("p-MEA" , 'default') , ("p-MEA" , 'w1.0') , ("Proton" , "default")]) :
    
        print(plot.bcolors.HEADER + f"cur : {i} with {j}" + plot.bcolors.ENDC)
        ax , data , lines = plot.plot_cluster(i,j , (2,1) , (True ,True ,True) , n = [0 , 500] , path = "plot_py/G_solv.txt" , preData=data)
        
        if i == 'Proton':
            ax[0].axhline(y = exp1, color = 'r', linestyle = '-' , alpha = 0.6)
            ax[0].axhline(y = exp2, color = 'r', linestyle = 'solid', alpha = 0.6)

        ax[0].legend(lines , replace_none(lines , [ None , None , 'G_solv' , None]) , loc='best',frameon=False)
        plt.figure(k+1)
        plt.suptitle(f"{i} with wall scale : {j}" , fontsize=25)
        # if i == "MEA":
        #     ax[0].set_ylim([-17.5 , -6])
        # elif i == "p-MEA":
        #     ax[0].set_ylim([-170 , -100])
        

        plt.savefig(f"figtmp/{i}_{j}.png")
        

    # for i in range(2):
        
    #     diff = dict()
    #     fig , ax = plt.subplots(2,1)

    #     for num,values in data[2+i].items():
    #         tmp = []
    #         for j in values:
    #             tmp.extend([j-k for k in data[4][num]])
    #         diff[num] = tmp
    #     x , y_avg , y_stdev = plot.plot_errorbar(diff , *ax )
    #     plot.plot_fitting( x , y_avg , y_stdev , ax[0])
        
    #     fig.set_figwidth(10)
    #     fig.set_figheight(6)
    #     fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    #     plt.figure(k+i+2)

    #     if i == 0:
    #         plt.suptitle("G$_{solv}$ (p-MEA_default - proton)" , fontsize=25)
            
    #     if i == 1:    
    #         plt.suptitle("G$_{solv}$ (p-MEA_w1.0 - proton)" , fontsize=25)

    #     if i == 0:
    #         plt.savefig("figtmp/p-MEA_default_M_proton.png")

    #     if i == 1:    
    #         plt.savefig("figtmp/p-MEA_w1.0_M_proton.png")
            
            
    for i in range(2):
        
        diff = dict()
        fig , ax = plt.subplots(2,1)

        for num,values in data[2+i].items():
            tmp = []
            for j in values:
                for m in data[i][num]:
                    tmp.extend([-np.log10(np.exp(-((41.85 + m)+(k + (-6.27) )- (49.04 + j))/1.987*1000/298)) for k in data[4][num]])
                    #tmp.extend([-np.log10(np.exp(-((m)+(k  )- ( j) + 893.0/4.18)/1.987*1000/298)) for k in data[4][num]])
            diff[num] = tmp
        x , y_avg , y_stdev , lines = plot.plot_errorbar(diff , *ax )
        plot.plot_fitting( x , y_avg , y_stdev , ax[0] , lines)

        ax[0].set_label('pKa')
        
        ax[0].set_ylabel('pK$_{a}$',fontsize=20)
            
        ax[0].legend(lines , replace_none(lines , [ None , None , 'pKa' , None]), loc='best',frameon=False)

        fig.set_figwidth(15)
        fig.set_figheight(9)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        plt.figure(k+i+2)

        if i == 0:
            #plt.suptitle("G$_{reaction,aq}$ (default)" , fontsize=25)
            plt.suptitle("pK$_{a}$ (default)" , fontsize=25)
            
        if i == 1:    
            #plt.suptitle("G$_{reaction,aq}$ (w1.0)" , fontsize=25)
            plt.suptitle("pK$_{a}$ (w1.0)" , fontsize=25)

        if i == 0:
            #plt.savefig("figtmp/G_reaction_default.png")
            plt.savefig("figtmp/pKa_default.png")

        if i == 1:    
            #plt.savefig("figtmp/G_reaction_w1.0.png")
            plt.savefig("figtmp/pKa_w1.0.png")





    plt.show()
    # read_clusters("plot_py/final_ensemble.xyz")
    
    
if __name__ == '__main__':
    main()
    

