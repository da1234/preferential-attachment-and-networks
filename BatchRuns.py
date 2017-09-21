

from PPA import BA as model 
import matplotlib.pyplot as plt







class BatchRuns(object):
    
    

    
    
    def __init__(self,Number_runs,m_list,N_list,pref=[True,False,False],rand=[False,True,False],ranW=[[False,True,False],1]):
        
        
        self.K_Maxs = []
        
        self.N_list=N_list
        self.Number_runs=Number_runs
        self.m_list=m_list
        self.pref = pref
        self.rand=rand
        self.ranW=ranW
        
        plt.figure()
        plt.grid()
        #

        for size in N_list:
            for m in m_list:
                G = model(size,m+1,m,BatchRunning=True,pref=self.pref,rand=self.rand,ranW=self.ranW,runs=Number_runs)
                self.K_Maxs.append(G.get_theoretical_K_max(pref=True,rand=False,ranW=False))
        plt.show()
        
     
        
        try:
            self.plot_kMax_T_N()
            
        except:
            
            print 'only do this with one m please, Darrell'
             
            
             
             
             
    def plot_kMax_T_N(self):
        
        plt.figure()
        plt.plot(self.N_list,self.K_Maxs)
        plt.title("Maximum K against System Size")
        plt.ylabel("degree, k")
        plt.xlabel("N")
        plt.legend(loc='lower left')
        plt.show()
        
        
        
        
        
        
                
testB =  BatchRuns(1,[8],[10000,1000])

                
                
                
        
        
        
        
        
        
        
        
        
