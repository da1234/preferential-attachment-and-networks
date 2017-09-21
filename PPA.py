import networkx as nx
import random as r
import matplotlib.pyplot as plt
from log_bin_CN_2016 import *
import numpy as np 
from scipy.optimize import curve_fit 
import time 
from scipy import stats
from scipy.stats import chisquare
import math 
import collections
### just for some error stuff 
        
class mInit_error(Exception):
    

    def __init__(self,msg):
        
        
                                
        self.msg = msg
        
    def __str__(self):
        
        return str(self.msg)
        
        
        
        
class No_Kmax_error(Exception):
    

    def __init__(self,msg):
        
        
                                
        self.msg = msg
        
    def __str__(self):
        
        return str(self.msg)




class BA:
    
    
    def __init__(self,max_time,initial_graph_size,no_edges_added,BatchRunning,pref=[True,True,False],rand=[False,True,False],ranW=[[False,True,False],0],runs = 1):
        
        
        
        self.start_time = time.time()
        self.BatchRunning = BatchRunning
        self.no_runs = runs
        
        
        if (initial_graph_size - no_edges_added) != 1:
            
            raise mInit_error('start with a graph size m+1')
        
        self.T = max_time
        self.N0 = initial_graph_size
        self.m = no_edges_added
        self.RandomWalk = ranW[1]
        self.RandomWalk_bool = ranW[0][0]
        self.pref = pref
        self.rand=rand
        self.master_degree_list = []
        self.master_prob_dic_list=[]
        self.master_binned_dic_list=[]
        self.master_count_dic_list=[]
        
        
        self.run()
        self.ini_prob_dicts()
        self.get_k_prob_plot(pref,rand,ranW[0])
    
        
  
            
   
    
    def run(self):
        
        for _ in range(self.no_runs):
        
            self.mGraph = nx.complete_graph(self.N0)
            self.edge_list = self.mGraph.edges()
            self.node_list = self.mGraph.nodes()
            self.critical_exp = 0
            self.drive(self.pref[0],self.rand[0],self.RandomWalk_bool)
            n_degree_list = nx.degree(self.mGraph).values()
            self.master_degree_list.append(n_degree_list)
        
        
   
    def drive(self,pref,rand,ranw):
        
        
        if pref and not rand and not ranw:
        
        
            for i in range(self.T):
                        
                self._addNode_with_pref(self.N0+i) 
                
        if not pref and rand and not ranw:
        
        
            for i in range(self.T):
                        
                self._addNode_with_rand(self.N0+i)
            
        if not pref and not rand and ranw:
            
            
            self.init_neighbours()
            
            for i in range(self.T):
            
                self._addNode_with_randW(self.N0+i)
        
            
        self.mGraph.add_edges_from(self.edge_list)
        
        
    
    def init_neighbours(self):
        
        
        self.neighbour_dic = {}
        
        for i in self.mGraph.nodes():
            
            self.neighbour_dic[i] = self.mGraph.neighbors(i)
        
    
    
    
    def _addNode_with_pref(self,new_node_id):
        
        edge_tuples = self.edge_list
        connected_nodes = []
        cnt = 0 
        
  
        
        while cnt <self.m:
     
            random_edge_pair = r.choice(edge_tuples)
            random_node = r.choice(random_edge_pair)
            
            if random_node not in connected_nodes:
                cnt+=1
                connected_nodes.append(random_node)
                n_edge = (random_node,new_node_id)
                
                self.edge_list.append(n_edge)
                
                
    def _addNode_with_rand(self,new_node_id):
        
        
        connected_nodes = []
        cnt = 0    
        node_list = self.node_list
        
  
        
        while cnt <self.m:
     
            random_node = r.choice(node_list)
            
            if random_node not in connected_nodes:
                cnt+=1
                connected_nodes.append(random_node)
                n_edge = (random_node,new_node_id)
                
                self.edge_list.append(n_edge)
                self.node_list.append(new_node_id)
                
                
    def _addNode_with_randW(self,new_node_id):
        
        
        connected_nodes = []
        cnt = 0    
        node_list = self.node_list
        L= self.RandomWalk
        n_neigbours = []
   
        while cnt <self.m:
     
            random_node = r.choice(node_list)
            
            for walk in range(L):    
                    
                random_node = r.choice(self.neighbour_dic[random_node])
            
            if random_node not in connected_nodes:
                cnt+=1
                
               
            
                connected_nodes.append(random_node)
                n_edge = (random_node,new_node_id)
                n_neigbours.append(random_node)
                self.neighbour_dic[random_node].append(new_node_id)
                self.edge_list.append(n_edge)
                self.node_list.append(new_node_id)
                
            
                self.neighbour_dic[new_node_id]=n_neigbours
        
        #if not len(n_neigbours) == 8 :
        #    print n_neigbours, new_node_id
        #    
  
    def get_theoretical_K_max(self,pref,rand,ranW):
        
      
        
        if pref and not rand and not ranW:
        
        
             return self.get_theoretical_K_max_pref()
                
        elif not pref and rand and not ranW:
            
            return self.get_theoretical_K_max_rand()
        
    
            
        elif not pref and not rand and ranW:
            
             raise No_Kmax_error('Theoretical Kmax for this model could not be determined')
            
        
        self.mGraph.add_edges_from(self.edge_list)
        
        
        
            
    def get_theoretical_K_max_pref(self):
        
        
        
        
        return ((-1./2.) + (0.25+self.m*(self.m+1)*self.T)**0.5)
                
                
    def get_theoretical_K_max_rand(self):
        
        
        return self.m + np.log(self.T)/np.log((self.m+1.)/self.m)
                
                
                
            
    
    def calculate_critical_exp(self):
        
        centres, counts = log_bin(nx.degree(self.mGraph).values(), bin_start=float(self.m), first_bin_width=1., a=1.5, datatype='float')
        
        def func(x,m,c):
            return m*x + c
        
        popt, pcov = curve_fit(func,np.log(np.array(centres)),np.log(np.array(counts)))
        
        print 'this is the gradient',popt[0]
        
        self.critical_exp = popt[0]
        
        self.end_time = time.time()
        
        total = self.end_time - self.start_time
        
        print "this is how long it took", total
        
        
        
    
    def get_k_prob_plot(self,pref,rand,ranW):
        
        if pref[0] and not rand[0] and not ranW[0]:
            
            self.get_k_prob_plot_pref(show_raw=pref[1],collapsed=pref[2])
        
        
        
    
                
        if not pref[0] and rand[0] and not ranW[0]:
            
            self.get_k_prob_plot_rand(show_raw=rand[1],collapsed=rand[2])
        
        
        
           
        if not pref[0] and not rand[0] and ranW[0]:
            
            self.get_k_prob_plot_randW(show_raw=ranW[1],collapsed=ranW[2])
            
    
    def ini_prob_dicts(self):
        
        
        #cnt = 0
        #prob_dictionary = {}
        #for degree in nx.degree(self.mGraph).values():
            
        for degree_list in self.master_degree_list:
            
            cnt = 0
            prob_dictionary = {}
            
            for degree in degree_list:
                
                cnt+=1
                prob_dictionary.setdefault(degree,0)
                prob_dictionary[degree]+=1
            
            count_dic = prob_dictionary.copy()   
            for size,count in prob_dictionary.items():
                
                prob_dictionary[size] /= float(cnt)
                
                
            
            centres, counts = log_bin(degree_list, bin_start=float(self.m), first_bin_width=1., a=1.2, datatype='integer')
            bin_dic = dict(zip(centres, counts))
            
            self.master_prob_dic_list.append(prob_dictionary)
            self.master_count_dic_list.append(count_dic)
            self.master_binned_dic_list.append(bin_dic)
        
        
        
        
    def get_probs(self):
        
        av_prob_dic = {}
        av_bin_dic = {}
        av_count_dic ={}
        
        
        ## averaging the values of the raw probability data ##
        for dic in self.master_prob_dic_list:
            
            for key in dic.keys():
                
                av_prob_dic.setdefault(key,0)
                av_prob_dic[key]+=dic[key]
                
        for size,prob in av_prob_dic.items():
            
            av_prob_dic[size] /= float(len(self.master_prob_dic_list))
            
    
                
                
                
                
          ## averaging the values of the binned probability data ##
        for dic in self.master_binned_dic_list:
            
            for key in dic.keys():
                
                av_bin_dic.setdefault(key,0)
                av_bin_dic[key]+=dic[key]
                
        for size,prob in av_bin_dic.items():
            
            av_bin_dic[size] /= float(len(self.master_binned_dic_list))
 
 
 
        ## averaging the values of the raw count data ##
        for dic in self.master_count_dic_list:
            
            for key in dic.keys():
                
                av_count_dic.setdefault(key,0)
                av_count_dic[key]+=dic[key]
                
        for size,prob in av_prob_dic.items():
            
            av_count_dic[size] /= float(len(self.master_count_dic_list))           
            
            
            
        av_bin_dic =   collections.OrderedDict(sorted(av_bin_dic.items()))
        av_prob_dic =   collections.OrderedDict(sorted(av_prob_dic.items()))
        av_count_dic =   collections.OrderedDict(sorted(av_count_dic.items()))
        

        
        return  av_prob_dic,av_bin_dic,av_count_dic  
        
        
        
        
        
        
    
    def get_k_prob_plot_pref(self,show_raw=False,collapsed=False):
        
        
        av_prob_dic,av_bin_dic,av_count_dic = self.get_probs()
        
        K_max = self.get_theoretical_K_max_pref()
        
        
        if not show_raw and not collapsed:
            
             if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(av_bin_dic.keys()),np.log10(av_bin_dic.values()),'o',label='binned data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),raw=av_count_dic)                
                self.plot_collapsed_probs_pref(av_bin_dic.keys(),av_bin_dic.values())
                plt.show()
                
             else:
    
                plt.grid()
                plt.plot(np.log10(av_bin_dic.keys()),np.log10(av_bin_dic.values()),'o',label='binned data for N is'+str(self.T)+'m is'+str(self.m))
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),collapsed=False,raw=av_count_dic)
                
        
       
                
        elif not show_raw and collapsed:
            
            print 'i got the keys'
            
            K_array=np.array(av_bin_dic.keys(),dtype=float)
            K_array_1 = K_array+1.
            K_array_2 = K_array+2.
            prob_array = np.array(av_bin_dic.values(),dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*K_array*K_array_1*K_array_2
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_pref(av_bin_dic.keys(),av_bin_dic.values())
                plt.show()
                
            else:
    
                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+""+'m is'+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_pref(av_bin_dic.keys(),av_bin_dic.values())
        
       
       
        elif show_raw and collapsed:
            

            K_array=np.array(av_bin_dic.keys(),dtype=float)
            K_array_1 = K_array+1.
            K_array_2 = K_array+2.
            prob_array = np.array(av_bin_dic.values(),dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*K_array*K_array_1*K_array_2
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_pref(av_bin_dic.keys(),av_bin_dic.values())
                plt.show()
                
            else:

                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+"and m is"+" "+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),collapsed=True,raw=av_count_dic)      
                self.plot_collapsed_probs_pref(av_bin_dic.keys(),av_bin_dic.values())
            
            
            
            
        else:
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),raw=av_count_dic)
                plt.show()
                
            else:
                
    
                plt.grid()
                plt.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data for N is'+str(self.T)+' '+'and m is'+' '+str(self.m))
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_pref(av_bin_dic.keys(),av_bin_dic.values(),raw=av_count_dic)
    
                
        
        
        
        
        
    def get_k_prob_plot_rand(self,show_raw=False,collapsed=False):
        
        
        av_prob_dic,av_bin_dic,av_count_dic = self.get_probs()
                
        K_max = self.get_theoretical_K_max_rand()
        
        centres = av_bin_dic.keys()
        counts = av_bin_dic.values()
                
        if not show_raw and not collapsed:
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(centres),np.log10(counts),'o',label='binned data')
                plt.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts,raw=av_count_dic)
                plt.show()
            
            
            else:
 
                plt.grid()
                plt.plot(np.log10(centres),np.log10(counts),'o',label='binned data for N is'+str(self.T)+'m is'+str(self.m))
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts,raw=av_count_dic)
                
               
        elif not show_raw and collapsed:
            K_array=np.array(centres,dtype=float)
            prob_array = np.array(counts,dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*(((self.m+1.)/self.m)**(K_array-self.m))
            
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower right')
                self.plot_exact_k_prob_rand(centres,counts,collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_rand(centres,counts)
                plt.show()
                
            else:
            
                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+""+'m is'+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower right')
                self.plot_exact_k_prob_rand(centres,counts,collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_rand(centres,counts)
                
            
        elif show_raw and collapsed:
            
            K_array=np.array(centres,dtype=float)
            prob_array = np.array(counts,dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*(((self.m+1.)/self.m)**(K_array-self.m))
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts)
                self.plot_collapsed_probs_rand(centres,counts,collapsed=True,raw=av_count_dic)
                plt.show()
                
            else:
    
                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+"and m is"+" "+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts,collapsed=True,raw=av_count_dic)
                self.plot_collapsed_probs_rand(centres,counts)
            
            
    
        else:
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts,raw=av_count_dic)
                plt.show()
            else:
                
                plt.grid()
                plt.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data for N is'+str(self.T)+' '+'and m is'+' '+str(self.m))
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_rand(centres,counts,raw=av_count_dic)
             
                
            
    
    def get_k_prob_plot_randW(self,show_raw=False,collapsed=False):
                
        
        av_prob_dic,av_bin_dic,av_count_dic = self.get_probs()
                
        K_max = self.get_theoretical_K_max_rand()
        
        centres = av_bin_dic.keys()
        counts = av_bin_dic.values()
        
              
            
        if not show_raw and not collapsed:
            
             if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(centres),np.log10(counts),'o',label='binned data')
                plt.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                plt.show()
                
            
             else:
         
                plt.grid()
                plt.plot(np.log10(centres),np.log10(counts),'o',label='binned data for N is'+str(self.T)+'m is'+str(self.m))
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                plt.show()
                
            
    
        elif not show_raw and collapsed:
            K_array=np.array(centres,dtype=float)
            prob_array = np.array(counts,dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*(((self.m+1.)/self.m)**(K_array-self.m))
            
           
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_collapsed_probs_randW(centres,counts)
                plt.show()
                
                
            else:
                #K_plot_Fig = plt.figure()
                #K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+""+'m is'+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')
                self.plot_exact_k_prob_randW(centres,counts,raw=av_count_dic)
                self.plot_collapsed_probs_randW(centres,counts)
        
    
                
        elif show_raw and collapsed:
            
            K_array=np.array(centres,dtype=float)
            prob_array = np.array(counts,dtype=float)
            n_centres = K_array/K_max
            n_counts = prob_array*(((self.m+1.)/self.m)**(K_array-self.m))
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data')
                K_Plot_ax.set_title("Collapsed prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                self.plot_collapsed_probs_randW(centres,counts)
                plt.show()
                
            else:
       
                plt.grid()
                plt.plot(np.log10(n_centres),np.log10(n_counts),'o-',label='Collapsed binned data for N is '+str(self.T)+"and m is"+" "+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')            
                self.plot_collapsed_probs_randW(centres,counts)
                plt.show()
            
        
        else:
            
            if not self.BatchRunning:
                K_plot_Fig = plt.figure()
                K_Plot_ax = K_plot_Fig.add_subplot(111)
                plt.grid()
                K_Plot_ax.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data')
                K_Plot_ax.set_title("Prob of node with degree K")
                K_Plot_ax.set_xlabel("degree, k")
                K_Plot_ax.set_ylabel("Probability")
                K_Plot_ax.legend(loc='lower left')
                plt.show()
                
            else:
                
                plt.grid()
                plt.plot(np.log10(av_prob_dic.keys()),np.log10(av_prob_dic.values()),'o',label='raw data for N is'+str(self.T)+' '+'and m is'+' '+str(self.m))
                plt.plot(np.log10(centres),np.log10(counts),'o-',label=' binned data for N is '+str(self.T)+"and m is"+" "+str(self.m))
                plt.title("Collapsed prob of node with degree K for")
                plt.title("Prob of node with degree K")
                plt.xlabel("degree, k")
                plt.ylabel("Probability")
                plt.legend(loc='lower left')

                
            
        
        
        
        
        
        
        
    def get_exact_solution_pref(self,k):
        
        m = self.m
        
        num = 2.*m*(m+1.)
        denom =k*(k+1)*(k+2)
        prob = num/denom
        return prob 
        
        
        
    def get_exact_solution_rand(self,K):
        
         
        m = self.m
        
        A = (1./(m+1.))
        B = (m/(m+1.))
        prob = A*(B**(K-m))
        return prob 
        
        
        
        
        
    
    def plot_exact_k_prob_pref(self,Ks,observed_probs,collapsed=False,raw=None):
   
        exact_probs = self.get_exact_solution_pref(np.array(Ks))
        K_array=np.array(Ks,dtype=float)
        
        if collapsed:
        
            k_max = self.get_theoretical_K_max_pref()
            K_array_1 = K_array+1.
            K_array_2 = K_array+2.
            n_centres = K_array/k_max
            theoretical_prob = exact_probs*K_array*K_array_1*K_array_2
        
            plt.plot(np.log10(n_centres),np.log10(theoretical_prob),'-',label='exact solution with m='+str(self.m)+'for N is'+str(self.T))
            plt.legend(loc='lower left')
            
            
        else:
            
            plt.plot(np.log10(K_array),np.log10(exact_probs),'-',label='exact solution with m='+str(self.m)+'for N is'+str(self.T))
            plt.legend(loc='lower left')
            
            
        ##for chi squared testing ##
           
        k_range = raw.keys()
        exact = self.get_exact_solution_pref(np.array(k_range))
        self.perform_chi_test(raw.values(),self.T*exact)
        
    
    
    
    def plot_exact_k_prob_rand(self,Ks,observed_probs,collapsed=False,raw=None):
        
        exact_probs = self.get_exact_solution_rand(np.array(Ks))
        
        if collapsed:
   
        
            k_max = self.get_theoretical_K_max_rand()
            theoretical_prob = exact_probs*(((self.m+1.)/self.m)**(np.array(Ks)-self.m))
            plt.plot(np.log10(Ks/k_max),np.log10(theoretical_prob),'-',label='exact solution with m='+str(self.m))
            plt.legend(loc='lower left')
            
        else:
            
            plt.plot(np.log10(Ks),np.log10(exact_probs),'-',label='exact solution with m='+str(self.m))
            plt.legend(loc='lower left')
        
        k_range = raw.keys()
        exact = self.get_exact_solution_rand(np.array(k_range))
        
        self.perform_chi_test(raw.values(),self.T*exact)
        
        
        
        
    def plot_exact_k_prob_randW(self,Ks,observed_probs,collapsed=False,raw=None):
        
        exact_probs = self.get_exact_solution_randW(np.array(Ks))
        
        if collapsed:
            #k_max = self.get_theoretical_K_max_rand()
        
            theoretical_prob = exact_probs*(((self.m+1.)/self.m)**(np.array(Ks)-self.m))
            plt.plot(np.log10(Ks/k_max),np.log10(theoretical_prob),'-',label='exact solution with m='+str(self.m)+"and N is"+""+str(self.T))
            plt.legend(loc='lower left')
            
        else:
            plt.plot(np.log10(Ks),np.log10(exact_probs),'-',label='exact solution with m='+str(self.m)+"and N is"+""+str(self.T))
            plt.legend(loc='lower left')
            
            
            
        
        #k_range = raw.keys()
        #exact = self.get_exact_solution_rand(np.array(k_range))
        
        #self.perform_chi_test(raw.values(),self.T*exact)
    
    
    
    
    
    def plot_collapsed_probs_pref(self,Ks,numerical_data):
        
        P_data = numerical_data
        P_theory = self.get_exact_solution_pref(np.array(Ks))
        
        #print type(P_data)
        
        ratio = P_data/P_theory
        
        if not self.BatchRunning:

            
            Col_plot_Fig = plt.figure()
            Col_Plot_ax = Col_plot_Fig.add_subplot(111)
            plt.grid()
            Col_Plot_ax.plot(np.log10(Ks),ratio,'o-',label='Ratio of Probs')
            Col_Plot_ax.set_title("Ratio of P_data and P_theory against degree K")
            Col_Plot_ax.axhline(y=1.0)
            Col_Plot_ax.set_xlabel("degree, k")
            Col_Plot_ax.set_ylabel("Ratio")
            Col_Plot_ax.legend(loc='lower left')
            plt.show()
            
        else:
            
            
            plt.grid()
            plt.plot(np.log10(Ks),ratio,'o-',label='Ratio of Probs')
            plt.title("Ratio of P_data and P_theory against degree K")
            plt.axhline(y=1.0)
            plt.xlabel("degree, k")
            plt.ylabel("Ratio")
            plt.legend(loc='lower left')
            
            
            
        
      
    def plot_collapsed_probs_rand(self,Ks,numerical_data):
        
        P_data = numerical_data
        P_theory = self.get_exact_solution_rand(np.array(Ks))
        

        ratio = P_data/P_theory
        
        
        if not self.BatchRunning:

                
            Col_plot_Fig = plt.figure()
            Col_Plot_ax = Col_plot_Fig.add_subplot(111)
            plt.grid()
            Col_Plot_ax.plot(Ks,ratio,'o-',label='Ratio of Probs')
            Col_Plot_ax.set_title("Ratio of P_data and P_theory against degree K")
            Col_Plot_ax.axhline(y=1.0)
            Col_Plot_ax.set_xlabel("degree, k")
            Col_Plot_ax.set_ylabel("Ratio")
            Col_Plot_ax.legend(loc='lower right')
            plt.show()
        
        else:
            
            plt.grid()
            plt.plot(np.log10(Ks),ratio,'o-',label='Ratio of Probs')
            plt.title("Ratio of P_data and P_theory against degree K")
            plt.axhline(y=1.0)
            plt.xlabel("degree, k")
            plt.ylabel("Ratio")
            plt.legend(loc='lower right')
            
            

    
    def perform_chi_test(self,observed_data,exact_data):
        
 
     
        a=chisquare(np.array(observed_data,dtype='float64'),f_exp=np.array(exact_data,dtype='float64'))
        
        print 'this is the type of a', type(a)
        print 'chisq and P_val are', a
     
    
        
        
    def draw_Graph(self):
        
        nx.draw_networkx(self.mGraph,pos = nx.circular_layout(self.mGraph),with_labels = False)
        plt.show()
        
        
    
    
    def get_Kmax(self):
        
        return 
        
        
    
        
    
        
        
testG = BA(10000,9,8,BatchRunning=False,pref=[True,False,True],rand=[False,False,True],ranW=[[False,False,False],1],runs =1)
    


    