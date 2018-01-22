# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:07:19 2018

@author: mondon
"""

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import cPickle
import numpy as np

class results_snfit():
    
    def read_snfit_results(self):
        """
        read snfit results for SNF and return salt parameters in list
        """
        snfit_res_path = '/users/divers/lsst/mondon/science/sugar_analysis_data/results/results_snfit.txt'
        snfit_res = open (snfit_res_path)
        
        self.sn_name = []
        self.z = []
        self.x0 =[]
        self.x0_err = []
        self.x1 =[]
        self.x1_err = []
        self.c = []
        self.c_err = []
        self.mb = []
        self.mb_err = []
        self.cov_x0_x1 = [] 
        self.cov_x0_c = []
        self.cov_x1_c = []
        self.cov_mb_x1 = []
        self.cov_mb_c = []
        
        snfit_lines = snfit_res.readlines()
        
        line_name = 0
        
        for line in range(len(snfit_lines)):
            
            snfit_lines[line] = snfit_lines[line].strip()
            snfit_lines[line] = snfit_lines[line].split(' ')
            snfit_lines[line]=[x.replace('\n','') for x in snfit_lines[line]]
            snfit_lines[line]=[x.replace("'",'') for x in snfit_lines[line]]
            
            if snfit_lines[line][0] == 'sn_name':            
                self.sn_name.append(snfit_lines[line][1])
                line_name = line
                
            if snfit_lines[line][0] == 'Redshift' and line == line_name + 35 :            
                self.z.append(float(snfit_lines[line][1]))    
                
            if snfit_lines[line][0] == 'X0' and line == line_name + 37:            
                self.x0.append(float(snfit_lines[line][1]))  
                self.x0_err.append(float(snfit_lines[line][2]))
    
            if snfit_lines[line][0] == 'X1' and line == line_name + 38 :            
                self.x1.append(float(snfit_lines[line][1]))  
                self.x1_err.append(float(snfit_lines[line][2]))
                
            if snfit_lines[line][0] == 'Color' and line == line_name + 36:            
                self.c.append(float(snfit_lines[line][1]))
                self.c_err.append(float(snfit_lines[line][2]))
    
            if snfit_lines[line][0] == 'RestFrameMag_0_B' and line == line_name + 39 :            
                self.mb.append(float(snfit_lines[line][1]))  
                self.mb_err.append(float(snfit_lines[line][2]))
    
            if snfit_lines[line][0] == 'CovX0X1' and line == line_name + 52:            
                self.cov_x0_x1.append(float(snfit_lines[line][1]))
    
            if snfit_lines[line][0] == 'CovColorX0' and line == line_name + 43:            
                self.cov_x0_c.append(float(snfit_lines[line][1]))           
                
     
            if snfit_lines[line][0] == 'CovColorX1' and line == line_name + 44:            
                self.cov_x1_c.append(float(snfit_lines[line][1]))
    
            if snfit_lines[line][0] == 'CovColorRestFrameMag_0_B' and line == line_name + 42:            
                self.cov_mb_c.append(float(snfit_lines[line][1]))
     
            if snfit_lines[line][0] == 'CovRestFrameMag_0_BX1' and line == line_name + 50:            
                self.cov_mb_x1.append(float(snfit_lines[line][1]))    
            
        self.sn_name = np.array(self.sn_name,str)
        self.z = np.array(self.z)
        self.x0 = np.array(self.x0)
        self.x0_err = np.array(self.x0_err)
        self.x1 = np.array(self.x1)
        self.x1_err = np.array(self.x1_err) 
        self.c = np.array(self.c)
        self.c_err = np.array(self.c_err)
        self.mb = np.array(self.mb)
        self.mb_err = np.array(self.mb_err)
        self.cov_x0_x1 = np.array(self.cov_x0_x1)
        self.cov_x0_c = np.array(self.cov_x0_c)
        self.cov_x1_c = np.array(self.cov_x1_c)
        self.cov_mb_x1 = np.array(self.cov_mb_x1)
        self.cov_mb_c = np.array(self.cov_mb_c)        
        
    def read_meta(self):
        """
        read meta.pkl and return salt2 parameters in array
        """
        meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        self.meta_sn_name_list = []
        self.meta_z = []
        self.meta_x0 =[]
        self.meta_x0_err = []
        self.meta_x1 =[]
        self.meta_x1_err = []
        self.meta_c = []
        self.meta_c_err = []
        self.meta_mb = []
        self.meta_mb_err = []
        self.meta_cov_x0_x1 = [] 
        self.meta_cov_x0_c = []
        self.meta_cov_x1_c = []
        self.meta_cov_mb_x1 = []
        self.meta_cov_mb_c = []  
        
        for meta_sn_name in meta.keys(): 
            
            if meta[meta_sn_name]['idr.subset'] != 'bad':
            
                self.meta_sn_name_list.append(meta_sn_name)
                self.meta_z.append(meta[meta_sn_name]['host.zcmb'])
                self.meta_x0.append(meta[meta_sn_name]['salt2.X0'])
                self.meta_x0_err.append(meta[meta_sn_name]['salt2.X0.err'])
                self.meta_x1.append(meta[meta_sn_name]['salt2.X1'])
                self.meta_x1_err.append(meta[meta_sn_name]['salt2.X1.err'])
                self.meta_c.append(meta[meta_sn_name]['salt2.Color'])
                self.meta_c_err.append(meta[meta_sn_name]['salt2.Color.err'])
                self.meta_mb.append(meta[meta_sn_name]['salt2.RestFrameMag_0_B'])
                self.meta_mb_err.append(meta[meta_sn_name]['salt2.RestFrameMag_0_B.err'])
                self.meta_cov_x0_x1.append(meta[meta_sn_name]['salt2.CovX0X1'])
                self.meta_cov_x0_c.append(meta[meta_sn_name]['salt2.CovColorX0'])
                self.meta_cov_x1_c.append(meta[meta_sn_name]['salt2.CovColorX1'])
                self.meta_cov_mb_x1.append(meta[meta_sn_name]['salt2.CovRestFrameMag_0_BX1'])
                self.meta_cov_mb_c.append(meta[meta_sn_name]['salt2.CovColorRestFrameMag_0_B'])
        
        
        self.meta_z = np.array(self.meta_z)
        self.meta_x0 = np.array(self.meta_x0)
        self.meta_x0_err = np.array(self.meta_x0_err)
        self.meta_x1 = np.array(self.meta_x1)
        self.meta_x1_err = np.array(self.meta_x1_err) 
        self.meta_c = np.array(self.meta_c)
        self.meta_c_err = np.array(self.meta_c_err)
        self.meta_mb = np.array(self.meta_mb)
        self.meta_mb_err = np.array(self.meta_mb_err)
        self.meta_cov_x0_x1 = np.array(self.meta_cov_x0_x1)
        self.meta_cov_x0_c = np.array(self.meta_cov_x0_c)
        self.meta_cov_x1_c = np.array(self.meta_cov_x1_c)
        self.meta_cov_mb_x1 = np.array(self.meta_cov_mb_x1)
        self.meta_cov_mb_c = np.array(self.meta_cov_mb_c)   
        
    def read_sncosmo(self):
        self.sncosmo_res = np.loadtxt('../sugar_analysis_data/results/res_salt2_SNF.txt',dtype='str')
        self.sncosmo_sn_name = np.array(self.sncosmo_res[:,0],str)
        self.sncosmo_z = np.array(self.sncosmo_res[:,1],float)
        self.sncosmo_x0 = np.array(self.sncosmo_res[:,13],float)
        self.sncosmo_x0_err = np.array(self.sncosmo_res[:,14],float)
        self.sncosmo_x1 = np.array(self.sncosmo_res[:,6],float)
        self.sncosmo_x1_err = np.array(self.sncosmo_res[:,7],float)
        self.sncosmo_c = np.array(self.sncosmo_res[:,8],float)
        self.sncosmo_c_err = np.array(self.sncosmo_res[:,9],float)
        self.sncosmo_mb = np.array(self.sncosmo_res[:,4],float)
        self.sncosmo_mb_err = np.array(self.sncosmo_res[:,5],float)
        self.sncosmo_cov_x1_c = np.array(self.sncosmo_res[:,12],float)
        self.sncosmo_cov_mb_x1 = np.array(self.sncosmo_res[:,10],float)
        self.sncosmo_cov_mb_c = np.array(self.sncosmo_res[:,11],float)       
        
    def plot_hist_snfit_meta(self):
        """
        plot histogramm of the difference between salt2 parameters using snfit and salt2 parameters in meta.pkl 
        """     
        
        self.read_meta()
        self.read_snfit_results()

            
        self.diff_x0 = self.meta_x0 - self.x0
        self.diff_x0_err = self.meta_x0_err - self.x0_err

        
        
        
        
        self.diff_x1 = self.meta_x1 - self.x1
        self.diff_x1_err = self.meta_x1_err - self.x1_err        
        self.diff_c = self.meta_c - self.c
        self.diff_c_err = self.meta_c_err - self.c_err         
        self.diff_mb = self.meta_mb - self.mb
        self.diff_mb_err = self.meta_mb_err - self.mb_err 
        self.diff_cov_x0_x1 = self.meta_cov_x0_x1 - self.cov_x0_x1
        self.diff_cov_x0_c = self.meta_cov_x0_c - self.cov_x0_c
        self.diff_cov_x1_c = self.meta_cov_x1_c -self.cov_x1_c
        self.diff_cov_mb_x1 = self.meta_cov_mb_x1 - self.cov_mb_x1
        self.diff_cov_mb_c = self.meta_cov_mb_c - self.cov_mb_c
 
        for i in range(len(self.diff_x1)):
            if self.diff_x1[i] > 0.0001:
                print self.diff_x1[i],  self.sn_name[i], self.meta_sn_name_list[i], self.x1[i], self.meta_x1[i]

        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x0,25,label='X0')
        ax0_2.hist(self.diff_x0_err,25,label='X0 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x1,25,label='X1')
        ax0_2.hist(self.diff_x1_err,25,label='X1 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_c,25,label='Color')
        ax0_2.hist(self.diff_c_err,25,label='Color error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()
        
    def plot_hist_snfit_sncosmo(self):
        """
        plot histogramm of the difference between salt2 parameters using snfit and salt2 parameters in meta.pkl 
        """
        
        self.read_sncosmo()
        self.read_snfit_results()
        self.diff_x0_sncosmo = []
        self.diff_x0_err_sncosmo = []
        self.diff_x1_sncosmo = []
        self.diff_x1_err_sncosmo = []        
        self.diff_c_sncosmo = []
        self.diff_c_err_sncosmo = []      
        self.diff_mb_sncosmo = []
        self.diff_mb_err_sncosmo = [] 
        self.diff_cov_x0_x1_sncosmo = []
        self.diff_cov_x0_c_sncosmo = []
        self.diff_cov_x1_c_sncosmo = []
        self.diff_cov_mb_x1_sncosmo = []
        self.diff_cov_mb_c_sncosmo = []

        for i in range (len(self.sn_name)):
            for j in range (len(self.sncosmo_sn_name)):
                if self.sn_name[i] == self.sncosmo_sn_name[j]:
                    if np.abs(self.x1[i] - self.sncosmo_x1[j]) < 0.05:
                        self.diff_x0_sncosmo.append(self.x0[i] - self.sncosmo_x0[j])
                        self.diff_x0_err_sncosmo.append(self.x0_err[i] - self.sncosmo_x0_err[j])
                        self.diff_x1_sncosmo.append(self.x1[i] - self.sncosmo_x1[j])
                        self.diff_x1_err_sncosmo.append(self.x1_err[i] - self.sncosmo_x1_err[j])      
                        self.diff_c_sncosmo.append(self.c[i] - self.sncosmo_c[j])
                        self.diff_c_err_sncosmo.append(self.c_err[i] - self.sncosmo_c_err[j])     
                        self.diff_mb_sncosmo.append(self.mb[i] - self.sncosmo_mb[j])
                        self.diff_mb_err_sncosmo.append(self.mb_err[i] - self.sncosmo_mb_err[j])
#                    self.diff_cov_x0_x1_sncosmo.append()
#                    self.diff_cov_x0_c_sncosmo.append()
#                    self.diff_cov_x1_c_sncosmo.append()
#                    self.diff_cov_mb_x1_sncosmo.append()
#                    self.diff_cov_mb_c_sncosmo.append()
                    else:
                        print self.x1[i] - self.sncosmo_x1[j],  self.sn_name[i],self.sncosmo_sn_name[j],  self.x1[i], self.sncosmo_x1[j]

        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x0_sncosmo,50,label='x0')
        ax0_2.hist(self.diff_x0_err_sncosmo,50,label='x0 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x1_sncosmo,50,label='X1')
        ax0_2.hist(self.diff_x1_err_sncosmo,50,label='X1 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_c_sncosmo,50,label='Color')
        ax0_2.hist(self.diff_c_err_sncosmo,50,label='Color error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.show()