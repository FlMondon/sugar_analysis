# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:07:19 2018

@author: mondon
"""

from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
import cPickle
import numpy as np
SUGAR_parameter_pkl = '../sugar/sugar/data_output/data_output/sugar_parameters.pkl'

class results_snfit():
    def __init__(self, width):
        self.width = width
        
    def read_snfit_results(self, snfit_res_path='../sugar_analysis_data/results/results_snfit_1A.txt'):
        """
        read snfit results for SNF and return salt parameters in list
        """

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
        self.snfit_chi2 = []
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

            if snfit_lines[line][0] == '@CHI2_LC' and line == line_name + 56:            
                self.snfit_chi2.append(float(snfit_lines[line][1]))  
            
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
        self.snfit_chi2 = np.array(self.snfit_chi2)
        
    def read_meta(self):
        """
        read meta.pkl and return salt2 parameters in array
        """
        meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        self.meta_sn_name_list = []
        self.meta_zcmb = []
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
        self.meta_zhl = []
        self.meta_zhl_err = []
        for meta_sn_name in meta.keys(): 
            
            if meta[meta_sn_name]['idr.subset'] != 'bad' and meta[meta_sn_name]['idr.subset'] != 'auxiliary':
            
                self.meta_sn_name_list.append(meta_sn_name)
                self.meta_zhl_err.append(meta[meta_sn_name]['host.zhelio.err'])
                self.meta_zhl.append(meta[meta_sn_name]['host.zhelio'])
                self.meta_zcmb.append(meta[meta_sn_name]['host.zcmb'])
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
        
        
        self.meta_zcmb = np.array(self.meta_zcmb)
        self.meta_zhl = np.array(self.meta_zhl)
        self.meta_zhl_err = np.array(self.meta_zhl_err)
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
        
    def read_sncosmo(self, path='../sugar_analysis_data/results/res_salt2_SNF_1A.txt'):
        self.sncosmo_res = np.loadtxt(path ,dtype='str')
        self.sncosmo_sn_name = np.array(self.sncosmo_res[:,0],str)
        self.sncosmo_z = np.array(self.sncosmo_res[:,2],float)
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
        self.sncosmo_chi2 = np.array(self.sncosmo_res[:,17],float)    
        
    def plot_hist_snfit_meta(self):
        """
        plot histogramm of the difference between salt2 parameters using snfit and salt2 parameters in meta.pkl 
        """     
        
        self.read_meta()
        self.read_snfit_results()

            
        self.diff_x0 = []
        self.diff_x0_err = []
        self.diff_x1 = []
        self.diff_x1_err = []        
        self.diff_c = []
        self.diff_c_err = []      
        self.diff_mb = []
        self.diff_mb_err = [] 
        self.diff_cov_x0_x1 = []
        self.diff_cov_x0_c = []
        self.diff_cov_x1_c = []
        self.diff_cov_mb_x1 = []
        self.diff_cov_mb_c = []

        for i in range (len(self.sn_name)):
            for j in range (len(self.meta_sn_name_list)):
                if self.sn_name[i] == self.meta_sn_name_list[j]:
                    if np.abs(self.x1[i] - self.meta_x1[j]) < 0.001:
                        self.diff_x0.append(self.x0[i] - self.meta_x0[j])
                        self.diff_x0_err.append(self.x0_err[i] - self.meta_x0_err[j])
                        self.diff_x1.append(self.x1[i] - self.meta_x1[j])
                        self.diff_x1_err.append(self.x1_err[i] - self.meta_x1_err[j])      
                        self.diff_c.append(self.c[i] - self.meta_c[j])
                        self.diff_c_err.append(self.c_err[i] - self.meta_c_err[j])     
                        self.diff_mb.append(self.mb[i] - self.meta_mb[j])
                        self.diff_mb_err.append(self.mb_err[i] - self.meta_mb_err[j])
#                    self.diff_cov_x0_x1.append()
#                    self.diff_cov_x0_c.append()
#                    self.diff_cov_x1_c.append()
#                    self.diff_cov_mb_x1.append()
#                    self.diff_cov_mb_c.append()
                    else:
                        print self.x1[i] - self.meta_x1[j],  self.sn_name[i],self.meta_sn_name_list[j],  self.x1[i], self.meta_x1[j]

        rcParams['font.size'] = 16.
        font = {'family': 'normal', 'size': 16}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=25)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)
        fig = plt.figure(figsize=(8.,8.))       
                
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
        pdffile = '../sugar_analysis_data/results/x0_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
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
        pdffile = '../sugar_analysis_data/results/x1_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
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
        pdffile = '../sugar_analysis_data/results/color_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()

        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_mb,50,label='mb')
        ax0_2.hist(self.diff_mb_err,50,label='mb error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        pdffile = '../sugar_analysis_data/results/mb_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
    def plot_hist_snfit_sncosmo(self):
        """
        plot histogramm of the difference between salt2 parameters using snfit and salt2 parameters in meta.pkl 
        """
        
        self.read_sncosmo(path='../sugar_analysis_data/results/res_salt2_SNF_'+str(self.width)+'.txt')
        self.read_snfit_results(snfit_res_path='../sugar_analysis_data/results/results_snfit_'+str(self.width)+'.txt')

#        self.read_sncosmo(path='../sugar_analysis_data/results/res_salt2_SNF_GF.txt')
#        self.read_snfit_results(snfit_res_path='../sugar_analysis_data/results/results_snfit_GF.txt')

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
        self.diff_chi2 = []
        for i in range (len(self.sn_name)):
            for j in range (len(self.sncosmo_sn_name)):
                if self.sn_name[i] == self.sncosmo_sn_name[j]:
                    if np.abs(self.x1[i] - self.sncosmo_x1[j]) < 0.02:
                        self.diff_x0_sncosmo.append(self.x0[i] - self.sncosmo_x0[j])
                        self.diff_x0_err_sncosmo.append(self.x0_err[i] - self.sncosmo_x0_err[j])
                        self.diff_x1_sncosmo.append(self.x1[i] - self.sncosmo_x1[j])
                        self.diff_x1_err_sncosmo.append(self.x1_err[i] - self.sncosmo_x1_err[j])      
                        self.diff_c_sncosmo.append(self.c[i] - self.sncosmo_c[j])
                        self.diff_c_err_sncosmo.append(self.c_err[i] - self.sncosmo_c_err[j])     
                        self.diff_mb_sncosmo.append(self.mb[i] - self.sncosmo_mb[j])
                        self.diff_mb_err_sncosmo.append(self.mb_err[i] - self.sncosmo_mb_err[j])
                        self.diff_chi2.append(self.snfit_chi2[i] - self.sncosmo_chi2[j])
#                    self.diff_cov_x0_x1_sncosmo.append()
#                    self.diff_cov_x0_c_sncosmo.append()
#                    self.diff_cov_x1_c_sncosmo.append()
#                    self.diff_cov_mb_x1_sncosmo.append()
#                    self.diff_cov_mb_c_sncosmo.append()
                    else:
                        print self.x1[i] - self.sncosmo_x1[j],  self.sn_name[i],self.sncosmo_sn_name[j],  self.x1[i], self.sncosmo_x1[j]

#        rcParams['font.size'] = 16.
#        font = {'family': 'normal', 'size': 16}
#        rc('axes', linewidth=1.5)
#        rc("text", usetex=True)
#        rc('font', family='serif')
#        rc('font', serif='Times')
#        rc('legend', fontsize=25)
#        rc('xtick.major', size=5, width=1.5)
#        rc('ytick.major', size=5, width=1.5)
#        rc('xtick.minor', size=3, width=1)
#        rc('ytick.minor', size=3, width=1)
#        fig = plt.figure(figsize=(8.,8.))       
#        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x0_sncosmo,50,label='x0_'+str(self.width))
        ax0_2.hist(self.diff_x0_err_sncosmo,50,label='x0 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
#        ax0_1.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        pdffile = '../sugar_analysis_data/results/x0_plot_'+str(self.width)+'.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x1_sncosmo,50,label='X1_'+str(self.width))
        ax0_2.hist(self.diff_x1_err_sncosmo,50,label='X1 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        pdffile = '../sugar_analysis_data/results/x1_plot_'+str(self.width)+'.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_c_sncosmo,50,label='Color_'+str(self.width))
        ax0_2.hist(self.diff_c_err_sncosmo,50,label='Color error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        pdffile = '../sugar_analysis_data/results/color_plot_'+str(self.width)+'.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()

        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_mb_sncosmo,50,label='mb_'+str(self.width))
        ax0_2.hist(self.diff_mb_err_sncosmo,50,label='mb error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        pdffile = '../sugar_analysis_data/results/mb_plot_'+str(self.width)+'.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()

        plt.hist(self.diff_chi2,50,label='chi2_'+str(self.width))
        pdffile = '../sugar_analysis_data/results/chi2_'+str(self.width)+'.pdf'
        plt.legend()
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()   
        
    def diff_sncosmo_dL(self):
        diff_mb_sncosmo = [] 
        outlayers = ['SN2012cu','SN2009hi','SNF20061022-014','SNNGC7589']
        self.read_sncosmo()
        sncosmo_sn_name_5A = self.sncosmo_sn_name
        sncosmo_mb_5A = self.sncosmo_mb
        self.read_sncosmo(path='../sugar_analysis_data/results/res_salt2_SNF_'+str(self.width)+'.txt')
#        self.read_sncosmo(path='../sugar_analysis_data/results/res_salt2_SNF_save.txt')

        sncosmo_mb_1A = self.sncosmo_mb
        sncosmo_sn_name_1A = self.sncosmo_sn_name 
        for i in range (len(sncosmo_sn_name_5A)):
            for j in range (len(sncosmo_sn_name_1A)):
                
                if sncosmo_sn_name_5A[i] == sncosmo_sn_name_1A[j] and sncosmo_sn_name_5A[i] not in outlayers:  
#                    if abs(sncosmo_mb_5A[i] - sncosmo_mb_1A[j]) < 0.01:
                        diff_mb_sncosmo.append(sncosmo_mb_5A[i] - sncosmo_mb_1A[j])
#                    else:
#                        print sncosmo_sn_name_5A[i]
        mean_hist = np.mean(diff_mb_sncosmo)
        sig_hist = np.std(diff_mb_sncosmo)
        print mean_hist, sig_hist

#        rcParams['font.size'] = 16.
#        font = {'family': 'normal', 'size': 16}
#        rc('axes', linewidth=1.5)
#        rc("text", usetex=True)
#        rc('font', family='serif')
#        rc('font', serif='Times')
#        rc('legend', fontsize=25)
#        rc('xtick.major', size=5, width=1.5)
#        rc('ytick.major', size=5, width=1.5)
#        rc('xtick.minor', size=3, width=1)
#        rc('ytick.minor', size=3, width=1)
#        fig = plt.figure(figsize=(8.,8.))
#        plt.text(0.0004, 30, '\mu = %f \n'%(mean_hist)+' \sigma = %f'%(sig_hist), fontweight = 'bold', fontsize = 20)
        plt.hist(diff_mb_sncosmo,40,label='diff mb th1A '+str(self.width)+' sncosmo')
#        pdffile = '../sugar_analysis_data/results/diff_mb_th1A_'+str(self.width)+'_sncosmo.pdf'

##        plt.hist(diff_mb_sncosmo,40,label='diff mb th 5A 1A sncosmo')
##        pdffile = '../sugar_analysis_data/results/diff_mb_th_5A_1A_sncosmo.pdf' 

#        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
#        plt.legend()
#        plt.savefig(pdffile, bbox_inches='tight')
        plt.show() 
        
    def diff_snfit_dL(self):
        diff_mb = [] 
        outlayers = ['SN2012cu','SN2009hi','SNF20061022-014','SNNGC7589']
#        self.read_snfit_results()
#        snfit_sn_name_5A = self.sn_name
#        snfit_mb_5A = self.mb
        self.read_sncosmo()
        snfit_sn_name_5A = self.sncosmo_sn_name
        snfit_mb_5A = self.sncosmo_mb
        
        self.read_snfit_results(snfit_res_path='../sugar_analysis_data/results/results_snfit_'+str(self.width)+'.txt')
#        self.read_snfit_results(snfit_res_path='../sugar_analysis_data/results/results_snfit_1A.txt')
        snfit_mb_1A = self.mb
        snfit_sn_name_1A = self.sn_name 
        for i in range (len(snfit_sn_name_5A)):
            for j in range (len(snfit_sn_name_1A)):
                if snfit_sn_name_5A[i] == snfit_sn_name_1A[j] and snfit_sn_name_5A[i] not in outlayers: 
#                    if abs(snfit_mb_5A[i] - snfit_mb_1A[j]) < 0.01:
                        diff_mb.append(snfit_mb_5A[i] - snfit_mb_1A[j])
#                    else:
#                        print snfit_sn_name_5A[i]
        mean_hist = np.mean(diff_mb)
        sig_hist = np.std(diff_mb)
        print mean_hist, sig_hist

#        rcParams['font.size'] = 16.
#        font = {'family': 'normal', 'size': 16}
#        rc('axes', linewidth=1.5)
#        rc("text", usetex=True)
#        rc('font', family='serif')
#        rc('font', serif='Times')
#        rc('legend', fontsize=25)
#        rc('xtick.major', size=5, width=1.5)
#        rc('ytick.major', size=5, width=1.5)
#        rc('xtick.minor', size=3, width=1)
#        rc('ytick.minor', size=3, width=1)
#        fig = plt.figure(figsize=(8.,8.))
#        plt.text(-0.009, 12, '\mu = %f \n'%(mean_hist)+' \sigma = %f'%(sig_hist), fontweight = 'bold', fontsize = 20)
        plt.hist(diff_mb,40,label='diff mb th1A '+str(self.width)+' snfit')
#        pdffile = '../sugar_analysis_data/results/diff_mb_th1A_'+str(self.width)+'_snfit.pdf'
        
##        plt.hist(diff_mb,40,label='diff mb th 5A 1A snfit')
##        pdffile = '../sugar_analysis_data/results/diff_mb_th_5A_1A_snfit.pdf'   
        
#        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
#        plt.legend()
#        plt.savefig(pdffile, bbox_inches='tight')
#        plt.show() 

    def plot_mean(self):
        rcParams['font.size'] = 16.
        font = {'family': 'normal', 'size': 16}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=25)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)
        fig = plt.figure(figsize=(8.,8.))
        x = [0, 1, 1.5, 5,10, 20, 40]
        y_fit = [0.0048,0.0031,  0.0025, 0.000087,-0.00022,-0.00008,0.001]
        y_fit_errors = [0.0016, 0.00090,  0.00083, 0.00059 ,0.00063,0.00095,0.0024]
        y_cosmo = [-0.00098,  0.0015, -0.00083, -0.00014,-0.00021, 0.00006, 0.00099]
        y_cosmo_errors = [0.0014, 0.00032,  0.000044,0.00011 ,0.00018,0.00065, 0.0022]
        x_1A_fit = [0]
        y_1A_fit = [0.0048]
        y_1A_fit_err = [0.0016]
        plt.errorbar(x,y_fit,yerr=y_fit_errors,fmt='o',label='Snfit 5A')
        plt.errorbar(x,y_cosmo,yerr=y_cosmo_errors,fmt='o',label='Sncosmo 5A')
        plt.errorbar(x_1A_fit, y_1A_fit, yerr=y_1A_fit_err, fmt='^', label='Snfit 1A')
        plt.plot([-5, 50], [0,0])
        plt.xlim(-5,50)
        plt.xlabel('width')
        plt.ylabel('mean')
        plt.legend()
        pdffile = '../sugar_analysis_data/results/mean.pdf' 
        plt.savefig(pdffile, bbox_inches='tight')

        plt.show()
        
    def HD_input_snfit_data(self):
        """
        Sample snfit results in HD input
        """

        dico = cPickle.load(open(SUGAR_parameter_pkl))
        self.read_snfit_results()
        self.read_meta()
        Filtre = np.array([True]*len(self.sn_name))
        self.zcmb = []
        self.z_err = []
        for j in range(len(self.sn_name)):
            if self.sn_name[j] in dico.keys() and self.sn_name[j] :

                for i in range (len(self.meta_sn_name_list)):
                    if self.sn_name[j] == self.meta_sn_name_list[i]:
                        
                        self.z_err.append(self.meta_zhl_err[i])
                        self.zcmb.append(self.meta_zcmb[i])
                        if np.abs(self.x1[j] - self.meta_x1[i]) > 0.001:
                            print 'problem with %s include in sample but difference between snfit and meta'%(self.sn_name[j])
            else:
                Filtre[j] = False

        for p in dico.keys():
            if p not in self.sn_name:
                print p
            
        self.x1 = self.x1[Filtre]
        self.x1_err = self.x1_err[Filtre] 
        self.c = self.c[Filtre]
        self.c_err = self.c_err[Filtre]
        self.mb = self.mb[Filtre]
        self.mb_err = self.mb_err[Filtre]
        self.cov_x0_x1 = self.cov_x0_x1[Filtre]
        self.cov_x0_c = self.cov_x0_c[Filtre]
        self.cov_x1_c = self.cov_x1_c[Filtre]
        self.cov_mb_x1 = self.cov_mb_x1[Filtre]
        self.cov_mb_c = self.cov_mb_c[Filtre]
        self.z = self.z[Filtre]
        self.zcmb = np.array(self.zcmb)
        self.z_err = np.array(self.z_err)

        self.cov_y = np.zeros((len(self.mb)*3,len(self.mb)*3))

        for i in range (len(self.mb)):
            self.cov_y[i*3,i*3] = self.mb_err[i]**2
            self.cov_y[i*3+ 1,i*3+ 1] = self.x1_err[i]**2
            
            self.cov_y[i*3+ 2,i*3+ 2] = self.c_err[i]**2
            self.cov_y[i*3+ 0,i*3+ 1] = self.cov_mb_x1[i]
            self.cov_y[i*3+ 1,i*3+ 0] = self.cov_mb_x1[i]
            self.cov_y[i*3+ 0,i*3+ 2] = self.cov_mb_c[i]
            self.cov_y[i*3+ 2,i*3+ 0] = self.cov_mb_c[i]
            self.cov_y[i*3+ 1,i*3+ 2] = self.cov_x1_c[i]      
            self.cov_y[i*3+ 2,i*3+ 1] = self.cov_x1_c[i]              
            
        self.salt_parm = np.array([self.mb,self.x1,self.c]).T
#        print len(self.salt_parm), len(self.cov_y), len(self.z), len(self.zcmb)
        return self.salt_parm, self.cov_y, self.z, self.zcmb, self.z_err

    def HD_input_sncosmo_data(self):
        """
        Sample snfit results in HD input
        """

        dico = cPickle.load(open(SUGAR_parameter_pkl))
        self.read_sncosmo()
        self.read_meta()
        self.read_snfit_results()
        Filtre = np.array([True]*len(self.sncosmo_sn_name))
        self.zcmb = []
        self.z_err = []
        for j in range(len(self.sncosmo_sn_name)):
            if self.sncosmo_sn_name[j] in dico.keys():

                for i in range (len(self.meta_sn_name_list)):
                    if self.sncosmo_sn_name[j] == self.meta_sn_name_list[i]:
                        
                        self.z_err.append(self.meta_zhl_err[i])
                        self.zcmb.append(self.meta_zcmb[i])
                        if np.abs(self.sncosmo_x1[j] - self.x1[i]) > 0.01:
                            print 'problem with %s include in sample but big difference between sncosmo and snfit'%(self.sncosmo_sn_name[j])
            else:
                Filtre[j] = False

        for p in dico.keys():
            if p not in self.sncosmo_sn_name:
                print p

        self.sncosmo_x1 = self.sncosmo_x1[Filtre]
        self.sncosmo_x1_err = self.sncosmo_x1_err[Filtre] 
        self.sncosmo_c = self.sncosmo_c[Filtre]
        self.sncosmo_c_err = self.sncosmo_c_err[Filtre]
        self.sncosmo_mb = self.sncosmo_mb[Filtre]
        self.sncosmo_mb_err = self.sncosmo_mb_err[Filtre]
        self.sncosmo_cov_x1_c = self.sncosmo_cov_x1_c[Filtre]
        self.sncosmo_cov_mb_x1 = self.sncosmo_cov_mb_x1[Filtre]
        self.sncosmo_cov_mb_c = self.sncosmo_cov_mb_c[Filtre]
        self.sncosmo_z = self.sncosmo_z[Filtre]
        self.zcmb = np.array(self.zcmb)
        self.z_err = np.array(self.z_err)

        self.sncosmo_cov_y = np.zeros((len(self.sncosmo_mb)*3,len(self.sncosmo_mb)*3))
        
        for i in range (len(self.sncosmo_mb)):
            self.sncosmo_cov_y[i*3,i*3] = self.sncosmo_mb_err[i]**2
            self.sncosmo_cov_y[i*3+ 1,i*3+ 1] = self.sncosmo_x1_err[i]**2
            
            self.sncosmo_cov_y[i*3+ 2,i*3+ 2] = self.sncosmo_c_err[i]**2
            self.sncosmo_cov_y[i*3+ 0,i*3+ 1] = self.sncosmo_cov_mb_x1[i]
            self.sncosmo_cov_y[i*3+ 1,i*3+ 0] = self.sncosmo_cov_mb_x1[i]
            self.sncosmo_cov_y[i*3+ 0,i*3+ 2] = self.sncosmo_cov_mb_c[i]
            self.sncosmo_cov_y[i*3+ 2,i*3+ 0] = self.sncosmo_cov_mb_c[i]
            self.sncosmo_cov_y[i*3+ 1,i*3+ 2] = self.sncosmo_cov_x1_c[i]      
            self.sncosmo_cov_y[i*3+ 2,i*3+ 1] = self.sncosmo_cov_x1_c[i]              
            
        self.salt_parm = np.array([self.sncosmo_mb,self.sncosmo_x1,self.sncosmo_c]).T
        print len(self.salt_parm), len(self.sncosmo_cov_y), len(self.sncosmo_z), len(self.zcmb)
        return self.salt_parm, self.sncosmo_cov_y, self.sncosmo_z, self.zcmb, self.z_err
