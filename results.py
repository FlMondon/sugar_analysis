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
import builtins_SNF as Build_SNF
from scipy import integrate

SUGAR_parameter_pkl = '../sugar/sugar/data_output/data_output/sugar_parameters.pkl'

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]

clight = 299792.458
H0 = 0.000070

class results_snfit():
    def __init__(self, width=10):
        self.width = width
        
    def read_snfit_results(self, snfit_res_path='../sugar_analysis_data/results/results_snfit.txt'):
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
        self.meta_idr = []
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
                self.meta_idr.append(meta[meta_sn_name]['idr.subset'])
        
        self.meta_idr = np.array(self.meta_idr)
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
                    if np.abs(self.mb[i] - self.meta_mb[j]) < 0.0001:
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


        fig = plt.figure(figsize=(8.,8.))       
                
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        f.subplots_adjust(hspace = 0.5)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x0,25,label='$\Delta$ X0')
        ax0_2.hist(self.diff_x0_err,25,label='$\Delta$ X0 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        pdffile = '../sugar_analysis_data/results/x0_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        f.subplots_adjust(hspace = 0.5)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_x1,25,label='$\Delta$ X1')
        ax0_2.hist(self.diff_x1_err,25,label='$\Delta$ X1 error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        pdffile = '../sugar_analysis_data/results/x1_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        
        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        f.subplots_adjust(hspace = 0.5)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_c,25,label='$\Delta$ Color')
        ax0_2.hist(self.diff_c_err,25,label='$\Delta$ Color error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        pdffile = '../sugar_analysis_data/results/color_plot_meta_snfit.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()

        gs = gridspec.GridSpec(2, 1) #subplots ratio
        f, (ax0_1, ax0_2) = plt.subplots(2, sharex=True)
        f.subplots_adjust(hspace = 0.5)
        ax0_1 = plt.subplot(gs[0, 0])
        ax0_2 = plt.subplot(gs[1, 0])
        
        ax0_1.hist(self.diff_mb,50,label='$\Delta$ mb')
        ax0_2.hist(self.diff_mb_err,50,label='$\Delta$ mb error')
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
        
        ax0_1.hist(self.diff_x0_sncosmo,50,label='$\Delta$ x0_'+str(self.width))
        ax0_2.hist(self.diff_x0_err_sncosmo,50,label='$\Delta$ x0 error')
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
        
        ax0_1.hist(self.diff_x1_sncosmo,50,label='$\Delta$ X1_'+str(self.width))
        ax0_2.hist(self.diff_x1_err_sncosmo,50,label='$\Delta$ X1 error')
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
        
        ax0_1.hist(self.diff_c_sncosmo,50,label='$\Delta$ Color_'+str(self.width))
        ax0_2.hist(self.diff_c_err_sncosmo,50,label='$\Delta$ Color error')
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
        
        ax0_1.hist(self.diff_mb_sncosmo,50,label='$\Delta$ mb_'+str(self.width))
        ax0_2.hist(self.diff_mb_err_sncosmo,50,label='$\Delta$ mb error')
        ax0_1.legend()
        ax0_2.legend()
        ax0_1.set_ylabel('N')
        ax0_2.set_ylabel('N')
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        pdffile = '../sugar_analysis_data/results/mb_plot_'+str(self.width)+'.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()

        plt.hist(self.diff_chi2,50,label='$\Delta$ chi2_'+str(self.width))
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
        plt.hist(diff_mb_sncosmo,40,label='SNCOSMO', range = (-0.02,0.015))
        pdffile = '../sugar_analysis_data/results/diff_mb_th1A_'+str(self.width)+'_sncosmo.pdf'

#        plt.hist(diff_mb_sncosmo,40,label='SNCOSMO',range = (-0.002,0.0015))
        plt.text(-0.018, 80, '\mu = %f \n'%(mean_hist)+' \sigma = %f'%(sig_hist), fontweight = 'bold', fontsize = 20)

        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        plt.legend()
#        pdffile = '../sugar_analysis_data/results/diff_mb_th_5A_1A_sncosmo.pdf' 
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show() 
        
    def diff_snfit_dL(self):
        diff_mb = [] 
        outlayers = ['SN2012cu','SN2009hi','SNF20061022-014','SNNGC7589']
        self.read_snfit_results()
        snfit_sn_name_5A = self.sn_name
        snfit_mb_5A = self.mb
#        self.read_sncosmo()
#        snfit_sn_name_5A = self.sncosmo_sn_name
#        snfit_mb_5A = self.sncosmo_mb
        
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
        plt.hist(diff_mb,60,label='SNFIT',range = (-0.02,0.015))
        pdffile = '../sugar_analysis_data/results/diff_mb_th1A_'+str(self.width)+'_snfit.pdf'
        
#        plt.hist(diff_mb,40, range = (-0.002,0.0015),label='SNFIT')
        plt.text(-0.018, 25, '\mu = %f \n'%(mean_hist)+' \sigma = %f'%(sig_hist), fontweight = 'bold', fontsize = 20)

#        pdffile = '../sugar_analysis_data/results/diff_mb_th_5A_1A_snfit.pdf'   
#        plt.xlim(-0.002,0.0015)
        plt.ticklabel_format(axis='x', style='scientific', scilimits=(-1, 2))
        plt.legend()
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show() 

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
        plt.plot([-5, 50], [0,0],label='Sncosmo 1A')
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
#        return self.salt_parm, self.cov_y, self.z, self.meta_zcmb, self.meta_zhl_err, self.sn_name, self.meta_idr
        return self.salt_parm, self.cov_y, self.z, self.zcmb, self.z_err

    def HD_input_sncosmo_data(self, sn_list):
        """
        Sample snfit results in HD input
        """

        dico = cPickle.load(open(SUGAR_parameter_pkl))
        self.read_sncosmo(path='../sugar_analysis_data/results/res_salt2_SNF_5_nomodelcov.txt')
        self.read_meta()
        self.read_snfit_results()
        Filtre = np.array([True]*len(self.sncosmo_sn_name))
        self.zcmb = []
        self.z_err = []
        for j, sn_name in enumerate(self.sncosmo_sn_name):
#            if self.sncosmo_sn_name[j] in dico.keys():
#
#                for i in range (len(self.meta_sn_name_list)):
#                    if self.sncosmo_sn_name[j] == self.meta_sn_name_list[i]:
#                        
#                        self.z_err.append(self.meta_zhl_err[i])
#                        self.zcmb.append(self.meta_zcmb[i])
#                        if np.abs(self.sncosmo_x1[j] - self.x1[i]) > 0.01:
#           i                 print 'problem with %s include in sample but big difference between sncosmo and snfit'%(self.sncosmo_sn_name[j])
#            else:
#                 Filtre[j] = False
            if sn_name in sn_list:
                Filtre[j] = True
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


class Sugar_photometry_plot():

  
        
    def read_sugar_phot(self, path='../sugar_analysis_data/results/res_sugar_SNF.txt'):
        self.sug_phot_res = np.loadtxt(path ,dtype='str')
        self.sug_phot_sn_name = np.array(self.sug_phot_res[:,0],str)
        self.sug_phot_z = np.array(self.sug_phot_res[:,2],float)
        self.sug_phot_Mgr = np.array(self.sug_phot_res[:,4],float)
        self.sug_phot_Mgr_err = np.array(self.sug_phot_res[:,5],float)
        self.sug_phot_q1 = np.array(self.sug_phot_res[:,6],float)
        self.sug_phot_q1_err = np.array(self.sug_phot_res[:,7],float)
        self.sug_phot_q2 = np.array(self.sug_phot_res[:,8],float)
        self.sug_phot_q2_err = np.array(self.sug_phot_res[:,9],float)
        self.sug_phot_q3 = np.array(self.sug_phot_res[:,10],float)
        self.sug_phot_q3_err = np.array(self.sug_phot_res[:,11],float)
        self.sug_phot_A = np.array(self.sug_phot_res[:,12],float)
        self.sug_phot_A_err = np.array(self.sug_phot_res[:,13],float)
#        self.sug_phot_t0 = np.array(self.sug_phot_res[:,24],float)
#        print self.sug_phot_t0[0]
#        self.sug_phot_t0_err = np.array(self.sug_phot_res[:,25],float)
#        self.sug_phot_chi2 = np.array(self.sug_phot_res[:,26],float)
        self.cov_mgrey_q1 = np.array(self.sug_phot_res[:,14],float)
        self.cov_mgrey_q2 = np.array(self.sug_phot_res[:,15],float)
        self.cov_mgrey_q3 = np.array(self.sug_phot_res[:,16],float)
        self.cov_mgrey_A =  np.array(self.sug_phot_res[:,17],float)
        self.cov_q1_q2 =  np.array(self.sug_phot_res[:,18],float)
        self.cov_q1_q3 = np.array(self.sug_phot_res[:,19],float)
        self.cov_q1_A = np.array(self.sug_phot_res[:,20],float)
        self.cov_q2_q3 = np.array(self.sug_phot_res[:,21],float)
        self.cov_q2_A = np.array(self.sug_phot_res[:,22],float)
        self.cov_q3_A = np.array(self.sug_phot_res[:,23],float)
        
    def read_sugar_phot_t0_fix(self, path='../sugar_analysis_data/results/res_sugar_SNF_t0_fix.txt'):       
        self.sug_phot_res_t0_fix = np.loadtxt(path ,dtype='str')
        self.sug_phot_q2_t0_fix = np.array(self.sug_phot_res[:,8],float)
        self.sug_phot_sn_name_t0_fix = np.array(self.sug_phot_res_t0_fix[:,0],str)
        self.sug_phot_chi2_t0_fix = np.array(self.sug_phot_res_t0_fix[:,24],float)
    
    def plot_chi2_t0_fix(self):
        '''
        '''
        self.read_sugar_phot_t0_fix()
        self.read_sugar_phot()
        SUGAR_parameter_pkl = '../sugar_model/sugar_parameters.pkl'
        dico = cPickle.load(open(SUGAR_parameter_pkl))
        self.diff_chi2_t0_fix = []
        for i in range(len(self.sug_phot_sn_name_t0_fix)):
            for j in range(len(self.sug_phot_sn_name)):
                if self.sug_phot_sn_name_t0_fix[i] == self.sug_phot_sn_name[j] and self.sug_phot_sn_name_t0_fix[i] in dico.keys():
                    self.diff_chi2_t0_fix.append(self.sug_phot_chi2_t0_fix[i]-self.sug_phot_chi2[j])
        plt.hist(self.diff_chi2_t0_fix,50,label='$\Delta$ chi2')
        pdffile = '../sugar_analysis_data/results/diff_chi2_t0_fix.pdf'
        plt.legend()
        plt.savefig(pdffile, bbox_inches='tight')        
        plt.legend()
        plt.show()
        plt.close()
                
                
    def read_sugar_spectro(self, path='../sugar_model/sugar_parameters.pkl'):

        SUGAR_parameter_pkl = path
        self.meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        self.dico = cPickle.load(open(SUGAR_parameter_pkl)) 
        self.sug_spec_sn_name = []
        self.sug_spec_Mgr = []
        self.sug_spec_q1 =  []
        self.sug_spec_q2 =  []
        self.sug_spec_q3 =  []
        self.sug_spec_A =  []
        self.sug_spec_t0 =  []
        self.sug_spec_cov_q =  []
        
        
        for sn_name in self.dico.keys():
            self.sug_spec_sn_name.append(sn_name)
            self.sug_spec_Mgr.append(self.dico[sn_name]['grey']+ self.distance_modulus_th(self.meta[sn_name]['host.zhelio'],self.meta[sn_name]['host.zcmb']))
            self.sug_spec_q1.append(self.dico[sn_name]['q1'])
            self.sug_spec_q2.append(self.dico[sn_name]['q2'])
            self.sug_spec_q3.append(self.dico[sn_name]['q3'])
            self.sug_spec_A.append(self.dico[sn_name]['Av'])
            self.sug_spec_t0.append(self.meta[sn_name]['salt2.DayMax'])
            self.sug_spec_cov_q.append(self.dico[sn_name]['cov_q']) 
    
    def errors_list(self):
        self.read_sugar_spectro()
        self.sug_spec_Mgr_err = []
        self.sug_spec_q1_err =  []
        self.sug_spec_q2_err =  []
        self.sug_spec_q3_err =  []
        self.sug_spec_A_err =  []

        for i in range(len(self.sug_spec_sn_name)):
            self.sug_spec_Mgr_err.append(self.sug_spec_cov_q[i][0,0])
            self.sug_spec_q1_err.append(self.sug_spec_cov_q[i][1,1])
            self.sug_spec_q2_err.append(self.sug_spec_cov_q[i][2,2])
            self.sug_spec_q3_err.append(self.sug_spec_cov_q[i][3,3])
            self.sug_spec_A_err.append(self.sug_spec_cov_q[i][4,4])
        
        
    def int_cosmo(self, z, Omega_M=0.3):     
        return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
        
    def luminosity_distance(self, zhl,zcmb):           
        integr = integrate.quad(self.int_cosmo, 0, zcmb)[0]
    
        return (1+zhl)*(clight/H0)*integr
     
    def distance_modulus_th(self, zhl,zcmb):      
        return 5.*np.log(self.luminosity_distance(zhl,zcmb))/np.log(10.)-5. 
      
    def plot_phot_spec(self):
        self.read_sugar_phot()
        self.read_sugar_spectro()
        self.errors_list()
        SUGAR_parameter_pkl = '../sugar_model/sugar_parameters.pkl'
        dico = cPickle.load(open(SUGAR_parameter_pkl))
        x_list = np.linspace(-10,10,5)
        self.diff_Mgr = []
        self.diff_q1 = [] 
        self.diff_q2 = [] 
        self.diff_q3 = [] 
        self.diff_A = [] 
        self.diff_t0 = [] 
        self.scatter_spec_Mgr = []
        self.scatter_spec_q1 = [] 
        self.scatter_spec_q2 = [] 
        self.scatter_spec_q3 = [] 
        self.scatter_spec_A = []
        self.scatter_spec_t0 = []
        self.scatter_phot_Mgr = []
        self.scatter_phot_q1 = [] 
        self.scatter_phot_q2 = [] 
        self.scatter_phot_q3 = [] 
        self.scatter_phot_A = [] 
        self.scatter_spec_Mgr_err = []
        self.scatter_spec_q1_err = [] 
        self.scatter_spec_q2_err = [] 
        self.scatter_spec_q3_err = [] 
        self.scatter_spec_A_err = []
        self.scatter_phot_Mgr_err = []
        self.scatter_phot_q1_err = [] 
        self.scatter_phot_q2_err = [] 
        self.scatter_phot_q3_err = [] 
        self.scatter_phot_A_err = [] 
        self.scatter_phot_t0 = []

        
        
        for i in range(len(self.sug_spec_sn_name)):
            for j in range(len(self.sug_phot_sn_name)):
                if self.sug_spec_sn_name[i] == self.sug_phot_sn_name[j] and self.sug_spec_sn_name[i] in dico.keys():
#                    if self.sug_spec_sn_name[i] == 'SN2008ec':
#                        print self.sug_spec_q2[i], self.sug_phot_q2[j]
                    self.diff_q1.append(self.sug_spec_q1[i]-self.sug_phot_q1[j])
                    self.diff_q2.append(self.sug_spec_q2[i]-self.sug_phot_q2[j])
                    self.diff_q3.append(self.sug_spec_q3[i]-self.sug_phot_q3[j])
                    self.diff_A.append(self.sug_spec_A[i]-self.sug_phot_A[j])
                    self.diff_Mgr.append(self.sug_spec_Mgr[i]-self.sug_phot_Mgr[j])
                    if self.sug_spec_q2[i]-self.sug_phot_q2[j] > 1.:
                        print self.sug_spec_sn_name[i]
                    self.diff_t0.append(self.sug_spec_t0[i]-self.sug_phot_t0[j])
                    self.scatter_spec_t0.append(self.sug_spec_t0[i])
                    self.scatter_spec_Mgr.append(self.sug_spec_Mgr[i])
                    self.scatter_spec_q1.append(self.sug_spec_q1[i]) 
                    self.scatter_spec_q2.append(self.sug_spec_q2[i]) 
                    self.scatter_spec_q3.append(self.sug_spec_q3[i]) 
                    self.scatter_spec_A.append(self.sug_spec_A[i]) 
                    self.scatter_phot_Mgr.append(self.sug_phot_Mgr[j])
                    self.scatter_phot_q1.append(self.sug_phot_q1[j]) 
                    self.scatter_phot_q2.append(self.sug_phot_q2[j]) 
                    self.scatter_phot_q3.append(self.sug_phot_q3[j]) 
                    self.scatter_phot_A.append(self.sug_phot_A[j]) 
                    self.scatter_phot_t0.append(self.sug_phot_t0[j])
                    self.scatter_spec_Mgr_err.append(self.sug_spec_Mgr_err[i])
                    self.scatter_spec_q1_err.append(self.sug_spec_q1_err[i]) 
                    self.scatter_spec_q2_err.append(self.sug_spec_q2_err[i]) 
                    self.scatter_spec_q3_err.append(self.sug_spec_q3_err[i]) 
                    self.scatter_spec_A_err.append(self.sug_spec_A_err[i]) 

                    self.scatter_phot_Mgr_err.append(self.sug_phot_Mgr_err[j])
                    self.scatter_phot_q1_err.append(self.sug_phot_q1_err[j]) 
                    self.scatter_phot_q2_err.append(self.sug_phot_q2_err[j]) 
                    self.scatter_phot_q3_err.append(self.sug_phot_q3_err[j]) 
                    self.scatter_phot_A_err.append(self.sug_phot_A_err[j])
                    
                    if abs(self.sug_spec_Mgr[i]-self.sug_phot_Mgr[j]) >= 0.2:
                        print self.sug_spec_sn_name[i]
        print len(self.scatter_phot_Mgr)
        for sn_name in  dico.keys():
            if sn_name not in self.sug_phot_sn_name:
                print sn_name
       
        #plot setup
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
        
        #plot q1 hist & scatter    
        plt.hist(self.diff_q1,50,label='$\Delta$ q1')
        pdffile = '../sugar_analysis_data/results/diff_q1.pdf'
        plt.savefig(pdffile, bbox_inches='tight')        
        plt.legend()
        plt.show()
        plt.close()
        
        plt.errorbar(self.scatter_spec_q1,self.scatter_phot_q1,xerr=self.scatter_spec_q1_err,yerr=self.scatter_phot_q1_err, color='red', fmt='.', mfc='red', zorder=1)
        plt.ylabel('phot',fontsize=25)
        plt.xlabel('spec',fontsize=25)
        plt.figtext(0.2, 0.8, 'q1', fontsize=25)
        pdffile = '../sugar_analysis_data/results/scatter_q1.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.legend()
        plt.show()
        plt.close()
        
        #plot q2 hist & scatter
        plt.hist(self.diff_q2,50,label='$\Delta$ q2')
        pdffile = '../sugar_analysis_data/results/diff_q2.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.legend()
        plt.show()
        plt.close()
 
        plt.errorbar(self.scatter_spec_q2,self.scatter_phot_q2,xerr=self.scatter_spec_q2_err,yerr=self.scatter_phot_q2_err, color='red', fmt='.', mfc='red', zorder=1)
        plt.plot(x_list,x_list)
        plt.xlim(-10,10)
        plt.ylim(-10,10)
        plt.ylabel('phot',fontsize=25)
        plt.xlabel('spec',fontsize=25)
        plt.figtext(0.2, 0.8, 'q2', fontsize=25)
        pdffile = '../sugar_analysis_data/results/scatter_q2.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.legend()
        plt.show()
        plt.close()
        
        #plot q3 hist & scatter
        plt.hist(self.diff_q3,50,label='$\Delta$ q3')
        pdffile = '../sugar_analysis_data/results/diff_q3.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.legend()
        plt.show()
        plt.close()

        plt.errorbar(self.scatter_spec_q3,self.scatter_phot_q3,xerr=self.scatter_spec_q3_err,yerr=self.scatter_phot_q3_err, color='red', fmt='.', mfc='red', zorder=1)
        plt.ylabel('phot',fontsize=25)
        plt.xlabel('spec',fontsize=25)
        plt.figtext(0.2, 0.8, 'q3', fontsize=25)
        pdffile = '../sugar_analysis_data/results/scatter_q3.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.legend()
        plt.show()
        plt.close()
        
        #plot A hist & scatter
        plt.hist(self.diff_A,50,label='$\Delta$ A')
        pdffile = '../sugar_analysis_data/results/diff_A.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.legend()
        plt.show()
        plt.close()

        plt.errorbar(self.scatter_spec_A,self.scatter_phot_A,xerr=self.scatter_spec_A_err,yerr=self.scatter_phot_A_err, color='red', fmt='.', mfc='red', zorder=1)
        plt.ylabel('phot',fontsize=25)
        plt.xlabel('spec',fontsize=25)
        plt.figtext(0.2, 0.8, 'A', fontsize=25)
        pdffile = '../sugar_analysis_data/results/scatter_A.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.legend()
        plt.show()
        plt.close()
        
        #plot Mgr hist & scatter
        plt.hist(self.diff_Mgr,50,label='$\Delta$ Mgr')
        pdffile = '../sugar_analysis_data/results/diff_Mgr.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.legend()
        plt.show()
        plt.close()    
        
        plt.errorbar(self.scatter_spec_Mgr,self.scatter_phot_Mgr,xerr=self.scatter_spec_Mgr_err,yerr=self.scatter_phot_Mgr_err, color='red', fmt='.', mfc='red', zorder=1)
        plt.ylabel('phot',fontsize=25)
        plt.xlabel('spec',fontsize=25)
        plt.figtext(0.2, 0.8, 'Mgr', fontsize=25)
        pdffile = '../sugar_analysis_data/results/scatter_Mgr.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.legend()
        plt.show()
        plt.close()
        
        #plot t0 scatter
        plt.hist(self.diff_t0,50,label='$\Delta$ t0')
        pdffile = '../sugar_analysis_data/results/diff_t0.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.legend()
        plt.show()
        plt.close()
        plt.scatter(self.scatter_spec_t0, self.scatter_phot_t0, c='r', marker='.')
        pdffile = '../sugar_analysis_data/results/scatter_t0.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        plt.close()
        
        #t0_q2
        plt.scatter(self.diff_q2, self.diff_t0, c='r', marker='.')
        print np.corrcoef(self.diff_q2,self.diff_t0)
        pdffile = '../sugar_analysis_data/results/diff_t0vsq2.pdf'
        plt.ylabel('$\Delta$ t0',fontsize=25)
        plt.xlabel('$\Delta$ q2',fontsize=25)
        plt.savefig(pdffile, bbox_inches='tight')
        plt.show()
        plt.close()
        
    def HD_input_sugar(self, path= '../sugar_analysis_data/results/res_sugar_SNF.txt'):
        SUGAR_parameter_pkl = '../sugar_model/sugar_parameters.pkl'
        dico = cPickle.load(open(SUGAR_parameter_pkl))     
        meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        self.read_sugar_phot(path=path)       
        self.HD_Mgr = []
        self.HD_q1 = [] 
        self.HD_q2 = [] 
        self.HD_q3 = [] 
        self.HD_A = []           
        self.zhl = []
        self.zcmb = []
        self.zerr = []
#        self.HD_cov = np.zeros([(len(dico.keys())-2)*5, (len(dico.keys())-2)*5])
        self.HD_cov = np.zeros([len(self.sug_phot_sn_name)*5, (len(self.sug_phot_sn_name)*5)])
#        self.HD_cov_Mgr = np.zeros([len(dico.keys())-2, len(dico.keys())-2])
        self.HD_cov_Mgr = np.zeros([len(self.sug_phot_sn_name), len(self.sug_phot_sn_name)])
        j = 0
        print len(self.sug_phot_sn_name)
        for i in range(len(self.sug_phot_sn_name)):  
#            if self.sug_phot_sn_name[i] in dico.keys() and self.sug_phot_sn_name[i] not in ['SNF20070331-025','SNF20080905-005']:
                
                self.HD_Mgr.append(self.sug_phot_Mgr[i])
                self.HD_q1.append(self.sug_phot_q1[i]) 
                self.HD_q2.append(self.sug_phot_q2[i]) 
                self.HD_q3.append(self.sug_phot_q3[i]) 
                self.HD_A.append(self.sug_phot_A[i])
                
                self.HD_cov[j*5,j*5] = self.sug_phot_Mgr_err[i]**2
                self.HD_cov[j*5 +1,j*5] = self.cov_mgrey_q1[i]
                self.HD_cov[j*5 +2,j*5] = self.cov_mgrey_q2[i]
                self.HD_cov[j*5 +3,j*5] = self.cov_mgrey_q3[i]
                self.HD_cov[j*5 +4,j*5] = self.cov_mgrey_A[i]
                self.HD_cov[j*5,j*5 +1] = self.cov_mgrey_q1[i]
                self.HD_cov[j*5 +1,j*5 +1] = self.sug_phot_q1_err[i]**2
                self.HD_cov[j*5 +2,j*5 +1] = self.cov_q1_q2[i]
                self.HD_cov[j*5 +3,j*5 +1] = self.cov_q1_q3[i]
                self.HD_cov[j*5 +4,j*5 +1] = self.cov_q1_A[i]
                self.HD_cov[j*5 ,j*5 +2] = self.cov_mgrey_q2[i]
                self.HD_cov[j*5 +1,j*5 +2] =  self.cov_q1_q2[i]
                self.HD_cov[j*5 +2,j*5 +2] = self.sug_phot_q2_err[i]**2
                self.HD_cov[j*5 +3,j*5 +2] = self.cov_q2_q3[i]
                self.HD_cov[j*5 +4,j*5 +2] = self.cov_q2_A[i]
                self.HD_cov[j*5,j*5 +3] = self.cov_mgrey_q3[i]
                self.HD_cov[j*5 +1,j*5 +3] = self.cov_q1_q3[i]
                self.HD_cov[j*5 +2,j*5 +3] = self.cov_q2_q3[i]
                self.HD_cov[j*5 +3,j*5 +3] = self.sug_phot_q3_err[i]**2
                self.HD_cov[j*5 +4,j*5 +3] = self.cov_q3_A[i]
                self.HD_cov[j*5,j*5 +4] = self.cov_mgrey_A[i]
                self.HD_cov[j*5 +1,j*5 +4] = self.cov_q1_A[i]
                self.HD_cov[j*5 +2,j*5 +4] = self.cov_q2_A[i]
                self.HD_cov[j*5 +3,j*5 +4] = self.cov_q3_A[i]
                self.HD_cov[j*5 +4,j*5 +4] = self.sug_phot_A_err[i]**2    
                
                self.HD_cov_Mgr[j,j] = self.sug_phot_Mgr_err[i]**2
#                print self.HD_cov_Mgr[j,j]
                self.zhl.append(meta[self.sug_phot_sn_name[i]]['host.zhelio'])
                self.zcmb.append(meta[self.sug_phot_sn_name[i]]['host.zcmb'])
                self.zerr.append(meta[self.sug_phot_sn_name[i]]['host.zhelio.err'])
                j +=1
        
        print len (self.zhl)
        print len (self.HD_Mgr)

        self.sug_param = np.array([np.array(self.HD_Mgr), np.array(self.HD_q1), np.array(self.HD_q2), np.array(self.HD_q3), np.array(self.HD_A)]).T
        self.zhl = np.array(self.zhl)
        self.zcmb = np.array(self.zcmb)
        self.zerr = np.array(self.zerr)
        
    def plot_errors(self):
        """
        """
        self.read_sugar_phot( path='../sugar_analysis_data/results/res_sugar_SNF_t0_fix.txt')
        SUGAR_parameter_pkl = '../sugar_model/sugar_parameters.pkl'
        dico = cPickle.load(open(SUGAR_parameter_pkl))       
        mask = []
        for i in range(len(self.sug_phot_sn_name)):
           if self.sug_phot_sn_name[i] in dico.keys(): 
               mask.append(True)
           else:
               mask.append(False)
        mask = np.array(mask)
        self.phot_Mgr_err = self.sug_phot_Mgr_err[mask]
        self.phot_q1_err = self.sug_phot_q1_err[mask]
        self.phot_q2_err = self.sug_phot_q2_err[mask]
        self.phot_q3_err = self.sug_phot_q3_err[mask]
        self.phot_A_err = self.sug_phot_A_err[mask]
        plt.hist(self.phot_q2_err,50,label='q2 err')