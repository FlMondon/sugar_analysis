#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 13:02:35 2018

@author: florian
"""

import numpy as np
p2 = '../sugar_model/'
p3 = '../sncosmo_jla/salt2-4/'
from operator import itemgetter
import sncosmo
import builtins_SNF as Build_SNF
from sncosmo.salt2utils import BicubicInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib import pyplot as plt
import cPickle as pkl
from matplotlib import rc, rcParams
import copy
import os



CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27       # [erg s]
# wavelength limits for sugar model
wl_min_sug = 3254.01639
wl_max_sug = 8649.03871
wl_min_salt = 2000.000000
wl_max_salt = 9200.000000
t_min_sug = -12
t_max_sug = 48

from astropy.table import Table
try:
    Build_SNF.register_SNf_bands_width(width=10)
    Build_SNF.mag_sys_SNF_width(width=10)
    Build_SNF.register_SUGAR()
except:
    print 'Filters and mag sys already registred'
class sugar_spectrum():
    
    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
                    m2file='sugar_template_2.dat',
                    m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat', 
                 parameters_init=np.array([0., 0., 0., 0., 37.]),
                 name=None, version=None):
        

        self._SCALE_FACTOR = 1.
        self._param_names = ['q1', 'q2', 'q3', 'A', 'Mgr']
        
    
        infile = open(p2 + 'SUGAR_model_v1.asci', 'r')
        lines = infile.readlines()
        infile.close()
        
        listik = []
        for line in lines:
            new_line = [float(i) for i in line.split()]
            listik.append(new_line)
        
        s = sorted(listik, key=itemgetter(0))
        
        names = ['0','1','2','3','4']
        for i in names:
            outfile = open(p2 + 'sugar_template_' + i + '.dat', 'w')
            for line in s:
                j = 2+int(i)
                outfile.write('%4.4f %8.8f %8.8f' %(line[0],line[1],line[j]))
                outfile.write('\n')
            outfile.close()        
            
        
        self.step = 15
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = parameters_init
        self._salt_params = [1., 0., 0.]
        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(p2 + names_or_objs[key])
            values *= self._SCALE_FACTOR
            # self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)
            self._model[key] = BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave
        self.dic_sts = None
    def model_spectrum_flux_m0(self, phase, wave):
        '''
        '''
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
#        return 10. ** (-0.4 * (m0 + 48.59)) / (wave** 2 / 299792458. * 1.e-10)
        return 10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4]  + 48.59))/ (wave** 2 / 299792458. * 1.e-10)


#        source_salt2 = sncosmo.get_source('salt2')
#        flux_salt2 = source_salt2._flux(phase, wave)
#        return flux_salt2  

    def spectrum_generator(self, parameters, phase=None, mean_errors=None):
        """
        simulate sugar spectrum given a set of parameters
        parameters : array of sugar parameters (np.array([q1, q2, q3, A, Mgr]))
        mean_error : array of errors mean for different phase 
        """
        self.dic_spectrum = {}
        par_init = copy.deepcopy(self._parameters)
        if  not isinstance(phase, np.ndarray):
            phase = self._phase
        self._parameters = parameters
        wave = np.linspace(3000,9000,10000)
        flux = self.model_spectrum_flux_m0(phase, wave)
        flux_err = np.zeros_like(flux)
        self.dic_spectrum['parameters'] = parameters
        self.dic_spectrum['wave'] = wave
        self.dic_spectrum['phase'] = phase
        if not isinstance(mean_errors,np.ndarray):
            self.dic_spectrum['fluxerr'] = None
            self.dic_spectrum['flux'] = flux
        else:
            for i in range(len(flux)):
                flux_err[i]=np.ones(len(flux[i]))*mean_errors[i]
                for j in range(len(flux[i])):
                    flux[i,j] = np.random.normal(flux[i,j],mean_errors[i],1)[0]
            self.dic_spectrum['fluxerr'] = flux_err 
            self.dic_spectrum['flux'] = flux            
        self._parameters = par_init 
        
    def error_generator_lc_factory(self, band):
        """
        """
        if band == 'new_fI_10' or band =='fI_10':
            error = np.abs(np.random.normal(0.0101, 0.0425, 1)[0])
        if band == 'fU_10':
            error = np.abs(np.random.normal(0.0146, 0.0741, 1)[0])
        if band == 'fB_10':
            error = np.abs(np.random.normal(0.0141, 0.06898, 1)[0])
        if band == 'fV_10':
            error = np.abs(np.random.normal(0.0139, 0.0646, 1)[0])
        if band == 'fR_10':
            error = np.abs(np.random.normal(0.0088, 0.0395, 1)[0])
        return error
    
    def AstropyTable_flux(self, error=False):
        """
        """
    
                
        band_used = ['new_fI_10','fB_10','fV_10','fR_10','fU_10']    
        time = []
        flux = []
        fluxerr = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')
        for i in range(len(self._phase)):
            for b in band_used:
            
                filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+b+'.dat')
                wlen = filt2[:,0]
                tran = filt2[:,1]
                self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
                
                #computation of the integral
                dt = 100000
                xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
                self.xs = xs
                dxs = (float(wl_max_sug-wl_min_sug)/(dt-1))
                spec_flux = self.model_spectrum_flux_m0(self._phase, xs)
#                plt.plot(xs,self.splB(xs))
#                plt.plot(xs,spec_flux[i,:])
#                print spec_flux[i,:] 
                plt.show()
                inte = np.sum(spec_flux[i,:]*(xs  / (CLIGHT*HPLANCK))*self.splB(xs)*dxs)
                    
                flux.append(inte)
                if error:
                    fluxerr.append(self.error_generator_lc_factory(b))
                else:
                    fluxerr.append(0.00000000001)
                time.append(self._phase[i])
                band.append(b)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
        data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
        return data

    def AstropyTable_flux_with_sim_spectra(self):
        """
        """
    
                
        band_used = ['new_fI_10','fB_10','fV_10','fR_10','fU_10']    
        time = []
        flux = []
        fluxerr = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')
        for i in range(len(self._phase)):
            for b in band_used:
            
                filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+b+'.dat')
                wlen = filt2[:,0]
                tran = filt2[:,1]
                self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
                
                #computation of the integral
                dt = 100000
                xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
                self.xs = xs
                dxs = (float(wl_max_sug-wl_min_sug)/(dt-1))
                
    def fit_lc_sugar(self, data ):
        """
        """
        source = sncosmo.get_source('sugar')
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
        model.set(mwebv=0)
        model.set(z=0)
                   
        #initial iteration x1 fix
        print 'initialisation'            

        model.set(q1=0.0)
        model.set(q2=0.0)
        model.set(q3=0.0)     
        model.set(t0=0.1)
        res, fitted_model = sncosmo.fit_lc(data, model, ['A','Mgr'], modelcov = False)

            
        print 'first iteration'
        chi2 = res.chisq
        print chi2
        chi2p = chi2*2
        m=0
    
        while chi2 < chi2p and m < 10:

            print m
            if m > 0:
                resp = res
                fitted_modelp = fitted_model
            m += 1
            t_peak = fitted_model.parameters[1]
            #print t_peak,fitted_model.parameters[4]

            t1 = t_peak + t_min_sug*(1 + model.get('z'))
            t2 = t_peak + t_max_sug*(1 + model.get('z'))
                        
            A=[]
            data_new = copy.deepcopy(data)
            for i in range(len(data_new)):                    
                if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                    A.append(i)
            A=np.array(A)
            for i in range(len(A)):
                #print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2)
                data_new.remove_row(A[i])
                A-=1   



            print "WARNING: Sugar don't have model covariance for the moment"
            modelcov = False
            model.set(q1=res.parameters[2])
            model.set(q2=res.parameters[3])
            model.set(q3=res.parameters[4])
            model.set(A=res.parameters[5])
            model.set(Mgr=res.parameters[6])
            #print res.parameters[6]
            model.set(t0=res.parameters[1])                        
            res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0','q1', 'q2', 'q3', 'A', 'Mgr'], modelcov  = modelcov)
            

#                        sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                    
#                        print res.chisq
    #                    print res

            chi2p = chi2
            chi2 = res.chisq
            print chi2p, chi2
        #final results
        res = resp
        fitted_model = fitted_modelp
        print res    
    def fit_lc(self, show_lc=True):
        '''
        '''

        data = self.AstropyTable_flux()
        
        source = sncosmo.get_source('salt2')
        source.EBV_snfit = 0.
        source.z_snfit = 0.
        source.Rv_snfit = 3.1
        dust = sncosmo.CCM89Dust()
       
        model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
        model.set(mwebv=0.)
        model.set(z=0.)
        model.set(t0=0.)
        res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0','x1','c'], modelcov = True)
        if show_lc :
            sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
        return res
#                print self.flux
    
    def sugar_dep_salt(self):
        """
        """
        x0 = []
        x0_err = []
        x1 = [] 
        x1_err = []
        color = [] 
        color_err = []
        t0 = []
        t0_err = []
        self.dic_sts = {}
        for p in range(len(self._parameters)):
            if p == 4:
                 param = np.linspace(35,45,self.step)
            else:
                 param = np.linspace(-4,4,self.step)
            for i in param:
                 print i
                 self._parameters[p] = i
                 res = self.fit_lc(show_lc=False)
                 x0.append(-2.5*np.log10(res['parameters'][2]))
                 t0.append(res['parameters'][1])
                 x1.append(res['parameters'][3])
                 color.append(res['parameters'][4])
                 x0_err.append(res['errors']['x0']*-2.5/(np.log(10)*res['parameters'][2]))
                 x1_err.append(res['errors']['x1'])
                 t0_err.append(res['errors']['t0'])
                 color_err.append(res['errors']['c'])
            self.dic_sts['x0_'+str(p)] = x0
            self.dic_sts['t0_'+str(p)] = t0
            self.dic_sts['x1_'+str(p)] = x1
            self.dic_sts['color_'+str(p)] = color
            self.dic_sts['x0_err_'+str(p)] = x0_err
            self.dic_sts['t0_err'+str(p)] = t0_err
            self.dic_sts['x1_err_'+str(p)] = x1_err
            self.dic_sts['color_err_'+str(p)] = color_err
            self._parameters = np.array([0., 0., 0., 0., 0.])
            x0 = []
            x0_err = []
            x1 = [] 
            x1_err = []
            color = [] 
            color_err = []
            t0 = []
            t0_err = []
        dic_file = open('../sugar_analysis_data/results/sug_to_salt.pkl','w')
        pkl.dump(self.dic_sts,dic_file)
        
    def plot_sds(self):
        """
        """
        if self.dic_sts == None:
             self.dic_sts = pkl.load(open('../sugar_analysis_data/results/sug_to_salt.pkl'))  
        # Plot Setup
        rcParams['font.size'] = 6.
        font = {'family': 'normal', 'size': 5}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=15)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)        
        
        for p in range(5):
            if p == 4:
                 param = np.linspace(35,45,self.step)
            else:
                 param = np.linspace(-4,4,self.step)
            
            plt.subplot(4,5,p+1)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['x0_'+str(p)],yerr=self.dic_sts['x0_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
#            print self.dic_sts['x0_'+str(p)]
            plt.scatter(param, self.dic_sts['x0_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('-2.5 log10(x0)',fontsize=15)
#            plt.xlabel(self._param_names[p],fontsize=25)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/x0_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
            
            plt.subplot(4,5,p+6)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['x1_'+str(p)],yerr=self.dic_sts['x1_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, self.dic_sts['x1_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('x1',fontsize=15)
#            plt.xlabel(self._param_names[p],fontsize=25)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/x1_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
            
            plt.subplot(4,5,p+11)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['color_'+str(p)],yerr=self.dic_sts['color_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, self.dic_sts['color_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('color',fontsize=15)
            
            plt.subplot(4,5,p+16)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['t0_'+str(p)],yerr=self.dic_sts['color_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, self.dic_sts['t0_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('t0',fontsize=15)
                
            plt.xlabel(self._param_names[p],fontsize=15)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/color_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
        pdffile = '../sugar_analysis_data/results/mapping_sug_to_salt.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.show()
        plt.close()
        

                 
#        param = np.linspace(-4,4,self.step) 
#        plt.scatter(param, self.dic_sts['t0_1'], color='red', marker='.', zorder=1)
#        pdffile = '../sugar_analysis_data/results/t0_q2.pdf'
#        plt.savefig(pdffile, bbox_inches='tight') 
#        plt.show()
        
    def plot_sds_t0_fix(self):
        """
        """
        if self.dic_sts == None:
             self.dic_sts = pkl.load(open('../sugar_analysis_data/results/sug_to_salt.pkl'))  
        # Plot Setup
        rcParams['font.size'] = 6.
        font = {'family': 'normal', 'size': 5}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=15)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)        
        
        for p in range(5):
            if p == 4:
                 param = np.linspace(35,45,self.step)
            else:
                 param = np.linspace(-4,4,self.step)
            
            plt.subplot(3,5,p+1)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['x0_'+str(p)],yerr=self.dic_sts['x0_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
#            print self.dic_sts['x0_'+str(p)]
            plt.scatter(param, self.dic_sts['x0_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('-2.5 log10(x0)',fontsize=15)
#            plt.xlabel(self._param_names[p],fontsize=25)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/x0_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
            
            plt.subplot(3,5,p+6)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['x1_'+str(p)],yerr=self.dic_sts['x1_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, self.dic_sts['x1_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('x1',fontsize=15)
#            plt.xlabel(self._param_names[p],fontsize=25)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/x1_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
            
            plt.subplot(3,5,p+11)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
#            plt.errorbar(param, self.dic_sts['color_'+str(p)],yerr=self.dic_sts['color_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, self.dic_sts['color_'+str(p)], color='red', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('color',fontsize=15)
            plt.xlabel(self._param_names[p],fontsize=15)
#            plt.legend()
#            pdffile = '../sugar_analysis_data/results/color_'+self._param_names[p]+'.pdf'
#            plt.savefig(pdffile, bbox_inches='tight')  
#            plt.show()
#            plt.close()
        pdffile = '../sugar_analysis_data/results/mapping_sug_to_salt.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.show()
        plt.close()
        

                 
#        param = np.linspace(-4,4,self.step) 
#        plt.scatter(param, self.dic_sts['t0_1'], color='red', marker='.', zorder=1)
#        pdffile = '../sugar_analysis_data/results/t0_q2.pdf'
#        plt.savefig(pdffile, bbox_inches='tight') 
#        plt.show()
    def read_salt2(self):
        """
        """
        SCALE_FACTOR = 1e-12
        model = {}
        
        for files in [p3+'salt2_template_0.dat', p3+'salt2_template_1.dat']:
            phase, wave, values = sncosmo.read_griddata_ascii(files)
            values *= SCALE_FACTOR
            model[files] = BicubicInterpolator(phase, wave, values)
        return model
    
    def flux_salt2(self, phase, wave):
        model = self.read_salt2()
        source_salt2 = sncosmo.get_source('salt2')
        m0 = model[p3+'salt2_template_0.dat'](phase, wave)
        m1 = model[p3+'salt2_template_1.dat'](phase, wave)
        
        return self._salt_params[0] * (m0 + self._salt_params[1] * m1)*10.** (-0.4 * source_salt2._colorlaw(wave) * self._salt_params[2])

#        
#        flux_salt2 = source_salt2._flux(phase, wave)
#        return flux_salt2
                

            
    def spectrum_sugar_salt(self):
        """
        """
        dt = 1000
        
        xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
        xso1 = np.linspace(float(wl_min_salt), float(wl_min_sug), dt)
        xso2 = np.linspace(float(wl_max_sug), float(wl_max_salt), dt)

        for p in range(len(self._parameters)-1):
            if p == 3:
                list_sug_param = [-1.,0.,1.]
            else:
                list_sug_param = [-4.,0.,4.]
            for l in list_sug_param:
                self._parameters[p] = l
                res = self.fit_lc(show_lc=False)
                
                self._salt_params = [res.parameters[2],res.parameters[3],res.parameters[4]]
                spec_flux_sugar = self.model_spectrum_flux_m0(self._phase, xs)
                spec_flux_salt1 = self.flux_salt2(self._phase, xso1)
                spec_flux_salt2 = self.flux_salt2(self._phase, xso2)
                for i in range((len(self._phase)-1)/5):
                    plt.plot(xs,spec_flux_sugar[i+4,:]*(4*i+1),color='r')
                    plt.plot(xso1,spec_flux_salt1[i+4,:]*(4*i+1),color='b')
                    plt.plot(xso2,spec_flux_salt2[i+4,:]*(4*i+1),color='b')
                
                plt.text(7000,max(spec_flux_sugar[11,:])*21, self._param_names[p]+' = '+str(l), fontweight = 'bold', fontsize = 20)
                pdffile = '../sugar_analysis_data/results/spectre_sug_salt2_self._'+self._param_names[p]+'_'+str(l)+'.pdf'
                plt.savefig(pdffile, bbox_inches='tight')  
                plt.show()
            self._parameters = np.array([0., 0., 0., 0., 0.])
        
    def multi_var_par(self):
        """
        """
        nq1 = 7
        nq2 = 7
        nq3 = 7
        na = 7
        q1 = np.linspace(-6, 6, nq1)    
        q2 = np.linspace(-3, 3, nq2)
        q3 = np.linspace(-1, 1, nq3)
        av = np.linspace(-1, 1, na)
        #grey = np.linspace(30, 40, binning)

        q1, q2, q3, av = np.meshgrid(q1, q2, q3, av)

        self.m_dic_sts = {}
        x0 = np.zeros_like(q1)
        x1 = np.zeros_like(q1)
        c = np.zeros_like(q1)
        t0 = np.zeros_like(q1)
        for i in range(nq1):
            print i
            for j in range(nq2):
                    for k in range(nq3):
                        for l in range(na):
                            #for m in range(binning):
                                self._parameters = np.array([q1[i,j,k,l],           
                                                           q2[i,j,k,l], 
                                                           q3[i,j,k,l], 
                                                           av[i,j,k,l],37.])
            
                                res = self.fit_lc(show_lc=False)
                                x0[i,j,k,l] = -2.5*np.log10(res['parameters'][2])
                                x1[i,j,k,l] = res['parameters'][3]
                                c[i,j,k,l] = res['parameters'][4]        
                                t0[i,j,k,l] = res['parameters'][1]   
        self.m_dic_sts['x0'] = x0
        self.m_dic_sts['x1'] = x1
        self.m_dic_sts['c'] = c
        self.m_dic_sts['t0'] = t0
        dic_file = open('../sugar_analysis_data/results/m_sug_to_salt_test.pkl','w')
        pkl.dump(self.m_dic_sts,dic_file)
#    
    def spectrum_sugar_saltv2(self):
        """
        """
        dt = 1000
        
        xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
        xso1 = np.linspace(float(wl_min_salt), float(wl_max_salt), dt)


#        for p in range(len(self._parameters)-1):
#            if p == 3:
#                list_sug_param = [-1.,0.,1.]
#            else:
#                list_sug_param = [-4.,0.,4.]
#            for l in list_sug_param:
        self._parameters[0] = 0
#        self._parameters[p] = l
        res = self.fit_lc(show_lc=False)
        
        self._salt_params = [res.parameters[2],res.parameters[3],res.parameters[4]]
        spec_flux_sugar = self.model_spectrum_flux_m0(self._phase, xs)
        spec_flux_salt1 = self.flux_salt2(self._phase, xso1) 
        q = [2,4,8,14] #matches with phase [-6,0,12,30]
        for i in range(len(q)): 
            print i
            plt.subplot(4,1,i+1)
            plt.plot(xs,spec_flux_sugar[q[i],:],color='r')
            plt.plot(xso1,spec_flux_salt1[q[i],:],color='b')
    
    
    def mc_generator(self):
        """
        """
        data = self.AstropyTable_flux(error=True)
        
        for i in range(len(data['flux'])):
            if data['band'][i] == 'new_fI_10' or data['band'][i] =='fI_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fU_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fB_10':
                data['flux'][i] = np.random.normal(data['flux'][i],data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fV_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fR_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]

        return data
        
    def plot_multi_sds(self):
        """
        """
        # Plot Setup
        rcParams['font.size'] = 6.
        font = {'family': 'normal', 'size': 5}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=15)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)        
        
        dic_sts_s = pkl.load(open('../sugar_analysis_data/results/sug_to_salt.pkl'))    
#        dic_sts_c = pkl.load(open('../sugar_analysis_data/results/sug_to_salt_withoutibands.pkl'))    
        for p in range(4):
            if p == 4:
                 param = np.linspace(35,45,self.step)
            else:
                 param = np.linspace(-4,4,self.step)
            plt.subplot(3,5,p+1)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
    #            plt.errorbar(param, self.dic_sts['x0_'+str(p)],yerr=self.dic_sts['x0_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
    #            print self.dic_sts['x0_'+str(p)]
            plt.scatter(param, dic_sts_s['x0_'+str(p)], color='red', marker='.', zorder=1)
#            plt.scatter(param, dic_sts_c['x0_'+str(p)], color='blue', marker='.', zorder=1)
            if p == 0:
                
                plt.ylabel('$-2.5 \ log_{10}(x_{0})$',fontsize=15)
    #            plt.xlabel(self._param_names[p],fontsize=25)
    #            plt.legend()
    #            pdffile = '../sugar_analysis_data/results/x0_'+self._param_names[p]+'.pdf'
    #            plt.savefig(pdffile, bbox_inches='tight')  
    #            plt.show()
    #            plt.close()
            plt.ylim(-35.,-25)
            plt.subplot(3,5,p+6)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
    #            plt.errorbar(param, self.dic_sts['x1_'+str(p)],yerr=self.dic_sts['x1_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, dic_sts_s['x1_'+str(p)], color='red', marker='.', zorder=1)
#            plt.scatter(param, dic_sts_c['x1_'+str(p)], color='blue', marker='.', zorder=1)
            if p == 0:
                
                plt.ylabel('x1',fontsize=15)
    #            plt.xlabel(self._param_names[p],fontsize=25)
    #            plt.legend()
    #            pdffile = '../sugar_analysis_data/results/x1_'+self._param_names[p]+'.pdf'
    #            plt.savefig(pdffile, bbox_inches='tight')  
    #            plt.show()
    #            plt.close()
            plt.ylim(-2,2)
            plt.subplot(3,5,p+11)
            plt.subplots_adjust(hspace = 0.,wspace = 0.6)
    #            plt.errorbar(param, self.dic_sts['color_'+str(p)],yerr=self.dic_sts['color_err_'+str(p)], color='red', fmt='.', mfc='red', zorder=1)
            plt.scatter(param, dic_sts_s['color_'+str(p)], color='red', marker='.', zorder=1)
#            plt.scatter(param, dic_sts_c['color_'+str(p)], color='blue', marker='.', zorder=1)
            if p == 0:
                plt.ylabel('color',fontsize=15)
            plt.ylim(-2,2)
            plt.xlabel(self._param_names[p],fontsize=15)        
        pdffile = '../sugar_analysis_data/results/mapping_sug_to_salt_comp.pdf'
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.show()
        plt.close()
    
    def plot_multi_var_par(self, plot_q1=True, plot_q2=False , plot_q3=True, plot_av=True):
        """
        """
        dic_msts = pkl.load(open('../sugar_analysis_data/results/m_sug_to_salt.pkl'))
        if self.dic_sts == None:
             self.dic_sts = pkl.load(open('../sugar_analysis_data/results/sug_to_salt.pkl'))          
        param = np.linspace(-4,4,self.step)
        nq1 = 7
        nq2 = 7 
        nq3 = 7
        nav = 7
        q1 = np.linspace(-6, 6, nq1)    
        q2 = np.linspace(-3, 3, nq2)
        q3 = np.linspace(-1, 1, nq3)
        av = np.linspace(-1, 1, nav)


        rcParams['font.size'] = 6.
        font = {'family': 'normal', 'size': 5}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=15)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)  
        
        q1, q2, q3, av = np.meshgrid(q1, q2, q3, av)
        
        if plot_q1:
            for j in range(nq2):
                for k in range(nq3):
                    for l in range(nav):
                        plt.subplot(4,1,1)
                        plt.scatter(param, self.dic_sts['x0_0'], color='red', marker='.', zorder=1)
                        plt.scatter(q1[0,:,0,0], dic_msts['x0'][j,:,k,l], marker='.')     
                        plt.ylabel('-2.5 log10(x0)',fontsize=15)
                        plt.figtext(0.55, 0.8, 'q2 = '+str(q2[j,0,0,0])+' q3 = '+str(q3[0,0,k,0])+' A = '+str(av[0,0,0,l]), fontsize=12)
                        plt.subplot(4,1,2)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['x1_0'], color='red', marker='.', zorder=1)
                        plt.scatter(q1[0,:,0,0], dic_msts['x1'][j,:,k,l], marker='.')
                        plt.ylabel('x1',fontsize=15)
                        plt.subplot(4,1,3)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['color_0'], color='red', marker='.', zorder=1)
                        plt.scatter(q1[0,:,0,0], dic_msts['c'][j,:,k,l], marker='.')
                        plt.ylabel('color',fontsize=15)
                        plt.subplot(4,1,4)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['t0_0'], color='red', marker='.', zorder=1)
                        plt.scatter(q1[0,:,0,0], dic_msts['t0'][j,:,k,l], marker='.')
                        plt.ylabel('t0',fontsize=15)
                        plt.xlabel('q1',fontsize=15)
                        pdffile = '../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q1_'+str(j)+str(k)+str(l)+'.png'
                        plt.savefig(pdffile, bbox_inches='tight',dpi=200)  
                        plt.show()
                        plt.close()
        if plot_q2:
            for j in range(nq1):
                for k in range(nq3):
                    for l in range(nav):
                        plt.subplot(4,1,1)
                        plt.scatter(param, self.dic_sts['x0_1'], color='red', marker='.', zorder=1)
                        plt.scatter(q2[:,0,0,0], dic_msts['x0'][:,j,k,l], marker='.')     
                        plt.ylabel('-2.5 log10(x0)',fontsize=15)
                        plt.figtext(0.55, 0.8, 'q1 = '+str(q1[0,j,0,0])+' q3 = '+str(q3[0,0,k,0])+' A = '+str(av[0,0,0,l]), fontsize=12)
                        plt.subplot(4,1,2)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['x1_1'], color='red', marker='.', zorder=1)
                        plt.scatter(q2[:,0,0,0], dic_msts['x1'][:,j,k,l], marker='.')
                        plt.ylabel('x1',fontsize=15)
                        plt.subplot(4,1,3)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['color_1'], color='red', marker='.', zorder=1)
                        plt.scatter(q2[:,0,0,0], dic_msts['c'][:,j,k,l], marker='.')
                        plt.ylabel('color',fontsize=15)
                        
                        plt.subplot(4,1,4)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['t0_1'], color='red', marker='.', zorder=1)
                        plt.scatter(q2[:,0,0,0], dic_msts['t0'][:,j,k,l], marker='.')
                        plt.ylabel('t0',fontsize=15)
                        plt.xlabel('q2',fontsize=15)
                        pdffile = '../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q2_'+str(j)+str(k)+str(l)+'.png'
                        plt.savefig(pdffile, bbox_inches='tight',dpi=200)  
                        plt.show()
                        plt.close()

        if plot_q3:
            for j in range(nq1):
                for k in range(nq2):
                    for l in range(nav):
                        plt.subplot(4,1,1)
                        plt.scatter(param, self.dic_sts['x0_2'], color='red', marker='.', zorder=1)
                        plt.scatter(q3[0,0,:,0], dic_msts['x0'][k,j,:,l], marker='.')     
                        plt.ylabel('-2.5 log10(x0)',fontsize=15)
                        plt.figtext(0.55, 0.8, 'q1 = '+str(q1[0,j,0,0])+' q2 = '+str(q2[k,0,0,0])+' A = '+str(av[0,0,0,l]), fontsize=12)
                        plt.subplot(4,1,2)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['x1_2'], color='red', marker='.', zorder=1)
                        plt.scatter(q3[0,0,:,0], dic_msts['x1'][k,j,:,l], marker='.')
                        plt.ylabel('x1',fontsize=15)
                        plt.subplot(4,1,3)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['color_2'], color='red', marker='.', zorder=1)
                        plt.scatter(q3[0,0,:,0], dic_msts['c'][k,j,:,l], marker='.')
                        plt.ylabel('color',fontsize=15)
                        
                        plt.subplot(4,1,4)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['t0_2'], color='red', marker='.', zorder=1)
                        plt.scatter(q3[0,0,:,0], dic_msts['t0'][k,j,:,l], marker='.')
                        plt.ylabel('t0',fontsize=15)
                        
                        plt.xlabel('q3',fontsize=15)
                        pdffile = '../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q3_'+str(j)+str(k)+str(l)+'.png'
                        plt.savefig(pdffile, bbox_inches='tight',dpi=200)  
                        plt.show()
                        plt.close()
                        
        if plot_av:
            for j in range(nq1):
                for k in range(nq2):
                    for l in range(nq3):
                        plt.subplot(4,1,1)
                        plt.scatter(param, self.dic_sts['x0_3'], color='red', marker='.', zorder=1)
                        plt.scatter(av[0,0,0,:], dic_msts['x0'][k,j,l,:], marker='.')     
                        plt.ylabel('-2.5 log10(x0)',fontsize=15)
                        plt.figtext(0.55, 0.8, 'q1 = '+str(q1[0,j,0,0])+' q2 = '+str(q2[k,0,0,0])+' q3 = '+str(q3[0,0,l,0]), fontsize=12)
                        plt.subplot(4,1,2)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['x1_3'], color='red', marker='.', zorder=1)
                        plt.scatter(av[0,0,0,:], dic_msts['x1'][k,j,l,:], marker='.')
                        plt.ylabel('x1',fontsize=15)
                        plt.subplot(4,1,3)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['color_3'], color='red', marker='.', zorder=1)
                        plt.scatter(av[0,0,0,:], dic_msts['c'][k,j,l,:], marker='.')
                        plt.ylabel('color',fontsize=15)

                        plt.subplot(4,1,4)
                        plt.subplots_adjust(hspace = 0.)
                        plt.scatter(param, self.dic_sts['t0_3'], color='red', marker='.', zorder=1)
                        plt.scatter(av[0,0,0,:], dic_msts['t0'][k,j,l,:], marker='.')
                        plt.ylabel('t0',fontsize=15)
                        
                        plt.xlabel('av',fontsize=15)
                        pdffile = '../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_av_'+str(j)+str(k)+str(l)+'.png'
                        plt.savefig(pdffile, bbox_inches='tight',dpi=200)  
                        plt.show()
                        plt.close()
                        
    def gif_multi_var_par(self):
        """
        """
        #gif q1
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q1_*33.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq1_q2.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q1_3*3.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq1_q3.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q1_33*.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq1_a.gif')
        
        #gif q2
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q2_*33.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq2_q1.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q2_3*3.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq2_q3.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q2_33*.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq2_a.gif')
        
        #gif q3
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q3_*33.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq3_q1.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q3_3*3.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq3_q2.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_q3_33*.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationq3_a.gif')
                
        #gif av
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_av_*33.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationav_q1.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_av_3*3.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationav_q2.gif')
        os.system('convert -delay 80 ../sugar_analysis_data/results/multi_var_sugar_to_salt/multi_var_mapping_sug_to_salt_comp_av_33*.png ../sugar_analysis_data/results/multi_var_sugar_to_salt/gif/animationav_q3.gif')
        
    def plot_example_filter(self):
        
        band_used = ['fB_10','fV_10','fR_10']     
        for b in band_used:
            
            filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+b+'.dat')
            wlen = filt2[:,0]
            tran = filt2[:,1]
            self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
            
            #computation of the integral
            dt = 100000
            xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
            self.xs = xs
            spec_flux = self.model_spectrum_flux_m0(self._phase, xs)

            plt.plot(self.xs,self.splB(self.xs))
             
        plt.plot(self.xs,spec_flux[3,:])



        
        
                