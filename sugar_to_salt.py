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

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27       # [erg s]
# wavelength limits for sugar model
wl_min_sug = 3254.01639
wl_max_sug = 8649.03871
wl_min_salt = 2000.000000
wl_max_salt = 9200.000000
from astropy.table import Table
try:
    Build_SNF.register_SNf_bands_width(width=10)
    Build_SNF.mag_sys_SNF_width(width=10) 
except:
    print 'Filters and mag sys already registred'
class sugar_spectrum():
    
    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
            		m2file='sugar_template_2.dat',
            		m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat', 
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
        self._parameters = np.array([0., 0., 0., 0., 0.])
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
    
    def AstropyTable_flux(self):
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
                dxs = (float(wl_max_sug-wl_min_sug)/(dt-1))
                spec_flux = self.model_spectrum_flux_m0(self._phase, xs)
#                plt.plot(xs,self.splB(xs))
#                plt.plot(xs,spec_flux[i,:])
#                print spec_flux[i,:] 
                plt.show()
                inte = np.sum(spec_flux[i,:]*(xs  / (CLIGHT*HPLANCK))*self.splB(xs)*dxs)
                    
                flux.append(inte)
                fluxerr.append(0.)
                time.append(self._phase[i])
                band.append(b)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
        data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
        return data
    
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
        res, fitted_model = sncosmo.fit_lc(data, model, ['x0','x1','c'], modelcov = True)
        if show_lc :
            sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
        return res
#                print self.flux
    
    def sugar_dep_salt(self):
        '''
        '''
        x0 = []
        x0_err = []
        x1 = [] 
        x1_err = []
        color = [] 
        color_err = []
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
                 x1.append(res['parameters'][3])
                 color.append(res['parameters'][4])
                 x0_err.append(res['errors']['x0']*-2.5*np.log10(res['parameters'][2])/ 1.0857362047581294)
                 x1_err.append(res['errors']['x1'])
                 color_err.append(res['errors']['c'])
            self.dic_sts['x0_'+str(p)] = x0
            self.dic_sts['x1_'+str(p)] = x1
            self.dic_sts['color_'+str(p)] = color
            self.dic_sts['x0_err_'+str(p)] = x0_err
            self.dic_sts['x1_err_'+str(p)] = x1_err
            self.dic_sts['color_err_'+str(p)] = color_err
            x0 = []
            x0_err = []
            x1 = [] 
            x1_err = []
            color = [] 
            color_err = []   
            self._parameters = np.array([0., 0., 0., 0., 0.])
        dic_file = open('../sugar_analysis_data/results/sug_to_salt.pkl','w')
        pkl.dump(self.dic_sts,dic_file)
        
    def plot_sds(self):
        
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
        plt.close()
    def read_salt2(self):
        '''
        '''
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
        '''
        '''
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
        
