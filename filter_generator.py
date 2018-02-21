# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:41:13 2018

@author: mondon
"""

import numpy as N 
from scipy import stats
import os , sys 
r_name = '../sugar_analysis_data/using_saltmine/data/Instruments/Florian/'



#fU = stats.norm.pdf(wl, loc=3701, scale= 350)
#fU /= fU.max()
#fB = stats.norm.pdf(wl, loc=4600, scale= 350)
#fB /= fB.max()
#fV = stats.norm.pdf(wl, loc=5700, scale= 350)
#fV /= fV.max()
#fR = stats.norm.pdf(wl, loc=6850, scale= 350)
#fR /= fR.max()

def filter_trans(wl, wlred, wlblue,width): 
    f = N.zeros(len(wl))
    for i in range(len(wl)):
        f[i] = (1+N.tanh((wl[i]-wlblue)/width))*(1+N.tanh(-(wl[i]-wlred)/width))/4
    return f
    
def filter_used(width, filters=['USNf','BSNf','VSNf','RSNf','ISNf']):
    """
    """
    
    wl = N.arange(2000,10000,0.1)

    for filter_name in filters:

         if filter_name == 'USNf':
            wlblue = 3300
            wlred = 4102
            fU = filter_trans(wl, wlred, wlblue,width)
            fileU = open(r_name+'fU_'+str(width)+'.dat','w')
            for j in range(len(wl)):
                fileU.write("%.3f %.3f \n"%(wl[j],fU[j]))
            fileU.close()
         elif filter_name == 'BSNf':
            wlblue = 4102
            wlred = 5100
            fB = filter_trans(wl, wlred, wlblue,width)
            fileB = open(r_name+'fB_'+str(width)+'.dat','w')
            for j in range(len(wl)):
                fileB.write("%.3f %.3f \n"%(wl[j],fB[j]))
            fileB.close()
         elif filter_name == 'VSNf':
            wlblue = 5200
            wlred = 6289
            fV = filter_trans(wl, wlred, wlblue,width)        
            fileV = open(r_name+'fV_'+str(width)+'.dat','w')
            for j in range(len(wl)):
                fileV.write("%.3f %.3f \n"%(wl[j],fV[j]))
            fileV.close()
         elif filter_name == 'RSNf':
            wlblue = 6289
            wlred = 7607
            fR = filter_trans(wl, wlred, wlblue,width)  
            fileR = open(r_name+'fR_'+str(width)+'.dat','w')
            for j in range(len(wl)):
                fileR.write("%.3f %.3f \n"%(wl[j],fR[j]))
            fileR.close()
         elif filter_name == 'ISNf':
            wlblue = 7607
            wlred = 9200
            fI = filter_trans(wl, wlred, wlblue,width)  
            fileI = open(r_name+'fI_'+str(width)+'.dat','w')
            for j in range(len(wl)):
                fileI.write("%.3f %.3f \n"%(wl[j],fI[j]))   
            fileI.close()
         else:
            raise ValueError('Error this band not implemented in filter generator')


########
##Main##
########
if __name__=='__main__':
   
   filter_used(500)
   filter_used(100)
   filter_used(40)
   filter_used(20)
   filter_used(10)
   filter_used(5)
   filter_used(5.5)
   filter_used(1.5)
   filter_used(1)
   filter_used(2)
   filter_used(2.5)
   filter_used(3)
   filter_used(0.1)
   