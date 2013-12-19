# -*- coding: utf-8 -*-
"""
Created on Oct  29 2013

@author: mikkel
"""
from __future__ import division
from matplotlib import pyplot as plt
plt.ion()
import numpy as np
from scipy import stats

def plot_ana_spk(anasigs, spktrains):
    #print first trace plus spiketrain
    #seg = block.segments[0]
    nb_figs = len(anasigs) + 1
    fig = plt.figure()
    
    ax2 = fig.add_subplot(nb_figs,1,nb_figs)
    ax2.set_title('spikes')
    for s,st in enumerate(spktrains):
        ax2.vlines(st, s + .5, s + 1.5, color='k')
        #ax2.plot(st, s*np.ones(st.size), linestyle = 'None', 
        #            marker = ',', color = 'k')
        #ax = plt.gca()
        #for ith, trial in enumerate(event_times_list):
        #plt.vlines(trial, ith + .5, ith + 1.5, color=color)
        #plt.ylim(.5, len(event_times_list) + .5)
        #return ax
    plt.ylim(.5, len(spktrains) + .5)
    for nb, asig in enumerate(anasigs):
        tmp_ax = fig.add_subplot(nb_figs,1,nb+1,sharex = ax2 )
        tmp_ax.set_title('axon: ' + str(nb+1))
        tmp_ax.plot(asig.times, asig)
    return fig
        
    #plt.show()

def spk_trig_avg_detekt(anasigs, spktrains):
    """ check if the membrane potential changes when units fire.
    For each unit: make a spike triggered average. Make a test to check if the
    distribution of membrane potentials changes after spikes"""
    
    for nb_sig, anasig in enumerate(anasigs):
        #get length of window
        #TODO: Dette er en meget grim måde at gøre det på...
        #TODO: Brug quantities class...    
        tmp = anasig[(anasig.times>(10-0.005))&(anasig.times<(10+0.005))]
        length = len(tmp)
        #xs = tmp.times
    
        res2 = np.zeros([len(spktrains),2])
        for s,st in enumerate(spktrains):
            spksig = st
            res = np.zeros([len(spksig), length])
            for nb, spk in enumerate(spksig):
                #include only when more than 10 spikes
                if len(spksig)>10:
                    t = spk.magnitude
                    tmp =  anasig[(anasig.times>(t-0.005))&(anasig.times<(t+0.005))]
                    if len(tmp)>=length:
                        res[nb,:] = tmp[0:length]-np.mean(tmp[0:length])
            v = np.var(np.mean(res,0))*(res.shape[0])
            res2[s,:] = [s,v]
            print str(s) + '/' + str(len(spktrains))
        #stats.ttest_ind(rvs1,rvs2, equal_var = False)
        plt.figure()
        plt.plot(res2[:,1])
    

def plot_spk_trig_avg(data):
    #calculate lengths
    anasig0 = data[0][0]
    spktrains0 = data[0][1]
    tmp = anasig0[(anasig0.times>(10-0.005))&(anasig0.times<(10+0.005))]
    length = len(tmp)
    len_spktrains = len(spktrains0)
    for start_nb in range(0,len_spktrains,25):
        slut_nb = min((start_nb+25),len_spktrains)
        fig = plt.figure()
        for plot_i, i in enumerate(range(start_nb,slut_nb)):  
            res = np.zeros([0, length])
            for anasig, spktrains in data:
                spksig = spktrains[i]
                if (len(spksig)>0):
                    for nb, spk in enumerate(spksig):
                        t = spk.magnitude
                        tmp =  anasig[(anasig.times>(t-0.005))&(anasig.times<(t+0.005))]
                        if len(tmp)>=length:
                            res = np.append(res,[(tmp[0:length]-np.mean(tmp[0:length]))],axis=0)   
            ax = fig.add_subplot(5,5,plot_i+1)
            ax.set_title(spksig.name + '(' + str(i) + ')')
            if (res.shape[0]>0):
                ax.plot(res.transpose(), color='blue', linewidth=0.5)
                ax.plot(np.mean(res,0), color='green', linewidth=2)
        plt.show()
        print slut_nb
        
def cross_corr(c, sr):
    #make time values
    length = c.shape[0]
    t_length = length/sr
    x = np.linspace(-t_length/2, t_length/2, length)
    plt.figure()
    plt.plot(x, c)
    