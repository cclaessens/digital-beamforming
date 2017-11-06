# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:05:49 2017

@author: chrischtel
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
from Beamformer import Beamformer

class SetOfBeamformers:
    def __init__(self, NRings, NAngles1, Array, fMix, SR, tmax):
        # antenna array
        self.Array = Array
                
        # dictionaries
        self.Beamformers = {}
        self.FormedSignal = {}
        self.FormedSpectrum = {}
        self.MaxPeak = {}
        self.MaxPeakFreq = {}        
        
        #other
        self.NAntennas = len(self.Array.Antennas)
        self.Nbeamformers = 0
        self.threshold = 10**5
        self.fMix = fMix
        self.SR = SR
        self.tmax = tmax
        
        print('Initializing beamformers...')
        start = time.time()
        
        # add a beamformer to the beamformer dictionary for every combination of r and phi
        self.targets = []
        self.Beamformers[str(self.Nbeamformers)]=Beamformer(0.0, 0.0, self.Array, self.fMix, self.SR, self.tmax)
        self.targets.append([0.0, 0.0])
        
        r = np.linspace(0, 4.0*(1-1/NRings), NRings)/100
       
        
        for i in range(NRings-1):
            
            nphi = (i+1)*NAngles1
            phi = np.linspace(0,2.0*np.pi*(1-1/nphi), nphi)
            
            for j in range(nphi):
                self.Nbeamformers +=1
                self.Beamformers[str(self.Nbeamformers)]=Beamformer(r[i+1], phi[j], self.Array, self.fMix, self.SR, self.tmax)
                
                self.targets.append([r[i+1], phi[j]])
        end = time.time()
        print('... took {} s'.format(np.round(end-start,3)))
        print('Got {} beamformers.'.format(self.Nbeamformers+1))
        
    def Plotter(self):
        # plot the array and target positions
        fig=plt.figure(figsize=(16, 7), dpi= 80, facecolor='w', edgecolor='k')
        plt.subplot(121)
        plt.scatter(self.Array.positions_xy.T[0], self.Array.positions_xy.T[1], s=50, c = self.Array.positions.T[1])
        for i in self.targets:
            plt.scatter(i[0]*np.cos(i[1]), i[0]*np.sin(i[1]), s=30, color='grey')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Array and beamformer targets')
                
        
    def DistributeSignal(self, signals):
        # call PhaseShift method of every beamformer and collect returns
        print('Distributing signal among beamformers...')
        start = time.time()   
        
        for i in self.Beamformers.keys():
            #print('Beamformer', i, ':')
            ShiftedSignals, self.FormedSignal[i], self.FormedSpectrum[i] = self.Beamformers[i].PhaseShift(signals)
        end = time.time()
        print('... took {} s.'.format(np.round(end-start,3)))
        
        
        
    def CollectAbsoluteMaxima(self): 
        # Find maximal signal and highlight it in plot
        for i in self.Beamformers.keys():
            self.MaxPeak[i] = np.max(np.abs(self.FormedSpectrum[i][1]))
            self.MaxPeakFreq[i] = self.FormedSpectrum[i][0][np.abs(self.FormedSpectrum[i][1])==np.max(np.abs(self.FormedSpectrum[i][1]))]
        
        
        highlight = max(self.MaxPeak, key=lambda i: self.MaxPeak[i])
        
        self.Plotter()
        r = self.Beamformers[highlight].r
        phi = self.Beamformers[highlight].phi
        print('Beamfomer {} looking at [r, phi/spi] = {}, {} has strongest signal.'.format(highlight, r, round(phi/(2*m.pi), 3)))
        print('The frequency of that signal is {} MHz.'.format(round(self.MaxPeakFreq[highlight][0]*10**-6)))
        
        plt.scatter(r*np.cos(phi), r*np.sin(phi), s=50, color='red')
        plt.subplot(122)
        plt.bar(range(len(self.MaxPeak)), self.MaxPeak.values(), 0.1, align='center')
        plt.xlabel('Beamformer number')
        plt.ylabel('Observed signal amplitude')
        #plt.xticks(range(len(self.MaxPeak)), self.MaxPeak.keys())
        plt.subplot(121)
        
    
    def FormerOverFreq(self, NFreqBins=10):
    # highlight all beamformer that found a signal above threshold for frequency ranges
    # map
    
        self.Plotter()
        
        freq = np.linspace(0, 10**9, NFreqBins)
        d = self.Beamformers
        d = {float(k) for k in d.keys()}    
        F, B = np.meshgrid(freq, np.array(list(d)))
        
        Z = np.zeros(np.shape(F))
        
        for f in range(NFreqBins):
            freqlimits = [f*10**9/NFreqBins, (f+1)*10**9/NFreqBins]
            iBins = np.where((self.FormedSpectrum['0'][0]>freqlimits[0]) & (self.FormedSpectrum['0'][0]<=freqlimits[1]))
            #print(freqlimits, np.shape(iBins), np.shape(self.FormedSpectrum['0'][0]))
            
            for i in self.Beamformers.keys():   
                #for multiple detections
                self.MaxPeak[i] = np.max(np.abs(self.FormedSpectrum[i][1][iBins]))
                self.MaxPeakFreq[i] = self.FormedSpectrum[i][0][np.abs(self.FormedSpectrum[i][1])==np.max(np.abs(self.FormedSpectrum[i][1][iBins]))]
                
                # for surf plot
                Z[int(i),f] = np.max(np.abs(self.FormedSpectrum[i][1][iBins]))
                
            threshold = self.threshold
            highlight = max(self.MaxPeak, key=lambda i: self.MaxPeak[i])
            
            if self.MaxPeak[highlight]>threshold:
                r = self.Beamformers[highlight].r
                phi = self.Beamformers[highlight].phi
                print('Beamfomer', highlight, 'looking at [r, phi/2pi]=', r, round(phi/(2*m.pi), 3), 'has strongest signal')
                print('The frequency of that signal is', round(self.MaxPeakFreq[highlight][0]*10**-6), 'MHz')
                plt.scatter(r*np.cos(phi), r*np.sin(phi), s=50, color='red')
                
                
        plt.subplot(122)  
        plt.pcolormesh(F*10**-6, B, (Z/np.max(Z)), cmap='jet')
        plt.colorbar()
        plt.xlabel('F [MHz]')
        plt.ylabel('Beamformer no')