# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:10:15 2017

@author: claessen
"""
import numpy as np
import math as m
import matplotlib.pyplot as plt


class Beamformer:
    def __init__(self, r, phi, array, fMix, SR, tmax):
        # r and phi is position of where to look at
        self.r = r
        self.phi = phi
        # beamformer gets array to access positions and antenna methods
        self.array = array
        
        self.NBins = int(tmax*SR)
        self.freqs = (np.fft.fftfreq(self.NBins))*SR
        self.fMix = fMix

        # for for loops
        
        self.SR = SR
        self.NAntennas = len(self.array.Antennas)

        # calculate weights during initialization. never have to be recalculated
        self.Weights()

    def TimeDelays(self):
        # time delay and distances from point to look at to antennas in array
        c0=3.0*10**8
        return self.Distances()/c0

    def Distances(self):
        r2 = self.array.positions.T[0]
        phi2 = self.array.positions.T[1]
        return np.sqrt(self.r**2+r2**2-2*self.r*r2*np.cos(self.phi-phi2))

    def Weights(self):
        # because the beamformer gets mixed down signals
        # it needs to know the mix down frequency to calculate the correct weights
        self.weights = np.zeros((self.NAntennas, self.NBins))
        self.weights = np.exp(1j*2.0*np.pi*(np.outer(self.TimeDelays(), self.freqs+self.fMix)))
        #plt.plot(self.freqs*10**-6, self.weights[0])
        #plt.xlabel('Frequency [MHz]')
        #plt.ylabel('real(weights)')

    def CalculateReceivedSignals(self, signals):
        raise("not implemented")


    def PhaseShift(self, ReceivedSignals=None, signals=None):
        # This method does the beamforming job
        # give it a list of signal objects and it will calculate what the antennas see (if not passed as argument)
        # and apply the beamforming ignoring its knowledge about the signal
        #print('Phase shifting signals')
        if ReceivedSignals!=None:
            self.ReceivedSignals = ReceivedSignals
        else:
            self.ReceivedSignals = self.CalculateReceivedSignals(signals)

        self.wsignal = np.zeros(np.shape(self.weights))*(1+1j)
        self.wspectrum = np.zeros(np.shape(self.weights))*(1+1j)

        for i in range(self.NAntennas):
            s = ReceivedSignals[i]
            #print(i, s0[0])
            self.fspectrum = np.fft.fft(s)
            self.wspectrum[i] = self.fspectrum*self.weights[i]
            self.wsignal[i]=np.fft.ifft(self.wspectrum[i])

        # sum over the shifted signals and return list with individual shifted signals and sums
        self.formedspectrum = np.sum(self.wspectrum, axis=0)/np.shape(s)[0]
        self.formedsignal = np.sum(self.wsignal, axis=0)

        return self.wsignal, self.formedsignal, [self.freqs, self.formedspectrum]

    def Plotter(self):
        t = np.linspace(0, 1/self.SR*(self.NBins-1), self.NBins)
        # plot seen signals, shifted signal, sum and spectrum
        fig=plt.figure(figsize=(16, 7), dpi= 80, facecolor='w', edgecolor='k')
        plt.subplot(121)


        for i in range(self.NAntennas):
            plt.plot(t*10**9, np.real(self.ReceivedSignals[i]))
            plt.xlim([0, 5])
            plt.xlabel('ns')
            plt.title('Seen by antennas')

        plt.subplot(122)
        for i in range(self.NAntennas):
            plt.plot(t*10**9, np.real(self.wsignal[i]))
            plt.xlim([0, 5])
            plt.xlabel('ns')
            plt.title('Phase shifted signals')


        fig=plt.figure(figsize=(16, 7), dpi= 80, facecolor='w', edgecolor='k')
        plt.subplot(121)

        plt.plot(t*10**9, np.real(self.formedsignal))
        plt.xlim([0, 5])
        plt.xlabel('ns')
        plt.title('Formed signal')

        plt.subplot(122)
        plt.plot(np.fft.fftshift(self.freqs*10**-6), (np.abs(np.fft.fftshift(self.formedspectrum))))
        plt.xlim([0, 1000])
        plt.xlabel('MHz')
        plt.title('Formed spectrum')