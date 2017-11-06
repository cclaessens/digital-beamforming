# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 13:49:04 2017

@author: claessen
"""
import numpy as np
import math as m

class Signal:
    # a signal object has an amplitude, a frequency and a position in r-phi
    def __init__(self, amp, f, r, phi, delay):
        self.amplitude = amp
        self.frequency = f
        self.r = r
        self.phi = phi
        self.position = [self.r, self.phi]
        self.delay = delay
        self.phase = 2*m.pi*self.frequency*self.delay

    def TimeDomainSignal(self, t):
        # returns the time domain signal
        return self.amplitude*np.exp(1j*(2.0*m.pi*self.frequency*t-self.phase))
                                     
    def FreqDomainSignal(self, t):
        # returns the frequency domain signal
        return [np.fft.fftshift(np.fft.fftfreq(np.shape(t)[0])*1/(t[1]-t[0])), np.fft.fftshift(np.fft.fft(self.TimeDomainSignal(t))/np.shape(t)[0])]


    def TimeDomainMixDown(self, t):
        return self.amplitude*np.exp(1j*(2.0*m.pi*self.mixfreq*t - self.phase))
                                     
    def FreqDomainMixDown(self, t):
        return [np.fft.fftshift(np.fft.fftfreq(np.shape(t)[0])*1/(t[1]-t[0])), np.fft.fftshift(np.fft.fft(self.TimeDomainMixDown(t))/np.shape(t)[0])]

    def MixDownSignal(self, F_mix, t):
        # call this signal to mix down with a frequency F_mix
        # returns the mixed down signal in time and frequency domain
        #print('Mix down with: ', F_mix)
        self.mixfreq = self.frequency - F_mix

        #print('frequency is now: ', self.mixfreq*10**-6, 'MHz')
        return self.TimeDomainMixDown(t), self.FreqDomainMixDown(t)

