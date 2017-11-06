# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 13:46:32 2017

@author: claessen
"""

import numpy as np
import math as m

#import Signal
from Signal import Signal



class Antenna:
    def __init__(self, r, phi):
        # so far every antenna only has a position in r-phi
        self.r = r
        self.phi = phi
        self.gain = 1.0

    def DistanceToSource(self, signal):
        return np.sqrt(self.r**2+signal.r**2-2*self.r*signal.r*np.cos(self.phi-signal.phi))

    def TimeDelay(self, signal):
        c0 = 3.0*10**8
        d = self.DistanceToSource(signal)
        return d/c0

   # def ReceivedSignalAfterMixDown(self, signal, t):
   #     # Assuming dipole radiation propagation ~1/r
   #     #print('propagation 1/r**2')
   #     self.signal = signal.TimeShiftedSignal(self.TimeDelay(signal), t)/(self.DistanceToSource(signal))
   #     return self.signal
    
    def ReceiveSignalBeforeMixDown(self, signal, t):
        # Assuming dipole radiation propagation ~1/r

        self.received_signal = Signal(signal.amplitude*1/self.DistanceToSource(signal)*self.gain, signal.frequency, 0, 0, delay=self.TimeDelay(signal))
        
        return self.received_signal.TimeDomainSignal(t)
    
    def MixdownReceivedSignal(self, fMix, t):
        return self.received_signal.MixDownSignal(fMix, t)[0]
    
    def MixAndDigitizeSignal(self, fMix, tmax, sr):
        t_sampled = np.arange(0, tmax,1/sr)
        self.digitized_signal = self.MixdownReceivedSignal(fMix, t_sampled)
        return t_sampled, self.digitized_signal


class AntennaArray:
    # a collection of antennas in a ring, initialized with radius and numbers os antennas
    def __init__(self, R, N):
        self.NAntennas = N
        self.Antennas, self.positions, self.positions_xy = [], [], []

        for i in range(N):
            self.Antennas.append(Antenna(R, 2.0*m.pi/N*i))
            self.positions.append([R, 2.0*m.pi/N*i])
            self.positions_xy.append([R*np.sin(2.0*m.pi/N*i), R*np.cos(2.0*m.pi/N*i)])

        self.positions = np.array(self.positions)
        self.positions_xy = np.array(self.positions_xy)

    def ReceiveSignals(self, signals, t, fMix, sr, tmax):
    # Calculate signals as seen by antennas and pass them to the beamformers so not everyboday recalculates
        #print('Calculating signal as seen by antennas')
        NSignalBins = np.shape(t)[0]
        #np.shape(self.Antennas[0].ReceivedSignalAfterMixDown(signals[0], t))[0]
        #print("NSignalBins", NSignalBins)

        t_sampled = np.arange(0, tmax,1/sr)
        ReceivedSignals = []
        
        # signals can be a list of electron signals
        for i in range(self.NAntennas):
            s = 0.0*t_sampled*(1+1j)
            for j in range(len(signals)):
                self.Antennas[i].ReceiveSignalBeforeMixDown(signals[j], t)
                t_sampled, sj = self.Antennas[i].MixAndDigitizeSignal(fMix, tmax, sr)
                s += sj
            ReceivedSignals.append(s)
        #print('Done')
        return ReceivedSignals