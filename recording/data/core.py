from neo import io
import numpy as np
from klustakwikio2 import KlustaKwikIO as KlustaIO
from neo.core import *
import os
from neo.io.hdf5io import NeoHdf5IO
import quantities as pq
import matplotlib.pyplot as plt
import math


class Data(object):
    """
    Responsible for I/O functions.
    2013 Mikkel Vestergaard
    """

    def __init__(self, rec):
        """
        Initialize the NEO Block.
        """
        self.rec = rec
        self.rec.block = Block()
        max_seg = max(self.rec.par['labbook'].keys())
        for idx in xrange(max_seg+1):
            if idx in self.rec.par['labbook']:
                seg = Segment(name=(idx), index=int(idx))
                self.rec.block.segments.append(seg)
            else:
                seg = Segment(name='empty')
                self.rec.block.segments.append(seg)
        
    def load_axon(self):
        """
        Load data recorded with Axon Digidata into the NEO
        """
        for data_nb in self.rec.par['labbook']:
            axon_nb = self.rec.par['labbook'][data_nb][0]
            if axon_nb == '':
                break
            fname = os.path.join(self.rec.par['root_raw'],self.rec.par['folder_axon'],
                    self.rec.par['base_axon'] + axon_nb + '.abf')
            r = io.AxonIO(filename=fname)
            bl = r.read_block(lazy=False, cascade=True)
            seg = bl.segments[0]
            #Find and adjust the offset so it fits with amplipex
            idx = np.where(seg.analogsignals[self.rec.par['axon_sync'][data_nb]]>3)[0][0]
            offset = seg.analogsignals[self.rec.par['axon_sync'][data_nb]].times[idx]
            for anasig in seg.analogsignals:
                anasig.t_start = -(offset-anasig.t_start)
            self.rec.block.segments[data_nb].merge(seg)
	    print 'merged ' + str(data_nb) + ' with axon_nb: ' + str(axon_nb)    

    def load_amplipex(self):
        """
        Load data recorded with amplipex into the NEO
        """
        for data_nb in self.rec.par['labbook']:
            ampli_nb = self.rec.par['labbook'][data_nb][1]
            if ampli_nb == '':
                break
            fname = os.path.join(self.rec.par['root_raw'],self.rec.par['folder_ampli'],
                    self.rec.par['base_ampli'] + ampli_nb + '.dat')
            r = io.AmpliIO(filename=fname)
            bl = r.read_block(lazy=False, cascade=True)
            seg = bl.segments[0]
            self.rec.block.segments[data_nb].merge(seg)
            print 'merged ' + str(data_nb) + ' with ampli_nb: ' + str(ampli_nb)
                  
    def load_klusta(self):
        """
        Load units from files saved in klustakwik format into the NEO
        """
        r2 = KlustaIO(os.path.join(self.rec.par['root'], 'sortering/sortering'), sampling_rate=39303) #TODO: read pars from somewhere
        blck = r2.read_block(lazy = False, cascade = True, waveform = False)
        for seg in blck.segments:
            for data_nb, axon_nb, ampli_nb in self.rec.par['labbook']:
                if seg.name == data_nb:
                    self.rec.block.segments[data_nb].merge(seg)
                    print 'merged ' + data_nb + ' with ampli_nb: ' + str(ampli_nb)
                    break

    def save_hdf5(self):
        """
        Save NEO into hdf structure. Not used because it is fast to load directly from .dat
        """
        iom = NeoHdf5IO(os.path.join(self.rec.par['root'],'neo.h5'))
        iom.save(self.block)

    def load_hdf5(self):
        """
        Load a NEO object from a hdf structure
        """
        pass
