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
        #Initialize neo sections based on labbook:
        #table for data... maybe based on XML...
        #self.labbook = [[0,'0009','001'],
        #    [1,'0011','002']] #[data_nb, axon_db, ampli_nb]
        for data_nb, axon_nb, ampli_nb in self.rec.par['labbook']:
            seg = Segment(name=(data_nb), index= int(data_nb))
            self.rec.block.segments.append(seg)
        
    def load_axon(self):
        """
        Load data recorded with Axon Digidata into the NEO
        """
        for data_nb, axon_nb, ampli_nb in self.rec.par['labbook']:
            fname = os.path.join(self.rec.par['root'],'axon',self.rec.par['base_axon'] + '_' + axon_nb + '.abf')
            r = io.AxonIO(filename=fname)
            bl = r.read_block(lazy=False, cascade=True)
            seg = bl.segments[0]
            #Find and adjust the offset so it fits with amplipex
            idx = np.where(seg.analogsignals[9]>3)[0][0]
            offset = seg.analogsignals[9].times[idx]
            for anasig in seg.analogsignals:
                anasig.t_start = -(offset-anasig.t_start-0.2*pq.s) #TODO: De 200 ms kommer fra klippemetoden. Skal laves smartere... 
            self.rec.block.segments[data_nb].merge(seg)
            print 'merged ' + data_nb + ' with axon_nb: ' + str(axon_nb)
            break

    def load_amplipex(self):
        """
        Load data recorded with amplipex into the NEO
        """
        for data_nb, axon_nb, ampli_nb in self.rec.par['labbook']:
            fname = os.path.join(self.rec.par['root'],'axon',self.rec.par['base_axon'] + '_' + axon_nb + '.abf')
            r = io.AmpliIO(filename=fname)
            bl = r.read_block(lazy=False, cascade=True)
            seg = bl.segments[0]
            self.rec.block.segments[data_nb].merge(seg)
            print 'merged ' + data_nb + ' with ampli_nb: ' + str(ampli_nb)
            break
        data = self.read_ch_ampli([self.nb_channels-1], dataset)

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

    def load_data_ampli(self):
        """
        Split all .dat files
        """
        data_dir = os.path.join(self.root,'data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
            folderPath = os.path.join(self.root, self.folder_amplipex)
            print 'Loading .dat files from: ' + folderPath
            for root, subFolders, files in os.walk(folderPath):
                del subFolders
                for filename in files:
                    #print filename
                    name, extension = os.path.splitext(filename)
                    file_nb = name.split('_')[-1]
                    if (extension==".dat"):
                        self._split_ampli_dat(filename,file_nb)
        else:
            print 'Data folder already exists'
    
    def _split_ampli_dat(self, filename, file_nb):
        """
        Splits .dat files into one file for each channel
        """
        try:
            file_dir = os.path.join(self.root,'data',file_nb)
            if not os.path.exists(file_dir):
                filepath = os.path.join(self.root, self.folder_amplipex, filename)
                #open and split
                fid_orig = open(filepath, "rb")
                os.mkdir(file_dir)
                fid_new = []
                for i in range(0,self.nb_channels):
                    fid_new.append(open(os.path.join(self.root,'data',file_nb, str(i) + '.dat'), 'wb'))
                for data in self._read_chunk_ampli_dat(fid_orig):
                    data = np.reshape(data,(-1,self.nb_channels))
                    for i in range(0,self.nb_channels):
                        #fid_new[i].write(data[:,i].tolist())
                        data[:,i].tofile(fid_new[i])
                print 'successfully split file: ' + filename
            else:
                print 'folder already created...'
                print 'failed split file: ' + filename
        except IOError, e:
            print e
            print 'Couldn\'t split file: ' + filename
    
    def read_ch_ampli(self, channels, dataset):
        """
        Load the specified channels from a dataset
        """
        #find the .dat file to load
        for files in os.listdir(os.path.join(self.root, self.folder_amplipex)):
            if files.endswith(dataset + ".dat"):
                filename = files        
        filepath = os.path.join(self.root, self.folder_amplipex, filename)
        fid = open(filepath, "rb")
        length = os.fstat(fid.fileno()).st_size/(2*self.nb_channels)
        #Read into the 'data' var
        data = np.zeros((length,len(channels)),dtype='int16')
        idx = 0
        for tmp_data in self._read_chunk_ampli_dat(fid):
            tmp_data = np.reshape(tmp_data,(-1,self.nb_channels))
            data[idx:(idx+tmp_data.shape[0]),:] = tmp_data[:,channels]
            idx += tmp_data.shape[0]
        return data
            
    def _read_chunk_ampli_dat(self, file_object):
        """
        Generator for reading the ampli dat file in chunks
        """
        while True:
            data = np.fromfile(file_object,'int16',1000000*self.nb_channels)       
            if not data.shape[0]>0:
                break
            yield data
            
    def read_ch_dat(self, channel, dataset):
        """obsolete??"""
        filepath = os.path.join(self.root, 'data', dataset, str(channel) + '.dat')
        fid = open(filepath, 'rb')
        data = np.fromfile(fid,'int16')
        return data
        
    def save_files(self):
        """obsolete??"""
        pass
    
    def save_combined_dat(self, channels, dataset):
        """obsolete??"""
        data = []
        for i in range(0,len(channels)):
            ch_data = self.read_ch_dat(channels[i], dataset)
            ch_data = np.reshape(ch_data,(-1,1))
            if i == 0:
                data = ch_data
            else:
                data = np.concatenate((data,ch_data),1)
        filepath = os.path.join(self.root,'data',dataset, 'combi.dat')
        fid = open(filepath,'wb')
        data.tofile(fid)
        
    def cut_files(self):
        """obsolete??"""
        pass

    def write_vars_to_file(self, fid, variables):
        for (name, val) in variables.items():
            fid.write("%s = %s\n" % (name, repr(val)))
            
    def read_raw_dat_names(self):
        raw_data = []
        folderPath = os.path.join(self.root, self.folder_amplipex)
        for root, subFolders, files in os.walk(folderPath):
            files.sort()
            del subFolders
            for filename in files:
                name, extension = os.path.splitext(filename)
                if (extension==".dat"):
                    raw_data.append(filename)
        return raw_data
    
    def cut_ampli(self):
        """
        Trims all data recorded with amplipex.
        """ 
        folderPath = os.path.join(self.root, self.folder_amplipex)
        print 'Cutting .dat files from: ' + folderPath
        for root, subFolders, files in os.walk(folderPath):
            del subFolders
            for filename in files:
                #print filename
                name, extension = os.path.splitext(filename)
                file_nb = name.split('_')[-1]
                if (extension==".dat"):
                    self._cut_ampli_set(file_nb)
                        
    def _cut_ampli_set(self, dataset):
        """
        Trims one dataset: Cut 200ms before and 200ms after the trigger
        """
        #identify datapoint before and after:
        data = self.read_ch_ampli([self.nb_channels-1], dataset)
        idx_start = np.where(data>5000)[0]
        if len(idx_start)>0:
            start = idx_start[0]
        else:
            start = 0
        idx_end = np.where(data<-5000)[0]
        if len(idx_end)>0:
            end = idx_end[0]
        else:
            end = data.shape[0]
        #find where to cut. leave 200 ms before and after trial
        delta_t = int(math.floor(0.2*self.sr))
        if start>delta_t:
            cut_start = start-delta_t
        else:
            cut_start = 0
        if end<(data.shape[0]-delta_t):
            cut_end = end+delta_t
        else:
            cut_end = 0
        #plot cutlines
        plt.plot(data)
        if cut_start>0:
            plt.plot([cut_start, cut_start], [-5000, 5000], 'k-')
        if cut_end>0:
            plt.plot([cut_end, cut_end],[-5000, 5000], 'k-')
        plt.show()
        #cut file
        if cut_start > 0 or cut_end > 0:
            for files in os.listdir(os.path.join(self.root, self.folder_amplipex)):
                if files.endswith(dataset + ".dat"):
                    filename = files        
            filepath = os.path.join(self.root, self.folder_amplipex, filename)
            os.rename(filepath,filepath + '_old')
            fid_old = open(filepath + '_old', "rb")
            fid_new = open(filepath, 'wb')
            fid_old.seek(cut_start*2*self.nb_channels)
            
            length = 2*(cut_end-cut_start)*self.nb_channels
            bufsize = 1024*1024
            while length:
                chunk = min(bufsize,length)
                data = fid_old.read(chunk)
                fid_new.write(data)
                length -= chunk
            fid_old.close()
            fid_new.close()

    

"""


>>> from neo.io.hdf5io import NeoHdf5IO
>>> iom = NeoHdf5IO('myfile.h5')
>>> iom

>>> iom.save(b)


>>> b.hdf5_path
'/block_0'



If you already have hdf5 file in NEO format, or you just created one, then you
may want to read NEO data (providing the path to what to read):

>>> b1 = iom.read_block("/block_0")
>>> b1
<neo.core.block.Block object at 0x34ee590>

or just use

>>> b1 = iom.get("/block_0")

You may notice, by default the reading function retrieves all available data,
with all downstream relations and arrays:

>>> b1._segments
[<neo.core.segment.Segment object at 0x34ee750>]
>>> b1._segments[0]._analogsignals[0].signal
array([  3.18987819e-01,   1.08448284e-01,   1.03858980e-01,
        ...
         3.78908705e-01,   3.08669731e-02,   9.48965785e-01]) * dimensionless

When you need to save time and performance, you may load an object without
relations

>>> b2 = iom.get("/block_0", cascade=False)
>>> b2._segments
[]

and/or even without arrays

>>> a2 = iom.get("/block_0/_segments/segment_0/_analogsignals/analogsignal_0",
lazy=True)
>>> a2.signal
[]

These functions return "pure" NEO objects. They are completely "detached" from
the HDF5 file - changes to the runtime objects will not cause any changes in the
file:

>>> a2.t_start
array(42.0) * ms
>>> a2.t_start = 32 * pq.ms
>>> a2.t_start
array(32.0) * ms
>>> iom.get("/block_0/_segments/segment_0/_analogsignals/analogsignal_0").t_start
array(42.0) * ms

However, if you want to work directly with HDF5 storage making instant
modifications, you may use the native PyTables functionality, where all objects
are accessible through "<IO_manager_inst>._data.root":

>>> iom._data.root
/ (RootGroup) 'neo.h5'
  children := ['block_0' (Group)]
>>> b3 = iom._data.root.block_0
>>> b3
/block_0 (Group) ''
  children := ['_recordingchannelgroups' (Group), '_segments' (Group)]

To understand more about this "direct" way of working with data, please refer to
http://www.pytables.org/
Finally, you may get an overview of the contents of the file by running

>>> iom.get_info()

"""