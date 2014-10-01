"""
Basic functions for reading and writing to files
@author mikkel
"""
import os
import numpy as np
from neo import io
from matplotlib import pyplot as plt
def split_ampli_dats(self):
    """
    Split all .dat files in folder
    TODO: fjern referencer til 'self'
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
    Splits one .dat file into one file for each channel
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

def read_ch_ampli(channels, filepath, nb_channels):
    """
    Load the specified channels from a dataset
    :param list channels: channels to load
    :param string filepath: the full path to the dataset
    :param int nb_channels: the total number of channels in file
    """
    #    #find the .dat file to load
    #    for files in os.listdir(os.path.join(self.root, self.folder_amplipex)):
    #        if files.endswith(dataset + ".dat"):
    #            filename = files
    #    filepath = os.path.join(self.root, self.folder_amplipex, filename)
    # load file
    fid = open(filepath, "rb")
    length = os.fstat(fid.fileno()).st_size/(2*nb_channels)
    # read into the 'data' var
    data = np.zeros((length,len(channels)),dtype='int16')
    idx = 0
    for tmp_data in _read_chunk_ampli_dat(fid,nb_channels):
        tmp_data = np.reshape(tmp_data,(-1,nb_channels))
        data[idx:(idx+tmp_data.shape[0]),:] = tmp_data[:,channels]
        idx += tmp_data.shape[0]
    return data
        
def _read_chunk_ampli_dat(file_object, nb_channels):
    """
    Generator for reading the ampli dat file in chunks
    """
    while True:
        data = np.fromfile(file_object,'int16',1000000*nb_channels)       
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

def write_vars_to_file(fid, variables):
    """ Write list as to (parameter) file"""
    for (name, val) in variables.items():
        fid.write("%s = %s\n" % (name, repr(val)))
        
def read_raw_dat_names(folderPath):
    """ Find all file names of .dat files"""
    raw_data = []
    #folderPath = os.path.join(self.root, self.folder_amplipex)
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

def find_sync_samples_ampli(data, sync_method=0):
    """ Find the number of samples between the trigger signals (ampli)"""
    if sync_method == 0:
        #Alex sync: one pulse marks beginning; one pulse marks end
        idx = np.where(data>5000)[0]
        start = idx[0]
        end = idx[np.where(idx>(start+100000))[0][0]] #100000= couple of secs before
            #checking for crossings
    elif sync_method == 1:
        #Peter sync: one pulse marks whole recording
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
    return (end-start)

def read_axon(filename):
    """ Read axon .abf file"""
    r = io.AxonIO(filename)
    bl = r.read_block(lazy=False, cascade=True)
    return bl

def find_sync_time_axon(data, sync_method=0):
    """Find time between the trigger signals (axon)"""
    if sync_method == 0:
        #Alex sync: one pulse marks beginning, one pulse marks end
        idx = np.where(data>3)[0]
        start = idx[0]
        end = idx[np.where(idx>(start+100000))[0][0]] #100000= couple of secs before
            #checking for crossings
        time = data.times[end]-data.times[start]
    elif sync_method == 1:
        #Peter sync: one pulse marks whole recording
        idx = np.where(data>3)[0]
        time = data.times[idx[-1]]-data.times[idx[0]]
    return time

def find_sample_rate_ampli(filepath_axon, filepath_ampli, ch_sync_axon, ch_sync_ampli, nb_channels):
    """ Find the sample rate of the amplipex file"""
    #load axon
    block = read_axon(filepath_axon)
    #load ampli
    ampli_data = read_ch_ampli([ch_sync_ampli],filepath_ampli,nb_channels)
    #plot raw data
    #plt.figure()
    #plt.plot(block.segments[0].analogsignals[ch_sync_axon])
    #plt.plot(ampli_data)
    #time axon
    time = find_sync_time_axon(block.segments[0].analogsignals[ch_sync_axon])
    #samples ampli
    nb_samples = find_sync_samples_ampli(ampli_data)
    #sample rate ampli
    #print(nb_samples)
    #print("time {} and samples {}".format(time,nb_samples))
    return nb_samples/time
