'''
Created on Oct 29, 2013

@author: mikkel
'''
from ..sorting import spike_sorting
from ..data import core
from ..data import basic
import os


class Recording(object):
    """
    Object representing one dataset
    """
    def __init__(self, par_file_path):
        """
        Initialise the parameters
        """
        self.par = {} 
        self.par['root'] = '/data/data/130904_peter'
        self.par['root_raw'] = '/data/data_raw/130904_peter'
        self.par['base_axon'] = '2013_09_04_'
	self.par['base_ampli'] = 'PeterP-2013-09-03_'
        self.par['labbook'] = {0:['0009','001'],
                1:['0011','002'],
                2:['0012','004'],
                3:['0013','005'],
                4:['0014','006'],
                5:['0015','007'],
                6:['0016','008'],
                7:['0017','009'],
                8:['0018','010'],
                9:['0019','011'],
                10:['0020','012'],
                11:['0021','013'],
                12:['0022','014'],
                13:['0023','015'],
                14:['0024','016']} #[data_nb, axon_db, ampli_nb]
        self.par['folder_ampli'] = "amplipex"
	self.par['folder_axon'] = 'axon'
        self.par['nb_channels'] = 193
        self.par['sr_ampli'] = 37132
	#eller self.par['sr_ampli'] = {0:39111,1:39250,2:39300}
	self.par['axon_sync'] = 7
	self.par['ampli_sync'] = 192
        self.par['klusta_sort_name'] = "sortering"
        self.par['probe_file'] = 'probedef.probe'
        self.par['probe_list'] = {1:0,2:64,3:128}
        par_file = self.read_par(par_file_path)
        self.par = dict(self.par.items() + par_file.items())
        #python 3: self.par = dict(list(self.par.items()) + list(par_file.items()))
        #make dictionaries out of single variables
	self.par['sr_ampli'] = self.make_par_dict('sr_ampli')
	self.par['axon_sync'] = self.make_par_dict('axon_sync')
	self.par['ampli_sync'] = self.make_par_dict('ampli_sync')
        #
        self.sort = spike_sorting.Spike_sorting(self)
        self.data = core.Data(self)
        
    def read_par(self, par_file_path):
        """
        Read the parameters from parameters file
        """
        parameters = {}
        parameters_name = os.path.join(par_file_path, "parameters.py")
        execfile(parameters_name, {}, parameters)
        return parameters

    def make_par_dict(self, name):
        """
        Convert variable to dictionary if it is not a dictionary already
        """
        if (isinstance(self.par[name],dict)):
	    return self.par[name]
        else:
	    tmp = {}
	    for data_nb in self.par['labbook']:
                tmp[data_nb] = self.par[name]
	    return tmp

    def list_sample_rate_ampli(self):
        """
        Print the sampling rate of amplipex for each trial
        """
        for data_nb in self.par['labbook']:
            axon_nb = self.par['labbook'][data_nb][0]
            ampli_nb = self.par['labbook'][data_nb][1]
	    filepath_axon = os.path.join(self.par['root_raw'],
	        self.par['folder_axon'],
		self.par['base_axon'] + axon_nb + '.abf')
	    filepath_ampli = os.path.join(self.par['root_raw'],
	        self.par['folder_ampli'],
		self.par['base_ampli'] + ampli_nb + '.dat')
            sr = basic.find_sample_rate_ampli(filepath_axon,
			    filepath_ampli,
			    self.par['axon_sync'][data_nb],
			    self.par['ampli_sync'][data_nb],
			    self.par['nb_channels'])
	    print(data_nb)
            print(sr)
