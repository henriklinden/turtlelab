'''
Created on Oct 29, 2013

@author: mikkel
'''
from ..sorting import spike_sorting
from ..data import core
import os


class Recording(object):
    
    def __init__(self, par_file_path):           
        self.par = {} 
        self.par['root'] = '/data/data/130904_peter'
        self.par['base_axon'] = '2013_09_04'
        self.par['labbook'] = [[0,'0009','001'],
                    [1,'0011','002'],
                    [2,'0012','004'],
                    [3,'0013','005'],
                    [4,'0014','006'],
                    [5,'0015','007'],
                    [6,'0016','008'],
                    [7,'0017','009'],
                    [8,'0018','010'],
                    [9,'0019','011'],
                    [10,'0020','012'],
                    [11,'0021','013'],
                    [12,'0022','014'],
                    [13,'0023','015'],
                    [14,'0024','016']] #[data_nb, axon_db, ampli_nb]
        self.par['folder_amplipex'] = "amplipex"
        self.par['nb_channels'] = 193
        self.par['sr_ampli'] = 37132
        self.par['klusta_sort_name'] = "sortering"
        self.par['probe_file'] = 'probedef.probe'
        self.par['probe_list'] = {1:0,2:64,3:128}
        par_file = self.read_par(par_file_path)
        self.par = dict(self.par.items() + par_file.items())
        #python 3: self.par = dict(list(self.par.items()) + list(par_file.items()))
        
        self.sort = spike_sorting.Spike_sorting(self)
        self.data = core.Data(self)
        
    def read_par(self, par_file_path):
        parameters = {}
        parameters_name = os.path.join(par_file_path, "parameters.py")
        execfile(parameters_name, {}, parameters)
        return parameters