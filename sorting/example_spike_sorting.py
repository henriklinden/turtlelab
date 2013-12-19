# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 17:29:01 2013

@author: mikkel
"""

import spike_sorting as sort

#TODO: Skal opdateres...
        
sorting = sort.spike_sorting()
sorting.root = "/data/data/130904_peter"
sorting.probe_list = {1:0,2:64,3:128}
sorting.load_data_ampli()
sorting.plot_raw_ampli([0,2,3,192],'010')
#sorting.cut_ampli()
#sorting.save_combined_dat(range(0,64), '004')
#sorting.save_klusta_files()
#sorting.run_spikedetekt()
sorting.run_klustakwik()


"""
spik_sorting.features = features;
spik_sorting.clusters = clusters;
spik_sorting.models = models;
spik_sorting.logp = logp;

%% Display the number of spikes in each cluster

%% Plot the spike forms of one cluster

%% Plot the spike forms of all clusters

%% Plot the ISI distribution of one cluster


%% Plot correlation diagram




%% Plot the features of all clusters 


eval = spik_evaluate_cluster(spik_sorting, 4,1);


"""