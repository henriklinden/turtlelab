# -*- coding: utf-8 -*-
"""
Created on Oct  29 2013

Functions for unit related plotting

@author: mikkel
"""
from __future__ import division
from matplotlib import pyplot as plt
plt.ion()
import numpy as np
from scipy import stats


def plot_raw_ampli(self, channels, dataset):
    """ Plot selected channels from amplipex dataset"""
    data = self.read_ch_ampli(channels, dataset)
    # Three subplots sharing both x/y axes
    f, axs = plt.subplots(len(channels), sharex=True, sharey=True)
    # loop through all plots
    for i in range(0,len(channels)):
        x = np.arange(0,data.shape[0])
        #axs[i].plot(x/self.sr,data)
        axs[i].plot(x,data[:,i])
        axs[i].set_title(str(channels[i]))
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.show()
        
def plot_raw(self, channels, dataset):
    """ Plot selected channels from dat-file"""
    # Three subplots sharing both x/y axes
    f, axs = plt.subplots(len(channels), sharex=True, sharey=True)
    # loop through all plots
    for i in range(0,len(channels)):
        data = self.read_ch_dat(channels[i], dataset)
        x = np.arange(0,data.shape[0])
        #axs[i].plot(x/self.sr,data)
        axs[i].plot(x,data)
        axs[i].set_title(str(channels[i]))
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.show()
    
def plot_filtered(self, channels, dataset):
    """ Plot filtered and raw data from selected channels from dat-file"""
    # Three subplots sharing both x/y axes
    f, axs = plt.subplots(len(channels),2, sharex=True)#, sharex=True, sharey=True)
    #loop through all plots
    for i in range(0,len(channels)):
        data = self.read_ch_dat(channels[i], dataset)
        data_filt = self.filtered(data)
        x = np.arange(0,data.shape[0])
        x = x/self.sr
        axs[i,0].plot(x,data)
        axs[i,0].set_title(str(channels[i]))
        axs[i,1].plot(x,data_filt)
        axs[i,1].set_title(str(channels[i]))
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.show()
    
def plot_axon(self):
    """ Plot selected channels from axon dataset"""
    pass
    """# Plexon files
distantfile = 'https://portal.g-node.org/neo/plexon/File_plexon_3.plx'
localfile = './File_plexon_3.plx'
urllib.urlretrieve(distantfile, localfile)

#create a reader
reader = neo.io.PlexonIO(filename = 'File_plexon_3.plx')
# read the block
bl = reader.read(cascade = True, lazy = False)
print bl
# acces to segments
for seg in bl.segments:
    print seg
    for asig in seg.analogsignals:
        print asig
    for st in seg.spiketrains:
        print st
        """
        
def filtered(self, data):
    """ Filter data for plotting"""
    stop_freq = np.array([500,10000])
    pass_freq = np.array([700,8000])
    norm_pass = pass_freq/(self.sr/2)
    norm_stop = stop_freq/(self.sr/2)
    (N, Wn) = signal.buttord(wp=norm_pass, ws=norm_stop, gpass=3, gstop=20, analog=0)
    (b, a) = signal.butter(N, Wn, btype='bandpass', analog=0, output='ba')
    data_filt = signal.filtfilt(b, a, data)
    return data_filt

def plot_info_spikes(self):
    pass

def plot_spikes_forms(self):
    pass
    
    """
    function spik_plot_spikeforms(spik_sorting,which_spikes,cluster_nb)
%SPIK_PLOT_SPIKEFORMS Plot the spike forms of one cluster

%cut out sniplet
cut_bef = spik_sorting.par.cut_bef; %samples before
cut_aft = spik_sorting.par.cut_aft; %samples after
ch = spik_sorting.ch;

spikes = spik_sorting.events_x(which_spikes);


handle_raw = figure();
handle_hist = figure('Position',[1 1 250 700]);

index = 1;
for k= ch
        t = [];
        vm = [];
        vm_raw = [];
        snippets = spik_sorting.snippet(k).data(which_spikes,:);
        %figure(handle_raw);
        %subplot(length(ch),1,index);
        %hold all
        %title(['Cluster nb: ' int2str(cluster_nb) ' , ch nb: ' int2str(k) ' , spikes: ' int2str(length(spikes))])
        
        tmp_length = size(spikes,2);
        if (tmp_length==0) 
            continue
        end
        idx = rand(tmp_length,1)>(1-1000/tmp_length);
        
        i = 1;
        for xs = spikes(idx)
            tmp = snippets(i,:);
            tmp = tmp-mean(tmp,2);
            %plot(tmp);
            %xlim([1 (cut_bef+1+cut_aft)])
            t = [t (-cut_bef):(cut_aft)];
            vm = [vm tmp];
            vm_raw = [vm_raw, tmp'];
            i = i+1;
            break
        end       
        X = [vm', t'];
        
        figure(handle_hist);
        %colormap hot
        n = hist3(X,[(cut_bef+1+cut_aft) (cut_bef+1+cut_aft)]);
        %xlabel('MPG'); ylabel('Weight');
        %set(gcf,'renderer','opengl');
        %set(get(gca,'child'),'FaceColor','interp','CDataMode',...
        %'auto');
        subplot(length(ch),1,index);
        hold on
        xlim([1 (cut_bef+1+cut_aft)])
        ylim([1 (cut_bef+1+cut_aft)])
        image(n);
        title(['Cluster nb: ' int2str(cluster_nb) ' , ch nb: ' int2str(k) ' , spikes: ' int2str(length(spikes))])
        

        figure(handle_raw);
        ax(index) = subplot(length(ch),1,index);
        hold on
        sp = spaps((1:length(vm_raw)),vm_raw,0.0000000001);  fnplt(sp);
        plot(vm_raw,'.');
        xlim([1 (cut_bef+1+cut_aft)])
        linkaxes(ax,'y');
        index = index+1;
end"""
    
def plot_isi(self):
    pass
    
    """
    function spik_plot_isi_dist(spik_sorting,which_spikes)
%SPIK_PLOT_ISI_DIST plot interspike interval distribution
%   2011 Mikkel Vestergaard


%% ISI
maxlag = 1000;
bins = round(maxlag/10);

spikes = spik_sorting.events_x(which_spikes)./(spik_sorting.dset.info.sr_multi/1000);

isis = diff(spikes);   
isis = isis(isis <= maxlag); %remove long isi 
figure
[n, xout] = hist(isis,linspace(0,maxlag,bins));

xlabel('Interspike interval')
ylabel(['No. of spikes ' int2str(length(spikes))])
xlim([0 maxlag])


ymax = max(n);
refrac = spik_sorting.par.refrac; %refractory period
remove_sp = spik_sorting.par.remove_sp; %time where two spikes are not possible due to the analyses

%indicate period where two spikes not possible due to the analyses
patch([0 remove_sp remove_sp 0 ], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none'); %gray box
%indicate refractory period
patch([remove_sp [refrac refrac] remove_sp ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none'); %red box
hold on
bar(xout,n);

        
end
"""
    
def plot_corr(self):
    pass

def plot_features(self):
    pass

def plot_raw_data_spikes(self):
    pass

def plot_filtered_data_spikes(self):
    pass
    
    """
    function spik_plot_spikes(dset,d_nb, ch)
% SPIK_PLOT_SPIKES obsolete function


%% setup / parameters

srCh = dset.info.sr_multi;

arg = dset.d{d_nb};

arg.x1 = 8;
arg.x2 = 25;

%x values
xCh = [(arg.x1*srCh+1):(arg.x2*srCh)];

%% load data

%load Files.txt
files = textread(arg.path.multi_index,'%s');

nb_files = length(ch);
    
data_raw = struct([]);
for k= ch
    file_to_read = files(k)
    fullpath = fullfile(arg.path.multi, char(file_to_read));
    data_raw(k).file = load(fullpath);
    
end

%% filter data
data = struct([]);
for k = ch
    data(k).filter = filtering(data_raw(k).file.data,'butter');
    noise = var(data(k).filter(3*srCh:7*srCh));
    
    data(k).filter = wiener2(data(k).filter,[1 20],noise);
end

% plot un/filtered data
    
figure;


nb_plots = nb_files; %total number of plots

index = 1;

for k= ch
    
    %plot raw
    ax(index)=subplot(nb_plots,2,index);
    hold on;
    plot((xCh/srCh),data_raw(k).file.data(xCh),'b');
    ylabel(mat2str(k));
    set(gca,'xtick',[],'ytick',[]);
    hold off;
    
    index = index +1;
    
    %plot filterd
    ax(index)=subplot(nb_plots,2,index);
    hold on;
    plot((xCh/srCh),data(k).filter(xCh),'b');
    ylabel(mat2str(k));
    set(gca,'xtick',[],'ytick',[]);
    hold off;
 
    index = index +1;
end
        
linkaxes(ax,'xy');
xlim([arg.x1 arg.x2]);

clear data_raw ax;

%% threshold

xmin_all = []; %all xmin

for k = ch
    %determine threshold level
    data(k).thres = -6*median(abs(data(k).filter(:)))/0.6745;
    %find indicies when voltage return from triggered level
    returns = (find(diff(data(k).thres>data(k).filter(:))==-1)+1)';
    data(k).returns = returns;
    xmin = [];
    %ymin = [];
    x1 = 1;
    for x2 = returns
        [y,x] = min(data(k).filter(x1:x2));
        %ymin = [ymin y];
        xmin = [xmin (x+x1-1)];
        x1 = x2;
    end
    data(k).xmin = xmin;
    xmin_all = [xmin_all xmin];
    %data(k).ymin = ymin;
end

%remove multi-triggered events
%OBS: burde nok fjerne den event med lavest y-værdi
xmin_all = sort(xmin_all);
xmin_all = xmin_all([1 (diff(xmin_all)>8)]>0);
    

%%    
figure;

nb_plots = nb_files; %total number of plots

index = 1;

for k= ch
       
    %plot filterd
    ax(index)=subplot(nb_plots,1,index);
    hold on;
    plot((xCh/srCh),data(k).filter(xCh),'b');
    %threshold
    line([arg.x1 arg.x2],[data(k).thres data(k).thres]);
    %local minima
    %plot((data(k).xmin/srCh),data(k).filter(data(k).xmin/srCh),'r*');
    plot((xmin_all/srCh),data(k).filter(xmin_all),'r*');
    
    %plot((data(k).returns/srCh),data(k).thres,'g*');
    ylabel(mat2str(k));
    set(gca,'ytick',[]);
    hold off;
 
    index = index +1;
end
        
linkaxes(ax,'xy');
xlim([arg.x1 arg.x2]);

%% Feature extraction (PCA)

%cut out sniplet
cut_bef = 15; %samples before
cut_aft = 20; %samples after

for k= ch
    sniplet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    i = 1;
    for xs = xmin_all
        sniplet(i,:) = data(k).filter((xs-cut_bef):(xs+cut_aft));
        i = i+1;
    end
    i
    
    data(k).sniplet = sniplet;
    %plot sniplets
    figure;
    plot(sniplet');
    
    %PCA
    [coefs,scores,variances,tsquare] = princomp(sniplet);

    data(k).scores = scores(:,1:3);
    
    percent_explained = (100*variances)/sum(variances);
    
    figure;
    pareto(percent_explained)
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    
    figure;
    biplot(coefs(:,1:3), 'scores',scores(:,1:3));
    axis([-.26 1 -.51 .51 -.61 .81]);
    view([30 40]);
end

%% save data

nbDimensions = 3;

%the fet file
fid = fopen('NAVN.fet.1', 'w');

% print number of feature columns
fprintf(fid, [int2str(nbDimensions*length(ch)+1) '\n']);

fclose(fid);

scores = [];

for k = ch
    scores = [scores data(k).scores];
end

feat_scaling_factor = 32000/max(max(max(scores)),-min(min(scores)));
scores = scores*feat_scaling_factor;
scores = round(scores);
scores = [scores xmin_all'];


    
dlmwrite('NAVN.fet.1', scores, '-append', 'delimiter', ' ','precision', '%i')

%the spk file

fid = fopen('NAVN.spk.1', 'w');

for k = ch
    wave_scaling_factor = 32000/max(max(max(data(k).sniplet)),-min(min(data(k).sniplet)));
    data(k).sniplet = round(data(k).sniplet*wave_scaling_factor);
end

t = 1;
%for all spikes/sniplets
for i = 1:length(xmin_all)
%for i = 1:3
    tmp = [];
    for k = ch
        tmp = [tmp; data(k).sniplet(i,:)];
        %data(k).sniplet(i,:)
    end
    for j = 1:numel(tmp)
        fwrite(fid, tmp(j),'int32');
        %tmp(j)
        t = t+1;
        %[num2str(tmp(j)) '\n']
    end
    %fprintf(fid, '%f\n',tmp);
end
t
fclose(fid);

fid = fopen('NAVN.dat', 'w');
tmp = [];
for k = ch
    tmp =  [tmp;data(k).filter];
end
wave_scaling_factor = 16000/max(max(max(tmp)),-min(min(tmp)));
tmp = round(tmp*wave_scaling_factor);
fwrite(fid, tmp,'int32');
fclose(fid);


end"""
