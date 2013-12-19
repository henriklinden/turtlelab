# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 18:50:42 2013

@author: mikkel
"""
from __future__ import division
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
import StringIO
import subprocess
import stat
import math
from scipy import signal

class Spike_sorting:    
    
    #parameters
    #thres
    #cut_bef
    #cut_after
    #refrac
           
    def __init__(self, rec):
        self.rec = rec
        self.root = rec.par['root']
        self.folder_amplipex = rec.par['folder_amplipex']
        self.nb_channels = rec.par['nb_channels']
        self.sr = rec.par['sr_ampli']
        self.klusta_sort_name = rec.par['klusta_sort_name']
        self.probe_file = rec.par['probe_file']
        self.probe_list = rec.par['probe_list']
        
    
    
    def save_klusta_files(self):
        filepath = os.path.join(self.root,self.folder_amplipex,self.klusta_sort_name + '.params')
        fid = open(filepath,'w')
        params = {}        
        params.update({'RAW_DATA_FILES':self.read_raw_dat_names()})
        params.update({'SAMPLERATE':self.sr})
        params.update({'NCHANNELS':self.nb_channels})
        params.update({'PROBE_FILE':self.probe_file})
        # Output directory, files are inserted in OUTPUT_DIR/OUTPUT_NAME
        #OUTPUT_DIR = None # the output directory, use params directory if None
        #OUTPUT_NAME = None # the filename for created directories, use params filename if None
        params.update({'OUTPUT_DIR':os.path.join(self.root,'sortering')})
        params.update({'OUTPUT_NAME':'sortering'})        
        #CHUNKS_FOR_THRESH = 15 # number of chunks used to determine threshold for detection
        params.update({'CHUNKS_FOR_THRESH':15})
        #THRESH_SD = 4.5 # threshold for detection. standard deviations of signal
        params.update({'THRESH_SD':4.5})
        # Options for computing in chunks
        #CHUNK_SIZE = 40000   # number of time samples used in chunk for filtering and detection
        #CHUNK_OVERLAP_SECONDS = 0.01 # overlap time (in seconds) of chunks, should be wider than spike width
        params.update({'CHUNK_SIZE':40000})
        params.update({'CHUNK_OVERLAP_SECONDS':0.1})
        self.write_vars_to_file(fid,params)        
        
        filepath = os.path.join(self.root,self.folder_amplipex,self.probe_file)
        fid = open(filepath,'w')
        probes = self.generate_probes()
        self.write_vars_to_file(fid, {'probes':probes})
    
    
        
    def generate_probes(self):
        probes = {}
        for probe_nb, probe_ch in self.probe_list.items():
            probes.update(self.generate_probe(probe_ch))  
        i = 1
        for shank_nb, shank_chs in probes.items():
            probes[i] = probes.pop(shank_nb)
            i += 1
        return probes        
        
    def generate_probe(self, ch_offset):
        shanks = []
        shanks.append([31,18,27,29,26,21,23,22])
        shanks.append([25,37,16,17,20,35,33,19])
        shanks.append([32,46,39,38,34,44,41,36])
        shanks.append([40,24,42,30,47,28,45,43])
        shanks.append([6,54,0,52,2,49,53,51])
        shanks.append([48,62,56,57,50,60,58,55])
        shanks.append([59,7,15,14,61,10,13,63])
        shanks.append([12,1,3,5,11,4,8,9])
        
        probe = {}
        j = 1
        for shank in shanks:
            ch_map = []
            for i in range(0,len(shank)-1):
                ch_map.append(((shank[i]+ch_offset),(shank[i+1]+ch_offset)))
                if i<(len(shank)-2):
                    ch_map.append(((shank[i]+ch_offset),(shank[i+2]+ch_offset)))
            probe.update({(str(ch_offset) + str(j)):ch_map})
            j += 1
        return probe
        
    def run_spikedetekt(self):
        filepath = os.path.join(self.root,self.folder_amplipex,self.klusta_sort_name + '.params')
        #brug Popen... 
        folder = os.path.join(self.root,self.folder_amplipex)
        #os.system('python ~/Documents/programmer/klusta-team/spikedetekt/scripts/detektspikes.py ' + filepath)
        cmd = 'python ~/Documents/programmer/klusta-team/spikedetekt/scripts/detektspikes.py ' + filepath        
        #subprocess.Popen(['terminator','-e', cmd],cwd = folder)
        filename = os.path.join(self.root,'runSpikeDetekt.sh')
        fid = open(filename,'w')
        fid.write('cd ' +folder + '\n')
        fid.write(cmd + '\n')
        fid.close()
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IEXEC)        
        
    def run_klustakwik(self):
        #fPath = os.path.join(self.root, self.folder_amplipex,self.probe_file)
        folderPath = os.path.join(self.root,'sortering')
        filename_sh = os.path.join(self.root,'runKlustaKwik.sh')
        fid = open(filename_sh,'w')
        fid.write('cd ' +folderPath + '\n')
        for root, subFolders, files in os.walk(folderPath):
            del subFolders
            for filename in files:
                print filename
                file_nb = filename.split('.')[-1]
                extension = filename.split('.')[-2]
                if (extension=="clu"):
                    cmd = '~/Documents/programmer/klusta-team/klustakwik/MaskedKlustaKwik ' \
                        + self.klusta_sort_name + ' ' + file_nb + ' ' + '-DropLastNFeatures 1 ' \
                        +'-MaskStarts 300 ' \
                        +'-MaxPossibleClusters 500 ' \
                        +'-PenaltyK 1.0 ' \
                        +'-PenaltyKLogN 0.0 ' \
                        +'-SplitFirst 20 ' \
                        +'-SplitEvery 100 ' \
                        +'-UseDistributional 1 '
                    #subprocess.Popen(['terminator','-e',cmd],cwd = folderPath)
                    fid.write(cmd + '\n')
        fid.close() 
        st = os.stat(filename_sh)
        os.chmod(filename_sh, st.st_mode | stat.S_IEXEC)
                    
         
    def read_config(self, fid):
        filepath = os.path.join(self.root,self.folder_amplipex,self.klusta_sort_name + '.params')
        config = ConfigParser.ConfigParser()
        config.read('example.cfg')
        ini_str = '[root]\n' + open(filepath, 'r').read()
        ini_fp = StringIO.StringIO(ini_str)
        config = ConfigParser.RawConfigParser()
        config.readfp(ini_fp)
        return config.get('root', 'RAW_DATA_FILES') 
            
    def make_split_file(self):
        filename_log = os.path.join(self.root,'sortering','sortering.log')
        fid_log = open(filename_log,'r')
        filename_splt = os.path.join(self.root,'sortering','split_file.txt')
        fid_splt = open(filename_splt,'w')
        ids = eval(self.read_config(fid_log))
        for dat_file in ids:
            print dat_file
            filepath = os.path.join(self.root, self.folder_amplipex, dat_file)
            fid = open(filepath, "rb")
            length = os.fstat(fid.fileno()).st_size/(2*self.nb_channels)
            fid.close()
            fid_splt.write(dat_file + ' ' + str(int(length)) + '\n')
            
        fid_log.close()
        fid_splt.close()
        
    
    
    
    
    def features_pca(self):
        pass
    
    def features_wavelets(self):
        pass
    
    def features_spike_detect(self):
        pass
    
    def sorting_etos(self):
        pass
    
    def sorting_klustakwik(self):
        pass
    
    def fix_wavedec(self):
        pass
    """function [c,l] = fix_wavedec(x,n)
% Does a haar wavelet decomposition with n scales.
% Avoids using wavelet toolbox


Lo_D = [ 0.7071 0.7071];
Hi_D = [-0.7071 0.7071];

s = size(x); x = x(:)'; 
c = []; 
l = [length(x)];

dwtEXTM = 'sym';
shift = 0;

for k = 1:n
    lf = length(Lo_D);
    lx = length(x);
    lenEXT = lf-1; lenKEPT = lx+lf-1;
   
    I = getSymIndices(lx,1);
    y  = x(I);
    
    x = convdown(y,Lo_D,lenKEPT,shift);
    d = convdown(y,Hi_D,lenKEPT,shift);
    
    c     = [d c];            % store detail
    l     = [length(d) l];    % store length
end

% Last approximation.
c = [x c];
l = [length(x) l];

if s(1)>1, c = c'; l = l'; end


%-----------------------------------------------------%
% Internal Function(s)
%-----------------------------------------------------%
function y = convdown(x,f,lenKEPT,shift)

y = conv2(x(:)',f(:)'); if size(x,1)>1 , y = y'; end

sx = length(y);
begInd = 1;
[first,last,ok] = GetFirstLast(sx,begInd,lenKEPT);
if ok , y = y(first(1):last(1)); end

y = y(2-rem(shift,2):2:end);

%-----------------------------------------------------%
%----------------------------------------------------------------------------%
function I = getSymIndices(lx,lf)

I = [lf:-1:1 , 1:lx , lx:-1:lx-lf+1];
if lx<lf
    K = (I<1);
    I(K) = 1-I(K);
    J = (I>lx);
    while any(J)
        I(J) = 2*lx+1-I(J);
        K = (I<1);
        I(K) = 1-I(K);
        J = (I>lx);
    end
end
%----------------------------------------------------------------------------%
%----------------------------------------------------------------------------%
function [first,last,ok] = GetFirstLast(sx,begInd,varargin)

oneDIM = isequal(begInd,1);
s = varargin{1}(:)';
if ~oneDIM
    K  = find(s>sx);
    s(K) = sx(K);
    m = find((s < 0) | (s ~= fix(s)));
    ok = isempty(m);
else
    ok = (s>=0) & (s<sx) & (s == fix(s));
end
if ok==0 , first = begInd; last = s; return; end

nbarg = length(varargin);
if nbarg<2, o = 'c'; else , o = lower(varargin{2}); end

err = 0;
if ischar(o(1))
    switch o(1)
        case 'c'
            d = (sx-s)/2;
            if nbarg<3
                if length(o)>1 , side = o(2:end); else , side = 'l'; end
            else
                side = varargin{3};
            end
            if oneDIM
                [first,last] = GetFirst1D(side,sx,d);
            else
                if length(side)<2 , side(2) = 0; end
                for k = 1:2
                    [first(k),last(k)] = GetFirst1D(side(k),sx(k),d(k));
                end
            end

        case {'l','u'} , first = begInd; last = s;
        case {'r','d'} , first = sx-s+1; last = sx;
        otherwise      , err = 1;
    end
else
    first = o; last = first+s-1;
    if ~isequal(first,fix(first)) | any(first<1) | any(last>sx)
        err = 1;
    end
end
if err
    errargt(mfilename,'invalid argument','msg');
    error('*');
end
%----------------------------------------------------------------------------%
function [first,last] = GetFirst1D(side,s,d)

switch side
  case {'u','l','0',0} , first = 1+floor(d); last = s-ceil(d);
  case {'d','r','1',1} , first = 1+ceil(d);  last = s-floor(d);
  otherwise    , first = 1+floor(d); last = s-ceil(d);  % Default is left side
end
%----------------------------------------------------------------------------%
"""



    def spik_add_clu(self):
        pass
    """
    % split cluster file from klusters
function spik_add_clu(name)
fid = fopen([name '.clu.1']);
clusters = textscan(fid, '%u');
fclose(fid);

clusters = cell2mat(clusters);

load([name '_spikes.mat'])
spikes(:,3) = clusters(2:end);
save([name '_spikes.mat'], 'spikes')

end"""

def feature_pca(self):
    pass

"""%find features by principle component analysis
%26/1/2011 Mikkel Vestergaard

function spik_sorting = spik_feature_pca(spik_sorting, show_plot)

%dset = spik_sorting.dset;
%d_nb = spik_sorting.d_nb;
ch = spik_sorting.ch;

%cut out sniplet
%cut_bef = 40; %samples before
%cut_aft = 30; %samples after

cut_bef = spik_sorting.par.cut_bef;
cut_aft = spik_sorting.par.cut_aft;

xmin_all = spik_sorting.events_x;

for k= ch
    k
    %snippet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    snippet = spik_sorting.snippet(k).data;
    i = 1;
    if (show_plot)
        figure
        hold on
        %plot(data(k).data);
    end
    for xs = xmin_all
        %snippet(i,:) = data(k).data((xs-cut_bef):(xs+cut_aft));
        i = i+1;
        if (show_plot)
            plot((xs-cut_bef):(xs+cut_aft),snippet(i,:),'r');
        end
    end
    i
    
    spik_sorting.snippet(k).data = snippet;
    %plot sniplets
    if (show_plot) %&(k==ch(1)))
        figure
        subplot(2,1,1);
        plot(snippet');
        subplot(2,1,2);
        temp = var(snippet);
        plot(temp);
    end
    
end

%PCA
%snippet = [];%zeros(size(xmin_all,1)*length(ch),(cut_bef+cut_aft+1));
snippet = [spik_sorting.snippet(ch).data];

[coefs,scores,variances,tsquare] = princomp(snippet);

spik_sorting.scores.data = scores(:,1:16);

if (show_plot)
    percent_explained = (100*variances)/sum(variances);
    figure
    pareto(percent_explained)
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    
    figure
    biplot(coefs(:,1:3), 'scores',scores(:,1:3));
    axis([-.26 1 -.51 .51 -.61 .81]);
    view([30 40]);
end"""

def features_wl(self):
    pass

"""
%find features by wavelets
%26/1/2011 Mikkel Vestergaard

function [spik_sorting] = spik_feature_wl(spik_sorting, show_plot)

%dset = spik_sorting.dset;
%d_nb = spik_sorting.d_nb;
ch = spik_sorting.ch;

cut_bef = spik_sorting.par.cut_bef; %samples before
cut_aft = spik_sorting.par.cut_aft; %samples after

xmin_all = spik_sorting.events_x;

temp_features = [];

for k= ch
    k
    snippet = spik_sorting.snippet(k).data;
    %snippet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    wl_features = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    i = 1;
    %if (show_plot)
        %figure
        %hold on
        %plot(data(k).filter);
    %end
    for xs = xmin_all
        %snippet(i,:) = data(k).filter((xs-cut_bef):(xs+cut_aft));
        %Wavelets
        %[c,l] = fix_wavedec(snippet(i,:),4);
        %snippet(i,:) = data(1,k).data((xs-cut_bef):(xs+cut_aft));
        
        %[c,l]=wavedec(snippet(i,:),4,'haar');                    
        c = waveletcdf97(snippet(i,:),(cut_bef+cut_aft+1));
        wl_features(i,:) = c;%(1:(cut_bef+cut_aft+1));
        i = i+1;
        if (show_plot)
            plot((xs-cut_bef):(xs+cut_aft),snippet(i,:),'r');
        end
    end
    i
    
    %spik_sorting.sniplet(k).data = snippet;
    temp_features(k).wl_features = wl_features;
    
end

% vælg features udfra lillifors KS-test
% for k=ch
%     sniplet = data(k).sniplet;
%     wl_features = data(k).wl_features;
%     %plot sniplets
%     if (show_plot) %&(k==ch(1)))
%         figure
%         subplot(2,1,1);
%         plot(sniplet');
%         subplot(2,1,2);
%         temp = var(sniplet);
%         plot(temp);
%     end
%     
%     %Lillifors KS-test
%     leng = size(wl_features,2);
%     kstats = zeros(1,leng);
%     for i = 1:leng
%         %[h,p,kstat,critval] = lillietest(wl_features(:,i));
%         %test_ks = wl_features(:,i);
%         
%         thr_dist = std(wl_features(:,i)) * 3;
%             thr_dist_min = mean(wl_features(:,i)) - thr_dist;
%             thr_dist_max = mean(wl_features(:,i)) + thr_dist;
%             aux = wl_features(find(wl_features(:,i)>thr_dist_min & wl_features(:,i)<thr_dist_max),i);
%             
%             if length(aux) > 10;
%                 [ksstat]=test_ks(aux);
%                 kstats(i)=ksstat;
%             else
%                 kstats(i)=0;
%             end
%         
%         
%         
%         %kstats(i) = kstat;
%     end
%     
%     [kstats,inx] = sort(kstats,'descend');
%     
%     %plot 10 best
%     figure
%     for i = 1:16
%         subplot(8,2,i);
%         hold on
%         temp = wl_features(:,inx(i));
%         hist(temp,30);
%         set(gca,'xtick',[],'ytick',[])
%         title([num2str(kstats(i)) ' , ' int2str(inx(i))])
%     end
%     
% end

% vælg features udfra anden test
spik_sorting.scores.data =[];
for k=ch
    snippet = spik_sorting.snippet(k).data;
    wl_features = temp_features(k).wl_features;
    %plot snippets
    if (show_plot) %&(k==ch(1)))
        figure
        subplot(2,1,1);
        plot(snippet');
        subplot(2,1,2);
        temp = var(snippet);
        plot(temp);
    end
    
    %Lillifors KS-test
    leng = size(wl_features,2);
    kstats = zeros(1,leng);

    for i = 1:leng
        %[h,p,kstat,critval] = lillietest(wl_features(:,i));
        %test_ks = wl_features(:,i);
        
        thr_dist = std(wl_features(:,i)) * 3;
        thr_dist_min = mean(wl_features(:,i)) - thr_dist;
        thr_dist_max = mean(wl_features(:,i)) + thr_dist;
        aux = wl_features(find(wl_features(:,i)>thr_dist_min & wl_features(:,i)<thr_dist_max),i);
        
        if length(aux) > 10;
            f1 = gmdistribution.fit(aux,1); % single component, unimodal
            f2 = gmdistribution.fit(aux,2); % two components, bimodal
            
            %[ksstat]=test_ks(aux);
            ksstat = f1.AIC-f2.AIC;
            %ksstat = f2.AIC;
            kstats(i)=ksstat;
        else
            kstats(i)=0;
        end
        
        
        
        %kstats(i) = kstat;
    end
    
    [kstats,inx] = sort(kstats,'descend');
    
    %plot 10 best
    if (show_plot) 
        figure
        for i = 1:16
            subplot(8,2,i);
            hold on
            temp = wl_features(:,inx(i));
            hist(temp,30);
            set(gca,'xtick',[],'ytick',[])
            title([num2str(kstats(i)) ' , ' int2str(inx(i))])
        end
    end
    
    spik_sorting.scores.data =[spik_sorting.scores.data wl_features(:,inx(1:4))];
end


end"""

def sort_wl(self):
    pass

"""function [data, xmin_all] = spik_sort_wl(dset,d_nb, ch, name, show_plot)
%SPIK_SORT_WL obsolete function
%% setup / parameters

srCh = dset.info.sr_multi;

%x values for graphs
%plot_x1 = 9; plot_x2 = 20;
%xCh = (plot_x1*srCh+1):(plot_x2*srCh);

%% load data



nb_files = length(ch);
    
data_raw = struct([]);
for k = ch
    data_raw(k).data = [];
end
data_breaks = [];

for j = d_nb
    arg = dset.d{d_nb};
    %load Files.txt
    files = textread(arg.path.multi_index,'%s');
    for k = ch
        file_to_read = files(k);
        fullpath = fullfile(arg.path.multi, char(file_to_read));
        tmp = load(fullpath);
        data_raw(k).data = [data_raw(k).data tmp.data];
        %data_raw(k).file = load(fullpath);
    end
    data_breaks = [data_breaks length(data_raw(k).data)];
end

clear tmp;

%% filter data
data = struct([]);
for k = ch
    data(k).filter = filtering(data_raw(k).data,'butter');
end

%wiener
for k = ch
    noise = var(data(k).filter(3*srCh:7*srCh));
    data(k).filter = wiener2(data(k).filter,[1 10],noise);
end
% plot un/filtered data

if (show_plot) 
figure

nb_plots = nb_files; %total number of plots

index = 1;

for k= ch
    
    %plot raw
    ax(index)=subplot(nb_plots,2,index);
    hold on;
    plot(data_raw(k).data,'b');
    ylabel(mat2str(k));
    %set(gca,'xtick',[],'ytick',[]);
    hold off;
    
    index = index +1;
    
    %plot filterd
    ax(index)=subplot(nb_plots,2,index);
    hold on;
    plot(data(k).filter,'b');
    ylabel(mat2str(k));
    %set(gca,'xtick',[],'ytick',[]);
    hold off;
 
    index = index +1;
end
        
linkaxes(ax,'xy');
clear ax;
end
%clear data_raw;

%% threshold

xmin_all = []; %all xmin
ymin_all = [];

for k = ch
    %determine threshold level
    data(k).thres = -7*median(abs(data(k).filter(:)))/0.6745;
    %find indicies when voltage return from triggered level
    returns = (find(diff(data(k).thres>data(k).filter(:))==-1)+1)';
    data(k).returns = returns;
    xmin = [];
    ymin = [];
    x1 = 1;
    for x2 = returns
        [y,x] = min(data(k).filter(x1:x2));
        ymin = [ymin y];
        xmin = [xmin (x+x1-1)];
        x1 = x2;
    end
    data(k).xmin = xmin;
    xmin_all = [xmin_all xmin];
    %data(k).ymin = ymin;
    ymin_all = [ymin_all ymin];
end

%remove multi-triggered events
[xmin_all,ix] = sort(xmin_all);
ymin_all = ymin_all(ix);

%OBS: burde nok fjerne den event med lavest y-værdi
%xmin_all = xmin_all([1 (diff(xmin_all)>8)]>0);

[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);
[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);

if (show_plot)
figure;

nb_plots = nb_files; %total number of plots

index = 1;

for k= ch
       
    %plot filterd
    ax(index)=subplot(nb_plots,1,index);
    hold on;
    plot(data(k).filter,'b');
    %threshold
    line([1 length(data(k).filter)],[data(k).thres data(k).thres]);
    %local minima
    plot(xmin_all,data(k).filter(xmin_all),'r*');
    %plot sniplet
    
    
    ylabel(mat2str(k));
    %set(gca,'xtick',[],'ytick',[]);
    hold off;
 
    index = index +1;
end
        
linkaxes(ax,'xy');

end
%% Feature extraction (PCA)

%cut out sniplet
cut_bef = 40; %samples before
cut_aft = 30; %samples after

for k= ch
    sniplet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    wl_features = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    i = 1;
    if (show_plot)
        figure
        hold on
        plot(data(k).filter);
    end
    for xs = xmin_all
        %sniplet(i,:) = data(k).filter((xs-cut_bef):(xs+cut_aft));
        %Wavelets
        %[c,l] = fix_wavedec(sniplet(i,:),4);
        sniplet(i,:) = data_raw(1,k).data((xs-cut_bef):(xs+cut_aft));
        
        %[c,l]=wavedec(sniplet(i,:),4,'haar');                    
        c = waveletcdf97(sniplet(i,:),71);
        wl_features(i,:) = c;%(1:(cut_bef+cut_aft+1));
        i = i+1;
        if (show_plot)
            plot((xs-cut_bef):(xs+cut_aft),data(k).filter((xs-cut_bef):(xs+cut_aft)),'r');
        end
    end
    i
    
    data(k).sniplet = sniplet;
    data(k).wl_features = wl_features;
    
end

%% vælg features udfra lillifors KS-test
for k=ch
    sniplet = data(k).sniplet;
    wl_features = data(k).wl_features;
    %plot sniplets
    if (show_plot) %&(k==ch(1)))
        figure
        subplot(2,1,1);
        plot(sniplet');
        subplot(2,1,2);
        temp = var(sniplet);
        plot(temp);
    end
    
    %Lillifors KS-test
    leng = size(wl_features,2);
    kstats = zeros(1,leng);
    for i = 1:leng
        %[h,p,kstat,critval] = lillietest(wl_features(:,i));
        %test_ks = wl_features(:,i);
        
        thr_dist = std(wl_features(:,i)) * 3;
            thr_dist_min = mean(wl_features(:,i)) - thr_dist;
            thr_dist_max = mean(wl_features(:,i)) + thr_dist;
            aux = wl_features(find(wl_features(:,i)>thr_dist_min & wl_features(:,i)<thr_dist_max),i);
            
            if length(aux) > 10;
                [ksstat]=test_ks(aux);
                kstats(i)=ksstat;
            else
                kstats(i)=0;
            end
        
        
        
        %kstats(i) = kstat;
    end
    
    [kstats,inx] = sort(kstats,'descend');
    
    %plot 10 best
    figure
    for i = 1:16
        subplot(8,2,i);
        hold on
        temp = wl_features(:,inx(i));
        hist(temp,30);
        set(gca,'xtick',[],'ytick',[])
        title([num2str(kstats(i)) ' , ' int2str(inx(i))])
    end
    
end

%% vælg features udfra anden test
for k=ch
    sniplet = data(k).sniplet;
    wl_features = data(k).wl_features;
    %plot sniplets
    if (show_plot) %&(k==ch(1)))
        figure
        subplot(2,1,1);
        plot(sniplet');
        subplot(2,1,2);
        temp = var(sniplet);
        plot(temp);
    end
    
    %Lillifors KS-test
    leng = size(wl_features,2);
    kstats = zeros(1,leng);

    for i = 1:leng
        %[h,p,kstat,critval] = lillietest(wl_features(:,i));
        %test_ks = wl_features(:,i);
        
        thr_dist = std(wl_features(:,i)) * 3;
            thr_dist_min = mean(wl_features(:,i)) - thr_dist;
            thr_dist_max = mean(wl_features(:,i)) + thr_dist;
            aux = wl_features(find(wl_features(:,i)>thr_dist_min & wl_features(:,i)<thr_dist_max),i);
            
            if length(aux) > 10;
                f1 = gmdistribution.fit(aux,1); % single component, unimodal
                f2 = gmdistribution.fit(aux,2); % two components, bimodal
                
                %[ksstat]=test_ks(aux);
                ksstat = f1.AIC-f2.AIC;
                %ksstat = f2.AIC;
                kstats(i)=ksstat;
            else
                kstats(i)=0;
            end
        
        
        
        %kstats(i) = kstat;
    end
    
    [kstats,inx] = sort(kstats,'descend');
    
    %plot 10 best
    if (show_plot) 
        figure
        for i = 1:16
            subplot(8,2,i);
            hold on
            temp = wl_features(:,inx(i));
            hist(temp,30);
            set(gca,'xtick',[],'ytick',[])
            title([num2str(kstats(i)) ' , ' int2str(inx(i))])
        end
    end
    
    data(k).scores = wl_features(:,inx(1:4));
end


%%
%PCA
% sniplet = [];%zeros(size(xmin_all,1)*length(ch),(cut_bef+cut_aft+1));
% sniplet = [data(ch).sniplet];
% 
%     [coefs,scores,variances,tsquare] = princomp(sniplet);
% 
%     %data(k).scores = scores(:,1:4);
%     
%     %percent_explained = (100*variances)/sum(variances);
%     
%     if (show_plot)
%     figure
%     pareto(percent_explained)
%     xlabel('Principal Component')
%     ylabel('Variance Explained (%)')
%     
%     figure
%     biplot(coefs(:,1:3), 'scores',scores(:,1:3));
%     axis([-.26 1 -.51 .51 -.61 .81]);
%     view([30 40]);
%     end

%% save data

nbDimensions = 4;

%the fet file
fid = fopen([name '.fet.1'], 'w');

% print number of feature (+ time) columns
fprintf(fid, [int2str(nbDimensions*length(ch)+1) '\n']);

fclose(fid);

scores = [];

for k = ch
    scores = [scores data(k).scores];
end
%scores = scores(:,1:16);

%[coefs,scores2,variances,tsquare] = princomp(scores);
%scores = scores2(:,1:16);

%feat_scaling_factor = 32000/max(max(max(scores)),-min(min(scores)));
feat_scaling_factor = 2147483648/max(max(max(scores)),-min(min(scores)));
scores = scores*feat_scaling_factor;
scores = round(scores);
scores = [scores xmin_all'];


    
dlmwrite([name '.fet.1'], scores, '-append', 'delimiter', ' ','precision', '%i')

%the spk file

fid = fopen([name '.spk.1'], 'w');

for k = ch
    wave_scaling_factor = 2147483648/max(max(max(data(k).sniplet)),-min(min(data(k).sniplet)));
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

fid = fopen([name '.dat'], 'w');
tmp = [];
for k = ch
    tmp =  [tmp;data(k).filter];
end
wave_scaling_factor = 2147483648/max(max(max(tmp)),-min(min(tmp)));
tmp = round(tmp*wave_scaling_factor);
fwrite(fid, tmp,'int32');
fclose(fid);

time_x1 = 0;
i = 1;
xmin_all = xmin_all';
spikes = zeros(length(xmin_all),3);
for time_x2 = data_breaks
    index = (time_x1<xmin_all)&(xmin_all<=time_x2);
    spikes(index,1) = d_nb(i);
    spikes(index,2) = xmin_all(index)-time_x1;
    time_x1 = time_x2;
    i = i+1;
end
save([name '_spikes.mat'], 'spikes')

%%

    function [xmins,ymins] = remove_doubles(xmin_all,ymin_all)
        
        old_xmin = 0;
temp_x = [];
temp_y = [];
step = 0;
for j = 1:length(xmin_all)
    if (j==1|step==1)
        old_xmin = xmin_all(j);
        step = 0;
        continue
    end
    if ((xmin_all(j)-old_xmin)>50)
        temp_x(end+1) = old_xmin;
        temp_y(end+1) = ymin_all(j-1);
        if (j==length(xmin_all))
            temp_x(end+1) = xmin_all(j);
            temp_y(end+1) = ymin_all(j);
        end
    elseif (ymin_all(j)>ymin_all(j-1))
        temp_x(end+1) = old_xmin;
        temp_y(end+1) = ymin_all(j-1);
        step = 1;
    elseif (j==length(xmin_all))
        temp_x(end+1) = xmin_all(j);
        temp_y(end+1) = ymin_all(j);
    end
    old_xmin = xmin_all(j);
end

xmins= temp_x;
ymins = temp_y;
    end
end"""

def threshold(self):
    pass

"""function [spik_sorting] = spik_threshold(spik_sorting, data, show_plot)
%SPIK_THRESHOLD find events by thresholding
%26/1/2011 Mikkel Vestergaard

%dset = spik_sorting.dset;
%d_nb = spik_sorting.d_nb;
ch = spik_sorting.ch;
thres = spik_sorting.par.thres;

cut_bef = spik_sorting.par.cut_bef; %samples before
cut_aft = spik_sorting.par.cut_aft; %samples after

xmin_all = []; %all xmin
ymin_all = [];
for k = ch
    k
    %determine threshold level
    thres_ch = thres*median(abs(spik_sorting.filtered(k).data(:)))/0.6745;
    spik_sorting.thres_ch(k).thres = thres_ch;
    %find indicies when voltage return from triggered level
    returns = (find(diff(thres_ch>spik_sorting.filtered(k).data(:))==-1)+1)';
    %data(k).returns = returns;
    xmin = [];
    ymin = [];
    x1 = 1;
    for x2 = returns
        [y,x] = min(spik_sorting.filtered(k).data(x1:x2));
        ymin = [ymin y];
        abs_x = (x+x1-1);
        %[~, inx] = min(spik_sorting.data(k).data((abs_x-10:(abs_x+10))));
        [~, inx] = min(spik_sorting.data(k).data((abs_x-3:(abs_x+3))));
        abs_x = abs_x-3+inx-1;
        %abs_x = abs_x-10+inx-1;
        xmin = [xmin abs_x];
        x1 = x2;
    end
    %data(k).xmin = xmin;
    xmin_all = [xmin_all xmin];
    %data(k).ymin = ymin;
    ymin_all = [ymin_all ymin];
end



%remove multi-triggered events
[xmin_all,ix] = sort(xmin_all);
ymin_all = ymin_all(ix);

%OBS: burde nok fjerne den event med lavest y-værdi
%xmin_all = xmin_all([1 (diff(xmin_all)>8)]>0);

[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);
[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);

spik_sorting.events_x = xmin_all;

for k= ch
    k
    snippet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    i = 1;
    %if (show_plot)
        %figure
        %hold on
        %plot(data(k).filter);
    %end
    for xs = xmin_all
        %sniplet(i,:) = data(k).filter((xs-cut_bef):(xs+cut_aft));
        %Wavelets
        %[c,l] = fix_wavedec(sniplet(i,:),4);
        snippet(i,:) = data(1,k).data((xs-cut_bef):(xs+cut_aft));
        i = i+1;
    end
    spik_sorting.snippet(k).data = snippet;
end

if (show_plot)
    figure;

    nb_plots = length(ch); %total number of plots

    index = 1;

    for k= ch

        %plot filterd
        ax(index)=subplot(nb_plots,1,index);
        hold on;
        plot(spik_sorting.filtered(k).data,'b');
        %threshold
        line([1 length(data(k).data)],[spik_sorting.thres_ch(k).thres spik_sorting.thres_ch(k).thres]);
        %local minima
        plot(xmin_all,spik_sorting.filtered(k).data(xmin_all),'r*');
        %plot sniplet


        ylabel(mat2str(k));
        %set(gca,'xtick',[],'ytick',[]);
        hold off;

        index = index +1;
    end

    linkaxes(ax,'xy');
end

function [xmins,ymins] = remove_doubles(xmin_all,ymin_all)
        
    old_xmin = 0;
    temp_x = [];
    temp_y = [];
    step = 0;
    for j = 1:length(xmin_all)
        if (j==1|step==1)
            old_xmin = xmin_all(j);
            step = 0;
            continue
        end
        if ((xmin_all(j)-old_xmin)>50)
            temp_x(end+1) = old_xmin;
            temp_y(end+1) = ymin_all(j-1);
            if (j==length(xmin_all))
                temp_x(end+1) = xmin_all(j);
                temp_y(end+1) = ymin_all(j);
            end
        elseif (ymin_all(j)>ymin_all(j-1))
            temp_x(end+1) = old_xmin;
            temp_y(end+1) = ymin_all(j-1);
            step = 1;
        elseif (j==length(xmin_all))
            temp_x(end+1) = xmin_all(j);
            temp_y(end+1) = ymin_all(j);
        end
        old_xmin = xmin_all(j);
    end

    xmins= temp_x;
    ymins = temp_y;
end

end"""

def threshold_wl(self):
    pass

"""function [spik_sorting] = spik_threshold_wl(spik_sorting, data, show_plot)
%SPIK_THRESHOLD_WL find in events by wavelets
%31/1/2011 Mikkel Vestergaard


%dset = spik_sorting.dset;
%d_nb = spik_sorting.d_nb;
ch = spik_sorting.ch;
sr_multi = spik_sorting.dset.info.sr_multi;

cut_bef = spik_sorting.par.cut_bef; %samples before
cut_aft = spik_sorting.par.cut_aft; %samples after

xmin_all = []; %all xmin


for k = ch
    k
    xmin_all = [xmin_all detect_spikes_wavelet(spik_sorting.data(k).data,sr_multi/1000,[0.5 1.0],6,'c',0.15,'bior1.5',0,1)];

end

for i = 1:length(xmin_all)
    x = xmin_all(i);
    [~, inx] = min(spik_sorting.data(14).data((x-30:(x+30))));
    xmin_all(i) = x-30+inx-1;
end

%remove multi-triggered events
[xmin_all,ix] = sort(xmin_all);
%ymin_all = ymin_all(ix);
ymin_all = zeros(size(xmin_all));
%OBS: burde nok fjerne den event med lavest y-værdi
%xmin_all = xmin_all([1 (diff(xmin_all)>8)]>0);

[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);
[xmin_all,ymin_all] = remove_doubles(xmin_all,ymin_all);

spik_sorting.events_x = xmin_all;

for k= ch
    k
    sniplet = zeros(size(xmin_all,1),(cut_bef+cut_aft+1));
    i = 1;
    %if (show_plot)
        %figure
        %hold on
        %plot(data(k).filter);
    %end
    for xs = xmin_all
        %sniplet(i,:) = data(k).filter((xs-cut_bef):(xs+cut_aft));
        %Wavelets
        %[c,l] = fix_wavedec(sniplet(i,:),4);
        sniplet(i,:) = data(1,k).data((xs-cut_bef):(xs+cut_aft));
        i = i+1;
    end
    spik_sorting.sniplet(k).data = sniplet;
end

if (show_plot)
    figure;

    nb_plots = length(ch); %total number of plots

    index = 1;

    for k= ch

        %plot filterd
        ax(index)=subplot(nb_plots,1,index);
        hold on;
        plot(data(k).data,'b');
        %threshold
        line([1 length(data(k).data)],[data(k).thres data(k).thres]);
        %local minima
        plot(xmin_all,data(k).data(xmin_all),'r*');
        %plot sniplet


        ylabel(mat2str(k));
        %set(gca,'xtick',[],'ytick',[]);
        hold off;

        index = index +1;
    end

    linkaxes(ax,'xy');
end

function [xmins,ymins] = remove_doubles(xmin_all,ymin_all)
        
    old_xmin = 0;
    temp_x = [];
    temp_y = [];
    step = 0;
    for j = 1:length(xmin_all)
        if (j==1|step==1)
            old_xmin = xmin_all(j);
            step = 0;
            continue
        end
        if ((xmin_all(j)-old_xmin)>50)
            temp_x(end+1) = old_xmin;
            temp_y(end+1) = ymin_all(j-1);
            if (j==length(xmin_all))
                temp_x(end+1) = xmin_all(j);
                temp_y(end+1) = ymin_all(j);
            end
        elseif (ymin_all(j)>ymin_all(j-1))
            temp_x(end+1) = old_xmin;
            temp_y(end+1) = ymin_all(j-1);
            step = 1;
        elseif (j==length(xmin_all))
            temp_x(end+1) = xmin_all(j);
            temp_y(end+1) = ymin_all(j);
        end
        old_xmin = xmin_all(j);
    end

    xmins= temp_x;
    ymins = temp_y;
end

end"""

def test_ks(self):
    pass
"""function [KSmax] = test_ks(x)
% 
% Calculates the CDF (expcdf)
%[y_expcdf,x_expcdf]=cdfcalc(x);

yCDF = [];
xCDF = [];
x = x(~isnan(x));
n = length(x);
x = sort(x(:));
% Get cumulative sums
yCDF = (1:n)' / n;
% Remove duplicates; only need final one with total count
notdup = ([diff(x(:)); 1] > 0);
x_expcdf = x(notdup);
y_expcdf = [0; yCDF(notdup)];

%
% The theoretical CDF (theocdf) is assumed to be normal  
% with unknown mean and sigma

zScores  =  (x_expcdf - mean(x))./std(x);

%theocdf  =  normcdf(zScores , 0 , 1);
mu = 0; 
sigma = 1; 
theocdf = 0.5 * erfc(-(zScores-mu)./(sqrt(2)*sigma));


%
% Compute the Maximum distance: max|S(x) - theocdf(x)|.
%

delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
deltacdf  =  abs([delta1 ; delta2]);

KSmax =  max(deltacdf);"""

def waveletcdf97(self):
    pass

"""function X = waveletcdf97(X, Level)
%WAVELETCDF97  Cohen-Daubechies-Feauveau 9/7 wavelet transform.
%   Y = WAVELETCDF97(X, L) decomposes X with L stages of the
%   Cohen-Daubechies-Feauveau (CDF) 9/7 wavelet.  For the
%   inverse transform, WAVELETCDF97(X, -L) inverts L stages.
%   Filter boundary handling is half-sample symmetric.
%
%   X may be of any size; it need not have size divisible by 2^L.
%   For example, if X has length 9, one stage of decomposition
%   produces a lowpass subband of length 5 and a highpass subband
%   of length 4.  Transforms of any length have perfect
%   reconstruction (exact inversion).
%
%   If X is a matrix, WAVELETCDF97 performs a (tensor) 2D wavelet
%   transform.  If X has three dimensions, the 2D transform is
%   applied along the first two dimensions.
%
%   Example:
%   Y = waveletcdf97(X, 5);    % Transform image X using 5 stages
%   R = waveletcdf97(Y, -5);   % Reconstruct from Y

% Pascal Getreuer 2004-2006

if nargin < 2, error('Not enough input arguments.'); end
if ndims(X) > 3, error('Input must be a 2D or 3D array.'); end
if any(size(Level) ~= 1), error('Invalid transform level.'); end

N1 = size(X,1);
N2 = size(X,2);

% Lifting scheme filter coefficients for CDF 9/7
LiftFilter = [-1.5861343420693648,-0.0529801185718856,0.8829110755411875,0.4435068520511142];
ScaleFactor = 1.1496043988602418;

S1 = LiftFilter(1);
S2 = LiftFilter(2);
S3 = LiftFilter(3);
ExtrapolateOdd = -2*[S1*S2*S3,S2*S3,S1+S3+3*S1*S2*S3]/(1+2*S2*S3);

LiftFilter = LiftFilter([1,1],:);

if Level >= 0   % Forward transform
   for k = 1:Level
      M1 = ceil(N1/2);
      M2 = ceil(N2/2);
      
      %%% Transform along columns %%%
      if N1 > 1         
         RightShift = [2:M1,M1];
         X0 = X(1:2:N1,1:N2,:);

         % Apply lifting stages
         if rem(N1,2)
            X1 = [X(2:2:N1,1:N2,:);X0(M1-1,:,:)*ExtrapolateOdd(1)...
                  + X(N1-1,1:N2,:)*ExtrapolateOdd(2)...
                  + X0(M1,:,:)*ExtrapolateOdd(3)]...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         else
            X1 = X(2:2:N1,1:N2,:) ...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         end

         X0 = X0 + filter(LiftFilter(:,2),1,...
            X1,X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 + filter(LiftFilter(:,3),1,...
            X0(RightShift,:,:),X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 + filter(LiftFilter(:,4),1,...
            X1,X1(1,:,:)*LiftFilter(1,4),1);

         if rem(N1,2)
            X1(M1,:,:) = [];
         end

         X(1:N1,1:N2,:) = [X0*ScaleFactor;X1/ScaleFactor];
      end

      %%% Transform along rows %%%
      if N2 > 1
         RightShift = [2:M2,M2];
         X0 = permute(X(1:N1,1:2:N2,:),[2,1,3]);

         % Apply lifting stages
         if rem(N2,2)
            X1 = permute([X(1:N1,2:2:N2,:),X(1:N1,N2-2,:)*ExtrapolateOdd(1)...
                  + X(1:N1,N2-1,:)*ExtrapolateOdd(2) ...
                  + X(1:N1,N2,:)*ExtrapolateOdd(3)],[2,1,3])...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         else
            X1 = permute(X(1:N1,2:2:N2,:),[2,1,3]) ...
               + filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
               X0(1,:,:)*LiftFilter(1,1),1);
         end

         X0 = X0 + filter(LiftFilter(:,2),1,...
            X1,X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 + filter(LiftFilter(:,3),1,...
            X0(RightShift,:,:),X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 + filter(LiftFilter(:,4),1,...
            X1,X1(1,:,:)*LiftFilter(1,4),1);

         if rem(N2,2)
            X1(M2,:,:) = [];
         end

         X(1:N1,1:N2,:) = permute([X0*ScaleFactor;X1/ScaleFactor],[2,1,3]);
      end

      N1 = M1;
      N2 = M2;
   end
else           % Inverse transform
   for k = 1+Level:0
      M1 = ceil(N1*pow2(k));
      M2 = ceil(N2*pow2(k));

      %%% Inverse transform along rows %%%
      if M2 > 1
         Q = ceil(M2/2);
         RightShift = [2:Q,Q];
         X1 = permute(X(1:M1,Q+1:M2,:)*ScaleFactor,[2,1,3]);

         if rem(M2,2)
            X1(Q,1,1) = 0;
         end

         % Undo lifting stages
         X0 = permute(X(1:M1,1:Q,:)/ScaleFactor,[2,1,3]) ...
            - filter(LiftFilter(:,4),1,X1,X1(1,:,:)*LiftFilter(1,4),1);
         X1 = X1 - filter(LiftFilter(:,3),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 - filter(LiftFilter(:,2),1,X1,...
            X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 - filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,1),1);

         if rem(M2,2)
            X1(Q,:,:) = [];
         end

         X(1:M1,[1:2:M2,2:2:M2],:) = permute([X0;X1],[2,1,3]);
      end

      %%% Inverse transform along columns %%%
      if M1 > 1
         Q = ceil(M1/2);
         RightShift = [2:Q,Q];
         X1 = X(Q+1:M1,1:M2,:)*ScaleFactor;

         if rem(M1,2)
            X1(Q,1,1) = 0;
         end

         % Undo lifting stages
         X0 = X(1:Q,1:M2,:)/ScaleFactor ...
            - filter(LiftFilter(:,4),1,X1,X1(1,:,:)*LiftFilter(1,4),1);
         X1 = X1 - filter(LiftFilter(:,3),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,3),1);
         X0 = X0 - filter(LiftFilter(:,2),1,X1,...
            X1(1,:,:)*LiftFilter(1,2),1);
         X1 = X1 - filter(LiftFilter(:,1),1,X0(RightShift,:,:),...
            X0(1,:,:)*LiftFilter(1,1),1);

         if rem(M1,2)
            X1(Q,:,:) = [];
         end

         X([1:2:M1,2:2:M1],1:M2,:) = [X0;X1];
      end
   end
end"""