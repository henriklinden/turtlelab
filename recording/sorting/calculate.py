from __future__ import division
from scipy import signal
import numpy as np
from pylab import rms_flat

def evaluate_cluster(self):
        pass
    
def harris_isolation_distance(self):
    pass
    """
    function output = spik_eval_harris_isola_dist(spik_sorting,cluster_nb,show)
%% SPIK_EVAL_HARRIS_DIST Evaluate the isolation distance (K Harris)
% 2011 Mikkel Vestergaard

features = spik_sorting.scores.data;
clusters = spik_sorting.clusters;
output = ClusterQuality(features,(find(clusters==cluster_nb)));

end


function [HalfDist] = ClusterQuality(Fet, MySpikes)
%functionen er faaet fra Pablo


% check there are enough spikes (but not more than half)
if length(MySpikes) < size(Fet,2) | length(MySpikes)>size(Fet,1)/2
    HalfDist = 0;
    return
end

% find spikes in this cluster
nSpikes = size(Fet,1);
nMySpikes = length(MySpikes);
InClu = ismember(1:nSpikes, MySpikes);


% mark other spikes
OtherSpikes = setdiff(1:nSpikes, MySpikes);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
m = mahal(Fet, Fet(MySpikes,:));
mMy = m(MySpikes); % mahal dist of my spikes
mOther = m(OtherSpikes); % mahal dist of others

%fprintf('done mahalanobis calculation...');
% calculate point where mD of other spikes = n of this cell
if (nMySpikes < nSpikes/2)
    [sorted order] = sort(mOther);
    HalfDist = sorted(nMySpikes);
else
    HalfDist = 0; % If there are more of this cell than every thing else, forget it.
end
    
%%%%%%%%%%%%%% plotting    %%%%%%%%%%%%%%%%%%%%%%%%%%%
mmin = min(m);
mmax = max(m);
        

% make histograms
%Bins = 0:1:300;
%hMy = hist(mMy, Bins);
%hOther = hist(mOther, Bins);

%clf

% % plot pdfs
% subplot(2,1,1)
% 
% semilogy(Bins, [hMy', hOther']);
% xlim([0 100]);
% ylabel('spike density');
% xlabel('Mahalanobis distance');
% legend('This cluster', 'Others');
% 
% % plot cdfs
% subplot(2,1,2)
% cdMy = cumsum(hMy);% / length(MySpikes);
% cdOther = cumsum(hOther);% / length(OtherSpikes);
% semilogy(Bins, [cdMy', cdOther']);
% ylabel('cumulative # spikes');
% xlabel('Mahalanobis distance');
%xlim([0 100]);

end"""
    
    def hill_censored(self):
        pass
    
    def hill_overlap(self):
        pass
    """
    function output = spik_eval_hill_overlap(spik_sorting,cluster_nb,show)
%% SPIK_EVAL_HILL_NEG_THRES calculate overlap of clusters (Hill et al. 2011)


features = spik_sorting.scores.data;
clusters = spik_sorting.clusters;

%%


nb_dim = length(models(1).center);
features2 = features(:,1:nb_dim);
for i = 1:size(features2,2)
    mmin = min(features2(:,i));
    mmax = max(features2(:,i));
    features2(:,i) = (features2(:,i)-mmin)/(mmax-mmin);
end


prob_c_f = zeros(size(features2,1),length(models));
prob_f_c = prob_c_f;
prob_c = zeros(1,length(models));
for j = 1:length(models)
    prob_c(j) = models(j).weight;
end

for i = 1:size(features2,1)
    fea = features2(i,:);
    for j = 1:length(models)
        if (j==1)
            prob_f_c(i,j) = (models(1).weight);
        else
            model = models(j);
            prob_f_c(i,j) = spik_prob_fea_clu(fea,model);
        end
    end
end

for i = 1:size(features2,1)
    prob_f = 0;
    for k = 1:length(models)
        prob_f = prob_f +prob_f_c(i,k)*prob_c(k);
    end
    for j = 1:length(models)
        if (j==1)
            prob_c_f(i,j) = prob_f_c(i,j)/prob_f;
        else
            prob_c_f(i,j) = prob_f_c(i,j)*prob_c(j)/prob_f;
        end
    end
end
[~,clusters2] = max(prob_c_f,[],2);
%% Cal false_pos and false_neg
nb = 6;
false_pos = 0;

temp_idx = find(clusters==nb);
for i = 1:max(clusters)
    if (i~=nb) 
        false_pos = false_pos + sum(prob_c_f(temp_idx,i))/length(temp_idx);
    end
end

temp_idx2 = find(clusters~=nb);
false_neg =sum(prob_c_f(temp_idx2,nb))/length(temp_idx);

%idx = find(clusters==nb);
false_pos
false_neg"""
    
    def hill_refractory(self):
        pass
    """
    function output = spik_eval_hill_refrac(spik_sorting,cluster_nb,show)
%% SPIK_EVAL_HILL_REFRAC # of false positive due to refractory violations (Hill et al. 2011)
% 2011 Mikkel Vestergaard


clusters = spik_sorting.clusters;
which_spikes = find(clusters==cluster_nb);

refrac = spik_sorting.par.refrac; %refractory period
remove_sp = spik_sorting.par.remove_sp; %time where two spikes are not possible due to the analyses
spikes = spik_sorting.events_x(which_spikes)./(spik_sorting.dset.info.sr_multi/1000);

isis = diff(spikes);   

refrac_vio = sum(isis<refrac); %number of refractory period violations

T = sum(isis(isis <= 5000)); %remove long isi and find total isi

[ev,lb,ub] = rpv_contamination(length(spikes), T, refrac-remove_sp, refrac_vio );
output = [ev, lb,ub];"""
    
    def hill_neg_thres(self):
        pass
    """
    function output = spik_eval_hill_neg_thres(spik_sorting,cluster_nb,show)
%% SPIK_EVAL_HILL_NEG_THRES # of false neg due to too high threshold (Hill et al. 2011)
% 2011 Mikkel Vestergaard


clusters = spik_sorting.clusters;
which_spikes = find(clusters==cluster_nb);

snippets = [];

for k = spik_sorting.ch
    snippet = spik_sorting.snippet(k).data(which_spikes,:);
    snippets = [snippets min(snippet,[],2)];
end

[n,xout] = hist(min(snippets,[],2),100);
figure
bar(xout,n);

output = NaN;"""
    
    def redish_lratio(self):
        pass
    
    def tolias_isolation_distance(self):
        pass
    
    def detect_spikes_wavelet(self):
        pass
    """
    function TE = detect_spikes_wavelet(...
    Signal, SFr, Wid, Ns, option, L, wname, PltFlg, CmtFlg)

% DETECT_SPIKES_WAVELET wavelet based algorithm for detection of transients
% from neural data.
%
%   TE=DETECT_SPIKES_WAVELET(Signal,SFr,Wid,Ns,option,L,wname,PltFlg,CmtFlg)
%
%   Signal - extracellular potential data to be analyzed 1 x Nt;
%
%   SFr - sampling frequency [kHz];
%
%   Wid - 1 x 2 vector of expected minimum and maximum width [msec] of transient 
%   to be detected Wid=[Wmin Wmax]. For most practical purposes Wid=[0.5 1.0];
%
%   Ns - (scalar): the number of scales to use in detection (Ns >= 2);
%
%   option - (string): the action taken when no coefficients survive hard 
%   thresholding 
%   'c' means conservative and returns no spikes if P(S) is found to be 0
%   'l' means assume P(S) as a vague prior (see the original reference)
%
%   L is the factor that multiplies [cost of comission]/[cost of omission].
%   For most practical purposes -0.2 <= L <= 0.2. Larger L --> omissions
%   likely, smaller L --> false positives likely. For unsupervised
%   detection, the suggested value of L is close to 0.  
%
%   wname - (string): the name of wavelet family in use
%           'bior1.5' - biorthogonal
%           'bior1.3' - biorthogonal
%           'db2'     - Daubechies
%           'sym2'    - symmlet
%           'haar'    - Haar function
%   Note: sym2 and db2 differ only by sign --> they produce the same
%   result;
%
%   PltFlg - (integer) is the plot flag: 
%   PltFlg = 1 --> generate figures, otherwise do not;
%  
%   CmtFlg - (integer) is the comment flag, 
%   CmtFlg = 1 --> display comments, otherwise do not;
%
%   TE is the vector of event occurrence times;
%
%   Reference: Z. Nenadic and J.W. Burdick, Spike detection using the 
%   continuous wavelet transform, IEEE T. Bio-med. Eng., vol. 52,
%   pp. 74-87, 2005.

%   Originally developed by:
%   Zoran Nenadic
%   California Institute of Technology
%   May 2003
%
%   Modified by:
%   Zoran Nenadic
%   University fo California, Irvine
%   February 2008


%admissible wavelet families (more wavelets could be added)
wfam = {'bior1.5','bior1.3','sym2','db2','haar'};

if sum(strcmp(wname,wfam)) == 0
    error('unknown wavelet family')
elseif CmtFlg == 1
    disp(['wavelet family: ' wname])
    to = clock;
end

%make sure signal is zero-mean
Signal = Signal - mean(Signal);

Nt = length(Signal);      %# of time samples

%define relevant scales for detection
W = determine_scales(wname,Wid,SFr,Ns);

%initialize the matrix of thresholded coefficients
ct = zeros(Ns,Nt);

%get all coefficients 
c = cwt(Signal,W,wname);  

%define detection parameter
Lmax = 36.7368;       %log(Lcom/Lom), where the ratio is the maximum 
                    %allowed by the current machine precision
L = L * Lmax;

%initialize the vector of spike indicators, 0-no spike, 1-spike
Io = zeros(1,Nt);

%loop over scales
for i = 1:Ns
    
    %take only coefficients that are independent (W(i) apart) for median
    %standard deviation
    
    Sigmaj = median(abs(c(i,1:round(W(i)):end) - mean(c(i,:))))/0.6745;
    Thj = Sigmaj * sqrt(2 * log(Nt));     %hard threshold
    index = find(abs(c(i,:)) > Thj);
    if isempty(index) & strcmp(num2str(option),'c')
        %do nothing ct=[0];
    elseif isempty(index) & strcmp(num2str(option),'l')
        Mj = Thj;
        %assume at least one spike
        PS = 1/Nt;
        PN = 1 - PS;
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];    %decision threshold
        DTh = abs(DTh) * (DTh >= 0);                 %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        if isempty(ind)
            %do nothing ct=[0];
        else
            ct(i,ind) = c(i,ind);
        end
    else
        Mj = mean(abs(c(i,index)));       %mean of the signal coefficients
        PS = length(index)/Nt;            %prior of spikes
        PN = 1 - PS;                        %prior of noise
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];   %decision threshold
        DTh = abs(DTh) * (DTh >= 0);         %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        ct(i,ind) = c(i,ind);
    end
    
    %find which coefficients are non-zero
    Index = ct(i,:) ~= 0;
    
    %make a union with coefficients from previous scales
    Index = or(Io,Index);
    Io = Index;
end

TE = parse(Index,SFr,Wid);

if PltFlg == 1
    close all
    figure(1)
    scale = 64./[max(abs(c),[],2) * ones(1,Nt)];
    temp = zeros(1,Nt);
    temp(TE) = 1;
    image(flipud(abs(c)) .* scale)
    colormap pink
    ylabel('Scales')
    Wt = [fliplr(W)];
    set(gca,'YTick',1:Ns,'YTickLabel',Wt,'Position',[0.1 0.2 0.8 0.6], ...
        'XTick',[])
    title(['|C| across scales: ' num2str(W)])
    ah2 = axes;
    set(ah2,'Position',[0.1 0.1 0.8 0.1])
    plot(temp,'o-m','MarkerSize',4,'MarkerFaceColor','m')
    set(gca,'YTick',[],'XLim',[1 Nt])
    xlabel('Time (samples)')
    ylabel('Spikes')
    
    figure(2)
    plot(Signal,'Color',[0.7 0.7 0.7],'LineWidth',2)
    hold on
    plot(ct','-o','LineWidth',1,'MarkerFaceColor','k', ...
        'MarkerSize',4)
    xlabel('Time (samples)')
    ylabel('Coefficients')
    set(gca,'XLim',[1 Nt])
end

if CmtFlg == 1
    disp([num2str(length(TE)) ' spikes found'])
    disp(['elapsed time: ' num2str(etime(clock,to))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Scale = determine_scales(wname,Wid,SFr,Ns)
  
%Ns - # of scales  

dt = 1/SFr;  %[msec]

%signal sampled @ 1 KHz  
Signal = zeros(1,1000);
%create Dirac function
Signal(500) = 1;
  
Width = linspace(Wid(1),Wid(2),Ns);

%infinitesimally small number
Eps = 10^(-15);

ScaleMax = 3;
ScaleMax = ScaleMax*SFr;

switch num2str(wname)
  
 case 'haar'
  for i = 1:Ns
    Scale(i) = Width(i)/dt - 1; 
  end
 case 'db2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'sym2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of positive slope zero crossings
    IndZeroCross = find(IndDer == 1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.3'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
   for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.5'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 otherwise
  error('unknown wavelet family')
end

NaNInd = isnan(Scale);

if sum(NaNInd) > 0
  warning(['Your choice of Wid is not valid given' ...
        ' the sampling rate and wavelet family'])
  if NaNInd(1) == 1
    disp(['Most likely Wid(1) is too small'])
  elseif NaNInd(Ns) == 1
    disp(['Most likely Wid(2) is too large'])
    disp(['Change the value on line: ''ScaleMax = 2'' to something larger'])
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn = parse(Index,SFr,Wid);

%This is a special function, it takes the vector Index which has 
%the structure [0 0 0 1 1 1 0 ... 0 1 0 ... 0]. This vector was obtained
%by coincidence detection of certain events (lower and upper threshold
%crossing for threshold detection, and the appearance of coefficients at
%different scales for wavelet detection). 
%The real challenge here is to merge multiple 1's that belong to the same
%spike into one event and to locate that event

Refract = 1.5 * Wid(2);    %[ms] the refractory period -- can't resolve spikes 
                           %that are closer than Refract;
Refract = round(Refract * SFr);

Merge = mean(Wid);      %[ms] merge spikes that are closer than Merge, since 
                        %it is likely they belong to the same spike

Merge = round(Merge * SFr);   


Index([1 end]) = 0;   %discard spikes located at the first and last samples

ind_ones = find(Index == 1);    %find where the ones are

if isempty(ind_ones)
    TE = [];
else
    temp = diff(Index);  %there will be 1 followed by -1 for each spike
    N_sp = sum(temp == 1); %nominal number of spikes
    
    lead_t = find(temp == 1);  %index of the beginning of a spike
    lag_t = find(temp == -1);  %index of the end of the spike
    
    for i = 1:N_sp
        tE(i) = ceil(mean([lead_t(i) lag_t(i)]));
    end
   
    i = 1;        %initialize counter
    while 0 < 1
        if i > (length(tE) - 1)
            break;
        else
            Diff = tE(i+1) - tE(i);
            if Diff < Refract & Diff > Merge
                tE(i+1) = [];      %discard spike too close to its predecessor
            elseif Diff <= Merge
                tE(i) = ceil(mean([tE(i) tE(i+1)]));  %merge
                tE(i+1) = [];                         %discard
            else
                i = i+1;
            end
        end
    end 
    TE = tE;
end

fcn = TE;"""

    def pxcorr(self):
        pass
    """
    function [C, lags] = pxcorr(x, varargin)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% pxcorr -  Efficient cross-correlation for point process data
%
%
% Description:
%   C = PXCORR(A,B,Fs) returns the cross-correlation between A and B at
%   sampling rate Fs for the point process whose times are given in the
%   length N(>1) vector A and the length M (>1) vector B.  C will be a row
%   vector with length given by
%                2*Fs*max(max(A)-min(B),max(B)-min(A))+1
%   i.e., twice the maximum time difference between events in A and B
%   measured in units of 1/Fs (plus 1 for 0 lag).  The vector C estimates
%                       C(t) = E[A(x-t)*B(x)]
%   The event times listed in A and B are assumed to be sorted and will be
%   rounded to bins of duration 1/Fs. 
%
%   The computation here is O(N*J), where 0<=J<=M is the mean number of
%   events in the vector B that are within MAXLAG (see below) of events in
%   the vector A (when MAXLAG is not specified, J=M).  Matlab's XCORR is
%   O(K log K), where K = max(max(A),max(B))*Fs -- note that if MAXLAG is
%   large and the mean interevent interval is small, XCORR can be faster.
%
%   PXCORR(A,Fs) returns the auto-correlation of the events listed in A 
%   and is equivalent to PXCORR(A,A,Fs)
%
%   PXCORR(...,MAXLAG) only considers lags upto MAXLAG; the length of C
%   will then be 2*Fs*MAXLAG + 1.
%
%   PXCORR(...,'sort') indicates the inputs are not yet sorted.
%
%   [C,LAGS] = PXCORR(...) also returns a vector of lag indices.

%%%%%%%%%%%%%%%%%%%%%%%%%% Argument Parsing %%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<2), error('Not enough arguments.');  end;
if (length(x)<=1), error('Point processes must consist of >1 event.'); end;
if (~isvector(x)), error('Matrix data are not supported.');  end;

presorted = 1;   % already sorted events?
if (strcmpi(varargin{end},'sort')),
    presorted = 0;
    varargin = varargin(1:end-1);   
    if(length(varargin)==0), error('Not enough arguments.'); end;
end;

if (length(varargin{1}) == 1), y = x;    % make auto-corr case look like X
else,
    y = varargin{1};
    varargin = varargin(2:end);
    if(length(varargin)==0), error('Not enough arguments.'); end;
end
if (length(y)<=1),  error('Point processes must consist of >1 event.');  end;
if (~isvector(x)), error('Matrix data are not supported.');  end;

x = x(:)';  y = y(:)';    % enforce row vectors
if (~presorted),  x = sort(x);  y = sort(y);  end;   % only do this if needed

nargs = length(varargin);
Fs = varargin{1};           
if(length(Fs)~=1), error('Fs must be a scalar.'); end;

if (nargs == 1)   % limiting the lag saves time
    maxlag = max(x(end)-y(1),y(end)-x(1));
elseif (nargs == 2)
    maxlag = varargin{2};           
    if(length(maxlag)~=1 || isinf(maxlag)), error('MAXLAG must be a finite scalar.'); end;
elseif ((nargs == 3) && (length(varargin{3})<=1))
    error('Point processes must consist of > 1 event.');
else
    error('Invalid syntax.');
end

x = round(x*Fs);  y = round(y*Fs);
maxlag = ceil(maxlag * Fs);

if (nargout == 2),  lags = [-maxlag:maxlag]./Fs;  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = zeros(1,2*maxlag+1);

% We do this efficiently by stepping along X and keeping track of those
% indices in Y that are near by (i.e., within a distance of MAXLAG).
limit = length(x);
a = 1;  c = 1;   

for b = 1:length(y)
    while((y(b)-x(a)) > maxlag),        % move left bracket until less than MAXLAG before x(b)
            a=a+1;   if (a > limit), return; end;
    end
    if (c < a), c = a; end;             % catch up the right bracket if it was passed by
    if (c <= limit)
        while((x(c)-y(b)) <= maxlag),   % move right bracket until more than MAXLAG after x(b)
            c=c+1;   if (c > limit), break; end;
        end
    end

    offset = -y(b)+maxlag+1;            % add 'em up
    for bb = a:(c-1)
        ind = x(bb)+offset;
        C(ind) = C(ind) + 1;
    end
end

% normalize to Hz
C = C *Fs / length(x);"""

def episode_events(self):
    pass

"""
function output =  spik_episode_events(spik_sorting)
%SPIK_EPISODE_EVENTS find events from a specific episode
clusters = spik_sorting.clusters;

epis = struct;
for i = spik_sorting.d_nb
    idx_d = find(spik_sorting.d_nb==i);

    if (idx_d==1)
        t1 = 1;
    else
       % t1 = sum(spik_sorting.data_breaks(1:idx_d-1))+1;
        t1 = spik_sorting.data_breaks(idx_d-1);
    end
    
    %t2 = sum(spik_sorting.data_breaks(1:idx_d));
    t2 = spik_sorting.data_breaks(idx_d);
    idx = (spik_sorting.events_x>=t1)&(spik_sorting.events_x<=t2);
    
    epis(i).events_x = spik_sorting.events_x(idx)-t1+1;
    epis(i).clusters = clusters(idx);
end

output = epis;

end
    """
    
def evaluate_cluster(self):
    pass

"""function output = spik_evaluate_cluster(spik_sorting,cluster_nb,show)
%SPIK_EVALUATE_CLUSTER evaluate quantatively the quality of a cluster
%   2011 Mikkel Vestergaard


%% refractory period violations

output.hill_refrac = spik_eval_hill_refrac(spik_sorting,cluster_nb,show);


%% undetected spikes

output.hill_thres = spik_eval_hill_thres(spik_sorting,cluster_nb,show);

%% overlap

%output.hill_overlap = spik_eval_hill_overlap(spik_sorting,cluster_nb,show);

%% censored

%output.hill_censored = spik_eval_hill_censored(spik_sorting,cluster_nb,show);

%% isolation distance Harris

output.harris_isola = spik_eval_harris_isola_dist(spik_sorting,cluster_nb,show);

%% Lratio Redish

%output.redish_lratio = spik_eval_redish_lratio(spik_sorting,cluster_nb,show);

%% Isolation quality Tolias

%output.tolias_isola = spik_eval_tolias_isola(spik_sorting,cluster_nb,show);
        
end

"""

def plot_corr(self):
    pass

"""
function spik_plot_corr(spikesA, spikesB)
%SPIK_PLOT_CORR plot auto/cross-correlations
%   2011 Mikkel Vestergaard

%% Autocorrelation
maxlag = 400;% data.autocorr_maxlag;
        

%[c, lags] = pxcorr(spikes,spikes,10, maxlag);
%c(find(lags==0)) = 0;
figure
%bar(lags, c,1)

% calculate autocorrelation
if length(spikesA) > 1
  %[cross,lags] = pxcorr(spikes,spikes, round(1000/data.corr_bin_size), maxlag);
  [cross,lags] = pxcorr(spikesA,spikesB, 20, maxlag);
else
    cross = 0;  lags = 0;
end
cross(find(lags==0)) = 0;


remove_sp = 30;
%ymax = max(n);
refrac = 10;
        
%indicate period where two spikes not possible due to the analyses
%patch([0 remove_sp remove_sp 0 ], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none'); %gray box
%indicate refractory period
%patch([remove_sp [refrac refrac] remove_sp ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none'); %red box



% place patches to represent shadow and refractory period
ymax = max(cross);
patch(remove_sp*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
patch([remove_sp [refrac refrac] remove_sp ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
patch(-[remove_sp [refrac refrac] remove_sp ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');

% plot autocorrelation histogram
hold on, bar(lags,cross,1.0); hold off;  
%set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

% set axes
%set(gca, 'XLim', maxlag*1000*[-1 1]);
%set(gca,'YLim',[0 ymax])
%xlabel( 'Time lag (msec)')
%ylabel({'Autocorrelation (Hz)',data.ystr})
        
end
"""

def prob_fea_clu(self):
    pass

"""function output = spik_prob_fea_clu(fea,model)
    center = model.center;
    sig =  triu(model.sigma)+triu(model.sigma,1)';
    %sig = model.sigma;
    %sig = model.sigma*model.sigma';
    weight = model.weight;
    nb_dim = length(center);
    fea = fea';
    det_sig = prod(diag(model.sigma));
    %output = (2*pi)^(-nb_dim/2)*det(sig)^(-0.5)*exp(-0.5*sum((model.sigma)/(fea-center)').^2);
    output = (2*pi)^(-nb_dim/2)*det(sig)^(-0.5)*exp(-0.5*(fea-center)'/sig*(fea-center));
    %output = (2*pi)^(-nb_dim/2)*det_sig^(-1)*exp(-0.5*sum(trisolve(model.sigma,(fea-center)).^2));

    
    [test_cho, ~] = lu(sig);
    Root = trisolve(test_cho,(fea-center));
   %Root = test_cho\(fea-center); 
   Mahal = 0;
    for i = 1:nb_dim
        Mahal = Mahal + Root(i)*Root(i);
    end
    
    temp = diag(model.sigma);
    LogRootDet = 0;
    for i = 1:nb_dim
        LogRootDet = LogRootDet + log(temp(i));
    end
    %output = Mahal/2 + LogRootDet- log(weight)+log(2*pi)*nb_dim/2;
   % output = ((fea-center)'/sig*(fea-center))/2 - log(det(sig)^(-0.5))- log(weight)+log(2*pi)*nb_dim/2;
   % output = exp(-output);
    %output = ((fea-center)'/(sig)*(fea-center))/2;
end

function output = trisolve(M, x)

D = length(x);
output = zeros(size(x));
sum = 0;
for i=1:D
    sum = x(i);
    j=i-1;
    while j >0 
        sum = sum - M(i,j)*output(j);
        j = j-1;
    end
    
    output(i) = sum/M(i,i);
end
end"""