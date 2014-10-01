##########################
# SpikeDetekt parameters #
##########################

experiment_name = 'test_hybrid_120sec'
raw_data_files = experiment_name + '.raw.kwd'
prb_file = '32chan1shankbuzsaki.prb'
nbits = 16
voltage_gain = 10.

sample_rate = 20000
nchannels = 32


# Filtering
# ---------
filter_low = 500. # Low pass frequency (Hz)
filter_high = 0.95 * .5 * sample_rate
filter_butter_order = 3  # Order of Butterworth filter.

filter_lfp_low = 0  # LFP filter low-pass frequency
filter_lfp_high = 300  # LFP filter high-pass frequency


# Chunks
# ------
chunk_size = int(1. * sample_rate)  # 1 second
chunk_overlap = int(.015 * sample_rate)  # 15 ms

# Spike detection
# ---------------
# Uniformly scattered chunks, for computing the threshold from the std of the
# signal across the whole recording.
nexcerpts = 50
excerpt_size = int(1. * sample_rate)
threshold_strong_std_factor = 7.5
threshold_weak_std_factor = 5.
detect_spikes = 'negative'
#precomputed_threshold = None

# Connected component
# -------------------
connected_component_join_size = int(.00005 * sample_rate)
    
# Spike extraction
# ----------------
extract_s_before = 16
extract_s_after = 16
waveforms_nsamples = extract_s_before + extract_s_after

# Features
# --------
nfeatures_per_channel = 3  # Number of features per channel.
pca_nwaveforms_max = 10000


#########################
# KlustaKwik parameters #
#########################
MaskStarts = 100
#MinClusters = 100 
#MaxClusters = 110
MaxPossibleClusters =  500
FullStepEvery =  10
MaxIter = 10000
RandomSeed =  654
Debug = 0
SplitFirst = 20 
SplitEvery = 100 
PenaltyK = 0
PenaltyKLogN = 1
Subset = 1
PriorPoint = 1
SaveSorted = 0
SaveCovarianceMeans = 0
UseMaskedInitialConditions = 1 
AssignToFirstClosestMask = 1
UseDistributional = 1
