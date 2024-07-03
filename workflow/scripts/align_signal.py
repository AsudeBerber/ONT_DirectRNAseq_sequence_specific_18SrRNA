__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

"""
this is a module called by both signal_summary.py and seq2mv_direct_RNA.py

for both single and multiple reads, the movetable is accessed and split into its parts, enabling the alignment of the signal to the query sequence/position in the reference
for the signal summary, the signal features per base are calculated and stored in a dictionary

to speed things up, only the features/signal locations for the needed base positions (window around given positions, e.g. 1842 +- 8bp) are returned
"""
import numpy as np

# for getting position in signal (goes 3' to 5') = position in sequence when starting from 3' end
def rev_locus(locus, read):
    locus_rev = read.query_length -1 - locus
    return locus_rev

# gives positions in query sequence (as well as position in reversed sequence (3'-> 5')) for given positions in reference after filtering out unusable positions (if read to short)
def get_loci(read, pairs, wd, motif_length, ref_pos):

    # reference locus, index gives index for reference pos. in alignment pairs (starts from first base in query, e.g. index = 0, ref_loci = 200 for (query=0, ref=200, seq="A"))
    ref_loci = []
    ref_loci_index = []
    # reversed loci, for operations from 3' end
    rev_loci = []

    #ref_pos: list of positions to look at (single read plot: only one position; signal_summary: i.e. all acetylated positions)
    #pairs[0]: query pos; [1]: ref pos; [2] ref base
    for i, pos in enumerate(ref_pos):
        # for empty array (= read doesn't cover requested positions)
        if pairs.shape == (0,):
            continue
        else:
            if (pos in pairs[:,1]): ref_loci.append(pos)
     
        # in cases where a read neither aligns to any of the ref positions, ref_loci_index would append [], which raises an index error --> solved by checking if array is empty
        index_pos = np.where(pairs[:,1] == pos)[0]
        if index_pos.shape == (0,):
            continue
        else:
            ref_loci_index.append(index_pos[0])
            
    # gets query pos for given pos in reference (using a dictionary is not working, as query/ref pos have to be found by using the other one)
    loci = [pairs[locus, 0] for locus in ref_loci_index]
    # Remove loci that are not present on the query or too close to the ends of the alignment, matches ref_loci, so that is wont include positions to removed loci
    loci = [locus for locus in loci if locus is not None and locus > wd and locus < read.query_length - wd - (motif_length)]
    ref_loci = [pairs[index, 1] for index in ref_loci_index if pairs[index, 0] != None and pairs[index,0] in loci]
    if len(loci) != len(ref_loci):
        raise Exception ("length of reference and query sequence index not matching")
    
    # rev_loci has reverse (3'-> 5') positions to needed query positions
    for locus in loci:
        rev_loci.append(rev_locus(locus, read))
    return loci, ref_loci, rev_loci

# this part is modified from https://github.com/WGLab/DeepMod2/blob/main/src/detect.py
# reads movetable and calculates features per read (for signal_summary) or build array with signal positions inferred from movetable (for seq2mv*.py)
def access_mv(signal, moves, offset, rev_loci, motif_length, extra_window, read, mode):
    """
    Normalises and collapses the signal based on the moves table. Outputs an array with the
    following values for each called based:
    4: log10 signal length
    5: mean signal intensity
    6: standard deviation of signal intensity
    7: median signal intensity
    8: median absolute deviation of signal intensity
    0-3: mean signal intensity for each quartile
    """
    # for explanation, s. https://github.com/hiruna72/squigualiser/blob/main/docs/move_table.md
    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)
    
    # dictionary in form (query position: [event 1,2,3...,9]); for query position in rev_loci +- extra window (e.g. rev_locus(1842) +- 9)
    dict_events = {}

    # takes reverse locus and adds positions around it (given by extra window)
    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)
    
    # for signal summary, extracts features
    if mode == "signal_stats":

         # normalise signal
        median = np.median(signal)
        mad = np.median(np.abs(signal-median))
        signal=(signal-median)/mad

        for locus in pos_get_signal:
                
            prev = move_index[locus]*stride+offset
            sig_end = move_index[locus+1]*stride+offset

            sig_len = sig_end-prev
            data_tmp= np.zeros((9))
            data_tmp[4]=np.log10(sig_len)
            data_tmp[5]=np.mean(signal[prev:sig_end])
            data_tmp[6]=np.std(signal[prev:sig_end])
            data_tmp[7]=np.median(signal[prev:sig_end])
            data_tmp[8]=np.median(np.abs(signal[prev:sig_end]-data_tmp[4]))
            
            # get the mean signal for each quarter of the base signal
            for j in range(4):
                tmp_cnt=0
                for t in range(j*sig_len//4,min(sig_len, (j+1)*sig_len//4)):
                    data_tmp[j]+=signal[t+prev]
                    tmp_cnt+=1
                data_tmp[j]=data_tmp[j]/tmp_cnt
            dict_events.update({locus:data_tmp})
        return dict_events

    # for plotting individual read
    elif mode == "single_read":
        # seq2mv: [signal start, signal end, position in query] --> is used for plotting single read plot
        seq2mv = []

        # get positions in signal where signal for individual base starts/ends; builds array
        # this array only includes the relevant positions needed for plotting, hence locus is needed to later find the corresponding query position
        for locus in pos_get_signal:
            prev = move_index[locus]*stride+offset
            sig_end = move_index[locus+1]*stride+offset
            sig_len = sig_end-prev
            seq2mv.append([prev, sig_end, locus])
        
        return seq2mv, rev_loci

    else: raise Exception("mode has to be either \'signal_stats\' or \'single_read\'")
