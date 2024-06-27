__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import numpy as np
import pdb
# for getting position in signal (goes 3' to 5')
def rev_locus(locus, read):
    locus_rev = read.query_length -1 - locus
    return locus_rev

def get_loci(read, pairs, wd, motif_length, ref_pos):
    """
    gets pos in query for reference position
    """    
    ref_loci = []
    ref_loci_index = []
    # reversed loci, for operations from 3' end
    rev_loci = []

    #pairs[0]: query pos; [1]: ref pos; [2] ref base
    for i, pos in enumerate(ref_pos):
        # for empty array
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
            
 
    loci = [pairs[locus, 0] for locus in ref_loci_index]
    # Remove loci that are not present on the query or too close to the ends of the alignment
    loci = [locus for locus in loci if locus is not None and locus > wd and locus < read.query_length - wd - (motif_length)]
    ref_loci = [pairs[index, 1] for index in ref_loci_index if pairs[index, 0] != None and pairs[index,0] in loci]
    if len(loci) != len(ref_loci):
        raise Exception ("length of reference and query sequence index not matching")
    
    for locus in loci:
        rev_loci.append(rev_locus(locus, read))
    return loci, ref_loci, rev_loci

# this part is modified from https://github.com/WGLab/DeepMod2/blob/main/src/detect.py
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

    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)
    breakpoint()
    
    dict_events = {}

    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)
    
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

    elif mode == "single_read":
        seq2mv = []
        breakpoint()
        for locus in pos_get_signal:
            prev = move_index[locus]*stride+offset
            sig_end = move_index[locus+1]*stride+offset
            sig_len = sig_end-prev
            seq2mv.append([prev, sig_end, locus])
        
        return seq2mv, rev_loci

    else: raise Exception("signal_stats has to be either True or False")
