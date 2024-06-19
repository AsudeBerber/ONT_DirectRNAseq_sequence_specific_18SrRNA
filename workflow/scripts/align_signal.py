__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

def get_events(signal, moves, offset, rev_loci, motif_length=1, extra_window=21, signal_stats = False):
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

    # normalise signal
    median = np.median(signal)
    mad = np.median(np.abs(signal-median))
    signal=(signal-median)/mad
    
    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)
    
    
    dict_events = {}

    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)
    
    if signal_stats == True:
        for locus in pos_get_signal:
                
            prev = move_index[locus]*stride+offset
            sig_end = move_index[locus+1]*stride+offset

        
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

    elif signal_stats == False:
        pass

    else: raise Exception("signal_stats has to be either True or False")