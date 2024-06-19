__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

def get_events(signal, moves, offset, rev_loci, motif_length, extra_window, signal_stats = False):
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
    
    for locus in pos_get_signal:
            
        prev = move_index[locus]*stride+offset
        sig_end = move_index[locus+1]*stride+offset

    if signal_stats == True:
        