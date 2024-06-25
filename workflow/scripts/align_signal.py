__author__ = "Jens Martin"
__email__ = "jens.martin@outlook.com"

import numpy as np
import pdb
def common():
    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)

    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)

    for locus in pos_get_signal:
        
        prev = move_index[locus]*stride+offset
        sig_end = move_index[locus+1]*stride+offset
        sig_len = sig_end-prev

# for getting position in signal (goes 3' to 5')
def rev_locus(locus, read):
    locus_rev = read.query_length -1 - locus
    return locus_rev

def get_loci(read, pairs, wd, motif_length):
    """
    find positions that match motif
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

    stride = moves.pop(0)
    move_index = np.where(moves)[0]
    rlen = len(move_index)
    
    
    dict_events = {}

    # last position is only to get sig_end
    pos_get_signal = [np.arange((locus - extra_window), locus + motif_length+extra_window) for locus in rev_loci]
    pos_get_signal = np.reshape(pos_get_signal, -1)
    
    if signal_stats == True:

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

    elif signal_stats == False:
        pass

    else: raise Exception("signal_stats has to be either True or False")


def seq_to_mv(reads_ids, region, sample, seq=None, mv=None, ts=0, fetch = True, pos=42):
    if fetch == True:
        seq, mv, ts, pos_read, ref_seq = bam_aligned(sample, reads_ids, region, pos)


    
    seq = seq[::-1] #sequence order is 5' -> 3', mv and signal are 3' -> 5': therefore sequence is turned around
    s = stride
    p = 1 #itinerates through movetable array
    x = 0 #number of additional strides (stride amount - 1)
    start = ts + 1

    seq2mv = np.array([[1, ts, "-", "-"]])
    for i, base in enumerate(seq):
        # print (base)
        while p < len(mv)-1 and mv[p + 1] == 0: #last movetable index(p): mv[p+1] doesn't exist, would cause error
            x = x + 1 #counts 0's
            p = p + 1 #0 found -> movetable index moves by one
            # print (f"{p},{x}\n")
        p = p + 1 

        end = start + (stride-1) + (x*s) #stride 1: start + 4 (as start number is already first position in stride), all further strides: additional +5
        x = 0 #resets number of additional strides
        ref_base = ref_seq [i]
        seq2mv = np.append(values=[[start, end, base, ref_base]], arr=seq2mv, axis=0) # appends 
        start = end+1 #next base starts 1 after end of previous base

    print ("sequence-to-signal alignment finished,", len(seq), "bases, signal length =", end)    
    return seq2mv, pos_read


def reverse_seq_mv(seq2mv):
    k = int(seq2mv[-1][1]) #signal range
    rev_seq2mv = np.flip(seq2mv,0)
    for base_data in rev_seq2mv:
        start_old = int(base_data[0])
        end_old = int(base_data[1])
        base_data[0] = k - end_old #start and end switch when string is read in other direction, former last base is now first and vice versa
        base_data[1] = k - start_old
    return rev_seq2mv


def bam_aligned(sample, read_ids, region, pos):

    samfile = ps.AlignmentFile(f"{sample}")

    max_reads = samfile.mapped
    i = 0

    # if index file present, fetch(region = region)
    for read in samfile.fetch(region = region):
        if read.query_name == read_ids:
            read_ID = read.query_name
            seq = read.query_sequence
            

            # Workaround in cases where two ts tags per read exists:
            # read.set_tag("ts", None) #first ts tag is transcript strand(+|-), has to be removed

            mv = read.get_tag("mv")
            ts = read.get_tag("ts")

            ref_seq = read.get_reference_sequence()

            aln_pairs = read.get_aligned_pairs(with_seq = True)
            
            # creates pairs of base positions (query, reference) -> looks up position in alignment sequence for corresponding reference base position
            for pair in aln_pairs:
                if pair [1] == pos: 
                    pos_read = pair[0]
                
            # print(read.get_tags())
            # print(f"ts:{ts}; {read_ID}") 
            return(seq, mv, ts, pos_read, ref_seq)
        else: 
            #removing the else part makes the code only 1s faster
            i = i+1
            k = i/500000
            if k.is_integer():
                # print("currently at " +read.query_name + "\n")
                print(f"reads compared: {i} of max. {max_reads}")
            continue