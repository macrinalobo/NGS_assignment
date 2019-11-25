import sys
import matplotlib.pyplot as plt
import math


def get_mid_5_mer_and_rev_compl(sequence):
    seq_len = len(sequence)
    if seq_len <= 5:
        return sequence
    elif seq_len % 2 != 0:
        midpoint_of_seq = int((seq_len - 1)/2)
        mid_sequence = sequence[midpoint_of_seq-2:midpoint_of_seq+3]
        return mid_sequence, get_reverse_complement(mid_sequence)
    else:  # never happens for this input since all sequence lengths are odd
        midpoint_of_seq = seq_len/2 - 1  # lower mid point considered if even
        lower_mid_sequence = sequence[midpoint_of_seq - 2:midpoint_of_seq + 3]
        upper_mid_sequence = sequence[midpoint_of_seq-1:midpoint_of_seq+4]
        # return comma separated 2 kmers for even length sequence
        return '{0},{1}'.format(lower_mid_sequence, upper_mid_sequence), '{0},{1}'.format(
            get_reverse_complement(lower_mid_sequence), get_reverse_complement(upper_mid_sequence))


def get_reverse_complement(sequence):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    # assume complement of 'any base' i.e. 'N' is 'N'
    complement_seq = []
    for base in sequence:
        if base not in complement_dict:  # For Ex: there are some 'Y' bases in the fasta
            # in this dataset it doesnt appear to happen for the mid 5 mer so ignore for now;
            # would need special handling for reverse complement of each non A,G,T,C entry
            print("this is a weird base {} in sequence {}".format(base,sequence))
            sys.exit(1)
        else:
            complement_seq.append(complement_dict[base])
    return ''.join(complement_seq[::-1])


def get_all_sorted_two_mers():
    # generates a list of all possible 2 base combinations of A,G,T,C sorted alphabetically
    bases = ['A', 'C', 'G', 'T']  # sorted alphabetically
    two_mers = []
    for base1 in bases:
        for base2 in bases:
            two_mers.append(base1 + base2)
    return two_mers


def get_vector_c(sequence, line_number):
    two_mer_list = get_all_sorted_two_mers()
    vector_c = [0] * len(two_mer_list)
    for count in range(len(sequence) - 1):
        # print(count)
        if sequence[count:count+2] not in two_mer_list:
            print("ignored {} in vector_C count for line {} at position {} of fasta".format(sequence[count:count+2],
                                                                                            line_number, count))
        else:
            vector_c_idx = two_mer_list.index(sequence[count:count+2])
            vector_c[vector_c_idx] = vector_c[vector_c_idx] + 1
    return '\t'.join([str(mer_count) for mer_count in vector_c])


def main():

    # QUESTION 1

    count = 0
    processed_read = {}
    read_lengths = []
    vector_c = []
    with open("Sequencing_Analysis_Assignment.fasta","r") as fasta, open("output_1a","w") as f_out:
        for line in fasta:
            count = count + 1
            line = line.strip()
            if line[0] == '>':  # identifier
                processed_read['read_id'] = line.split(' ')[0][1:]
                # print(processed_read['read_id'])
            else:  # base sequence
                read_lengths.append(len(line))
                processed_read['mid_5mer'], processed_read['reverse_comp_mid_5mer'] = get_mid_5_mer_and_rev_compl(line)
                vector_c.append(get_vector_c(line, count))
                f_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(processed_read['read_id'],
                            read_lengths[-1], processed_read['mid_5mer'],
                            processed_read['reverse_comp_mid_5mer'], vector_c[-1]))

    # PLOTS
    num_bins = int(math.sqrt(len(read_lengths))) + 1
    fig1 = plt.figure("histogram of read lengths")
    plt.style.use('ggplot')
    plt.hist(read_lengths, bins=num_bins)
    plt.xlabel("read length")
    plt.ylabel("number of reads")
    plt.title("histogram of read lengths")
    # plt.draw()
    fig1.savefig("read_length_histogram.png")

    # plot of cdf
    fig2 = plt.figure("cumulative distribution of 2mer frequencies")
    ax = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
    plt.xlabel("two-mer")
    xtick_labels = get_all_sorted_two_mers()
    ax.set_xticks(range(len(xtick_labels)))
    ax.set_xticklabels(xtick_labels)
    plt.ylabel("cumulative frequency")
    plt.title("cumulative distribution of 2mer frequencies for each sequence")

    for genomic_sig_c_str in vector_c:
        genomic_sig_c_list = genomic_sig_c_str.split('\t')
        pseudo_counts_p = [float(two_mer_cnt) + 1 for two_mer_cnt in genomic_sig_c_list]
        sum_l = sum(pseudo_counts_p)
        seq_2mer_freq_f = [px/sum_l for px in pseudo_counts_p]
        sequence_cdf = [seq_2mer_freq_f[0]]
        for idx in range(1, len(seq_2mer_freq_f)):
            sequence_cdf.append(seq_2mer_freq_f[idx] + sequence_cdf[idx-1])
        # print(seq_2mer_freq_f, sequence_cdf)
        plt.plot(sequence_cdf)

    # plt.draw()
    fig2.savefig("cumulative_distribution_of_2mer_frequencies.png")


if __name__ == '__main__':
    main()





