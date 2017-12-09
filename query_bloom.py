from __future__ import print_function
import sys
import time
from pybloomfilter import BloomFilter

filename = '.bloom'
kmer_length = 0
error_rate = 0
bloom = None
trues = 0
falses = 0


def main():
    global bloom, filename, trues, falses

    text = sys.stdin.readlines()
    # selects the reads from the fastq file
    reads = [line.strip() for line in text[1::4]]

    start = time.time()
    for read in reads:
        search_read(read)

    ms = 1000 * (time.time() - start)
    frac = 1.0 * trues / (trues + falses)

    with open('results/' + filename, 'w') as file:
        file.write('Searching reads took an average of %.3f milliseconds.'
                   % (ms / len(reads)))
        file.write('\n%.2f%% of k-mers were found in the BloomFilter.'
                   % (frac * 100))

    print('Finished analysis.')


def search_read(read):
    """
    This splits the read into its k-mer and checks to see
    if each is in the bloom filter
    If so, add to trues
    Otherwise add to falses
    """
    global bloom, trues, falses

    kmers = []
    index = len(read) - kmer_length
    while index >= 0:
        kmer = read[index:index + kmer_length]
        if kmer in bloom:
            trues += 1
        else:
            falses += 1
        index -= 1
    return kmers


def setup():
    global bloom, kmer_length, filename

    if len(sys.argv) < 2:
        return False

    try:
        bloom = BloomFilter.open(sys.argv[1])
        kmer_length = int(sys.argv[2])
    except ValueError:
        return False

    fname = sys.argv[1].split('/')[-1]
    ind = fname.index('.bloom')
    if ind < 1:
        return False
    filename = fname[:ind] + '.dat'

    return True


def print_usage_error():
    print('USAGE: $ python2 query_bloom.py example.bloom 20 < reads.fastq',
          '\n\twhere example.bloom is the already built bloom filter,'
          '\n\t20 is the kmer length, and'
          '\n\treads.fastq is a .fastq file.')


if __name__ == "__main__":

    if setup():
        main()
    else:
        print_usage_error()
