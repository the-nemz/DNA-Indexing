from __future__ import print_function
import sys
import time
from pybloomfilter import BloomFilter

filename = '.bloom'
kmer_length = 0
error_rate = 0
bloom = None


def main():
    global bloom, filename

    text = sys.stdin.readlines()
    lines = [line.strip() for line in text[1:]]
    dna = ''.join(lines)

    start = time.time()
    kmers = make_kmers(dna)
    load_bloom(kmers)
    ms = 1000 * (time.time() - start)
    print('Loading BloomFilter took %.2f milliseconds.' % (ms))

    print('Finished creating', filename)
    print('It has %d bits and %d hashes.' % (bloom.num_bits, bloom.num_hashes))


def make_kmers(dna):
    """
    Builds the k-mers of the DNA
    """
    kmers = []
    index = len(dna) - kmer_length
    while index >= 0:
        kmer = dna[index:index + kmer_length]
        kmers.append(kmer)
        index -= 1
    return kmers


def load_bloom(kmers):
    """
    Inserts all of the k-mers into the bloom filter
    """
    global bloom, filename

    filename = '%d_kmer_%d_rate.bloom' % (kmer_length, int(100 * error_rate))
    print(len(kmers) // 2)
    bloom = BloomFilter(len(kmers) // 2, error_rate, filename)

    for kmer in kmers:
        bloom.add(kmer)


def setup():
    global kmer_length, error_rate

    if len(sys.argv) < 3:
        return False

    try:
        kmer_length = int(sys.argv[1])
        error_rate = float(sys.argv[2])
    except ValueError:
        return False

    return True


def print_usage_error():
    print('USAGE: $ python2 kmer_bloom.py 20 .01 < example.fa',
          '\n\twhere 20 is the kmer length,'
          '\n\t.01 is the error rate, and'
          '\n\texample.fa is a .fa file.')


if __name__ == "__main__":

    if setup():
        main()
    else:
        print_usage_error()
