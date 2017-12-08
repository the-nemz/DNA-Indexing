from __future__ import print_function
import sys
import time
from pybloomfilter import BloomFilter

filename = 'plain.bloom'
bloom = BloomFilter(500, .2, filename)


def main():
    global bloom, filename

    text = sys.stdin.readlines()
    # reads = [line.strip() for line in text[1::4]]
    reads = [line.strip() for line in text[1:]]

    load_bloom(reads)

    print('Finished creating', filename)


def load_bloom(reads):

    start = time.time()
    for read in reads:
        bloom.add(read)

    ms = 1000 * (time.time() - start)
    print('Loading BloomFilter took %.2f milliseconds.' % (ms))


def print_usage_error():
    print('USAGE: $ python2 < example.fa',
          '\n\twhere example.fa is a .fa file.')


if __name__ == "__main__":

    if len(sys.argv) == 1:
        main()
    else:
        print_usage_error()

