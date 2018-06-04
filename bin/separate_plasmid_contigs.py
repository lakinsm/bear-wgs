#!/usr/bin/env python3

import sys
import re


length_regex = re.compile(r'length=(\d+)')


def parse_fasta(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return


def parse_mummer_out(infile, minlen):
    ret = set()
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        contig_name = ""
        contig_len = 0
        for line in data:
            if not line:
                continue
            if line.startswith('>'):
                continue
            else:
                contig_name = line.split()[0]
                print(contig_name)
                contig_len = int(length_regex.search(line).group(1))
                if contig_len > minlen:
                    ret.add(contig_name)
    return ret


def separate_contigs(S, A, plasmid_out, contig_out):
    with open(plasmid_out, 'w') as pout, open(contig_out, 'w') as cout:
        for header, seq in S.items():
            if header in A:
                pout.write('>{}\n{}\n'.format(header, seq))
            else:
                cout.write('>{}\n{}\n'.format(header, seq))


if __name__ == '__main__':
    spades = {k: v for k, v in parse_fasta(sys.argv[1])}
    alignments = parse_mummer_out(sys.argv[2], int(sys.argv[5]))
    separate_contigs(spades, alignments, sys.argv[3], sys.argv[4])
