#!/usr/bin/env python3

import sys


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


def fix_header_and_write(sequences):
    for k, v in sequences.items():
        sys.stdout.write('>{}\n{}\n'.format(
            k.replace(' ', '_'),
            v
        ))


if __name__ == "__main__":
    S = {k: v for k, v in parse_fasta(sys.argv[1])}
    fix_header_and_write(S)
