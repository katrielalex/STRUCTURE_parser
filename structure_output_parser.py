#! /usr/bin/env python
"""Written by Katriel Cohn-Gordon (2016).

STRUCTURE is population genetics software which has a terrible,
non-machine-readable output format.

This script attempts to translate its output format into something
that can actually be used by non-humans.
"""

import argparse
import pandas as pd
import re

try:
    import io
except ImportError:
    # Python 2
    import cStringIO as io

description = "Parse the (horrible) output of STRUCTURE into a machine-readable format"
parser = argparse.ArgumentParser(description)
parser.add_argument("structure_output_file", help="Output STRUCTURE file")

class StructureResults(dict):
    def __init__(self):
        self._facts = {}

    def add(self, description, fact):
        self._facts[description] = fact

class StructureResultsBlock(object):
    def __init__(self, name, regex, n_lines_after_regex):
        self.name = name
        self.matcher = re.compile(regex)
        self.n_lines_after_regex = n_lines_after_regex
        self.raw_lines = []

    def match(self, value):
        return self.matcher.search(value)

    def append(self, line):
        self.raw_lines.append(line)

    def to_dataframe(self):
        data = io.StringIO("".join(self.raw_lines))
        df = pd.read_table(data, delim_whitespace=True, header=None)
        return df

    def is_end(self, line):
        return re.match('^-*$', line)

def readStructureValuesFrom(lines):
    results = StructureResults()
    numeric_values = {name: re.compile(value) for name, value in {
        'indivs' : r'^(\d+) individuals.*$',
        'loci'   : r'^(\d+) loci.*$',
        'k'      : r'^(\d+) populations assumed.*$',
        'burnin' : r'^(\d+) Burn\-in period.*$',
        'reps'   : r'^(\d+) Reps.*$',
        # nan, inf
        'lnprob' : r'^Estimated Ln Prob of Data\s+=\s+([\d.$nainf-]+).*$',
        'meanln' : r'^Mean value of ln likelihood\s+=\s+([\d.$nainf-]+).*$',
        'varln'  : r'^Variance of ln likelihood\s+=\s+([\d.$nainf-]+).*$',
    }.items()}  # Pinched from StructureHarvester
    blocks = [
        StructureResultsBlock("ClusterMembership", "Proportion of membership", 5),
        StructureResultsBlock("AFDivergence", "Net nucleotide distance", 3),
        StructureResultsBlock("Heterozygosity", "expected heterozygosity", 1),
    ]

    currentBlock = None
    skipLines = 0

    for line in lines:
        # Skip lines if we were told to do so
        if skipLines:
            skipLines -= 1
            continue
        
        # If we're in the middle of a block, just read the next line blindly
        if currentBlock:
            if currentBlock.is_end(line.strip()):
                results[currentBlock.name] = currentBlock.to_dataframe()
                currentBlock = None
            else:
                currentBlock.append(line)
            continue
            
        # Check if we're reading a numeric value
        for name, value in numeric_values.items():
            match = value.match(line)
            if match:
                results[name] = match.group(1)
                break

        # Check if we're reading something that should start a block
        for block in blocks:
            if block.match(line):
                currentBlock = block
                skipLines = block.n_lines_after_regex
                continue

    return results

    
def main(args):
    with open(args.structure_output_file) as f:
        results = readStructureValuesFrom(f)

    print(results)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)

