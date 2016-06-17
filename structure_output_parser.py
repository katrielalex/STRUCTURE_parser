#! /usr/bin/env python
"""Written by Katriel Cohn-Gordon (2016).

STRUCTURE is population genetics software which has a terrible,
non-machine-readable output format.

This script attempts to translate its output format into something
that can actually be used by non-humans.
"""

import argparse
import re
import pandas as pd

try:
  import io
except ImportError:
  # Python 2
  import cStringIO as io

description = "Parse the (horrible) output of STRUCTURE into something nicer."
parser = argparse.ArgumentParser(description)
parser.add_argument("structure_output_file", help="Output STRUCTURE file")


class StructureResultsBlock(object):
  """Represents, effectively, a dataframe embedded in the STRUCTURE output."""

  def __init__(self, name, regex, n_lines_after_regex):
    self.name = name
    self.matcher = re.compile(regex)
    self.n_lines_after_regex = n_lines_after_regex
    self.raw_lines = []

  def match(self, value):
    return self.matcher.search(value)

  def append(self, line):
    self.raw_lines.append(line)

  def finish(self):
    data = io.StringIO("".join(self.raw_lines))
    df = pd.read_table(data, delim_whitespace=True, header=None)
    return df

  def is_end(self, line):
    return re.match("^-*$", line)


class InferredAncestryPop0Block(StructureResultsBlock):
  """Special-cased block for the inferred ancestry of population 0."""

  def __init__(self):
    StructureResultsBlock.__init__(self,
                                   "InferredAncestry",
                                   "Inferred ancestry",
                                   2)

  def append(self, line):
    self.raw_lines.append(line.split(":")[1].lstrip())

  def is_end(self, line):
    return line.split()[3] != "0"


def read_structure_from(lines):
  """Attempt to parse some sense into the STRUCTURE output file |lines|."""
  results = {}
  numeric_values = {name: re.compile(value) for name, value in {
      "indivs": r"^(\d+) individuals.*$",
      "loci": r"^(\d+) loci.*$",
      "k": r"^(\d+) populations assumed.*$",
      "burnin": r"^(\d+) Burn\-in period.*$",
      "reps": r"^(\d+) Reps.*$",
      # nan, inf
      "lnprob": r"^Estimated Ln Prob of Data\s+=\s+([\d.$nainf-]+).*$",
      "meanln": r"^Mean value of ln likelihood\s+=\s+([\d.$nainf-]+).*$",
      "varln": r"^Variance of ln likelihood\s+=\s+([\d.$nainf-]+).*$",
  }.items()}  # Pinched from StructureHarvester
  blocks = [
      StructureResultsBlock("ClusterMembership", "Proportion of membership", 5),
      StructureResultsBlock("AFDivergence", "Net nucleotide distance", 3),
      StructureResultsBlock("Heterozygosity", "expected heterozygosity", 1),
      InferredAncestryPop0Block(),
  ]

  current_block = None
  skip_lines = 0

  for line in lines:
    # Skip lines if we were told to do so
    if skip_lines:
      skip_lines -= 1
      continue

    # If we're in the middle of a block, just read the next line blindly
    if current_block:
      if current_block.is_end(line.strip()):
        results[current_block.name] = current_block.finish()
        current_block = None
      else:
        current_block.append(line)
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
        current_block = block
        skip_lines = block.n_lines_after_regex
        continue

  return results


def main(args):
  with open(args.structure_output_file) as f:
    results = read_structure_from(f)

  print(results['InferredAncestry'])

if __name__ == "__main__":
  main(parser.parse_args())

