#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pathlib
import pandas as pd

embl_file = sys.argv[1]

with open("ipd.emb", "w") as output_handle:
    with open(embl_file, "r") as input_file:
        for line in input_file:
            if line.startswith('ID'):
                output_handle.write('ID   XXX; XXX; linear; XXX; XXX; XXX;\n')
            else:
                output_handle.write(line)
