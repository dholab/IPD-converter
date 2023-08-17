#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pathlib
import pandas as pd

annotated_genbank = sys.argv[1]

# Specify full species name for use in DE and organism tags
species = sys.argv[2]

with open("convert-gbk.emb", "w") as output_handle:
    for record in SeqIO.parse(annotated_genbank, "genbank"):
        SeqIO.write(record, output_handle, "embl")
