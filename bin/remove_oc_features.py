#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pathlib
import pandas as pd

embl_file = SeqIO.parse(sys.argv[1], "embl")

with open("remove_oc.emb", "w") as output_handle:
    for record in embl_file:
        for feature in record.features:
            # check if type is source
            if feature.type == 'source':
                # remove Geneious ID and MES cluster ID
                feature.qualifiers.pop("mes_cluster")
                feature.qualifiers.pop("geneious_id")
            # if feature.type == 'exon':
            #     # remove notes from Geneious
            #     feature.qualifiers.pop("note")

        # write edited sequence to file
        SeqIO.write(record, output_handle, "embl")
