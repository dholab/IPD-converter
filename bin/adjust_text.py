#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pathlib
import pandas as pd

embl_file = SeqIO.parse(sys.argv[1], "embl")
fasta_file = sys.argv[2]
genus = str(sys.argv[3])
species = str(sys.argv[4])
scientific_name = f"{genus} {species}"

with open("error.emb", "w") as error_handle:
	with open("oc.emb", "w") as output_handle:
		for record in embl_file:
			# extract information on matching MES cluster
			# if not matching MES cluster, alert user and send sequence to error file
			# iterate over putative alleles
			# check if putative allele sequence matches EMBL file sequence
			mes_clusters_file = SeqIO.parse(fasta_file, "fasta")
			
			# add source annotation
			source_feature = SeqFeature(FeatureLocation(0, len(record)), type="source")
			source_feature.qualifiers["organism"] = scientific_name
			source_feature.qualifiers["mol_type"] = 'genomic DNA'
			source_feature.qualifiers["geneious_id"] = record.id

			# get MES cluster that matches annotated sequence - lookup by sequence
			for sequences in mes_clusters_file:
				if record.seq == sequences.seq:
					# set mes_cluster property to mes_cluster ID
					source_feature.qualifiers["mes_cluster"] = sequences.id

					# extract animal name from mes cluster
					animal_name = source_feature.qualifiers["mes_cluster"].split('_')[0]

					# add source to record
					record.features.append(source_feature)

			# if there is not a matching source annotation, meaning no MES cluster has same sequence as annotated genbank
			# save sequence to error file and go to next record
			if "mes_cluster" not in source_feature.qualifiers:
				SeqIO.write(record, error_handle, "embl")
				continue
			
			# IPD allele name (extract from description)
			# allele_index = allele_name = record.description.split('-')[-1]
			allele_name = record.description.split('|')[0] # + "_" + allele_index

			# change organism to species, which is what IPD wants
			record.annotations['organism'] = scientific_name
			
			# move allele name from description to accession
			record.annotations['accession'] = allele_name
		
			# change description to full-length name requested by IPD using regular expressions
			record.description = scientific_name + ' gene for MHC class I antigen (' + animal_name + ' gene), isolate ' + animal_name + ', allele ' + allele_name
			
			# add decoration to CDS annotation
			# remove standard_name lint from Gneeious from CDS anntoation
			for feature in record.features:
				# check if type is CDS
				if feature.type == 'CDS':  
					feature.qualifiers['standard_name'] = "CDS"
					feature.qualifiers['gene'] = animal_name
					feature.qualifiers['allele'] = allele_name
			
			# write edited sequence to file
			SeqIO.write(record, output_handle, "embl")
