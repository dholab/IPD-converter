#!/usr/bin/env nextflow

nextflow.enable.dsl = 2




// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_gbk = Channel
		.fromPath( "${params.input_data}/*.gb*" )

	ch_fasta = Channel
		.fromPath( "${params.input_data}/*.fasta" )
	
	ch_curated_embl = Channel
		.fromPath( params.curated_embl )
	
	// Workflow steps 
	FIX_GENBANK_DEFLINES ( 
		ch_gbk
	)
	
	FIX_FASTA_DEFLINES ( 
		ch_fasta
	)
	
	CONVERT_GENBANK_TO_EMBL ( 
		FIX_GENBANK_DEFLINES.out
	)
	
	ADJUST_TEXT ( 
		CONVERT_GENBANK_TO_EMBL.out,
		ch_fasta
	)
	
	REMOVE_OC_FEATURES ( 
		ADJUST_TEXT.out.oc
	)
	
	OBSCURE_ID ( 
		REMOVE_OC_FEATURES.out
	)
	
	SPLIT_EMBL_FILE ( 
		OBSCURE_ID.out.ipd
			.mix (
				ch_curated_embl
			)
	)
	
	TAR_EMBL_FILES ( 
		SPLIT_EMBL_FILE.out.collect()
	)
	
	
}
// --------------------------------------------------------------- //




// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.indiv_embl = params.results + "/" + params.experiment_number + "-ipd"
// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process FIX_GENBANK_DEFLINES {
	
	// In this process, the workflow corrects any definition sections in 
	// the GenBank file where the allele/locus name is missing. Typically,
	// this manifests itself as a series of possible alleles that are not
	// preceded by the locus name and a pipe symbol. For most alleles, no
	// definition sections will need to be corrected, but if you are exp-
	// orting the GanBank file from Geneious, occasionally, incomplete 
	// definition sections will arise as an artifact.
	
	// publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path gbk
	
	output:
	path "*_corrected.gb"
	
	script:
	"""
	genbank_fixer.R ${gbk}
	"""
}

process FIX_FASTA_DEFLINES {
	
	// This process simply removes the pipe symbol and any text that fol-
	// lows it, leaving only the allele/locus name 
	
	// publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path fasta
	
	output:
	path "*.fasta"
	
	script:
	"""
	fasta_fixer.R ${fasta}
	"""
}

process CONVERT_GENBANK_TO_EMBL {
	
	// Here we convert genbank input file to embl file. EMBL file will be 
	// adjusted in subsequent steps
	
	// publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path gbk
	
	output:
	path "*.emb"
	
	script:
	"""
	convert_genbank_to_embl.py ${gbk} "${params.species}"
	"""
}

process ADJUST_TEXT {
	
	// Here we adjust the text in various fields to match IPD requirements.
	// This is based on an example IPD EMBL file provided to Roger and John 
	// Caskey. We save the initial, adjusted EMBL file as oc.emb, since it is 
	// the one that the O'Connor Host Genomics Team and others will want to 
	// use for QC.
	
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path embl
	path fasta
	
	output:
	path "oc.emb", emit: oc
	path "error.emb", emit: error
	
	script:
	"""
	adjust_text.py ${embl} ${fasta} ${params.species}
	"""
}

process REMOVE_OC_FEATURES {
	
	// IPD does not want Geneious and MES contig identifiers in source anno-
	// tations
	
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path embl
	
	output:
	path "*.emb"
	
	script:
	"""
	remove_oc_features.py ${embl}
	"""
}

process OBSCURE_ID {
	
	// IPD has specific requires for the IPD line formatting, so we remove our
	// own information here so they can fill it in later.
	
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path embl
	
	output:
	path "ipd.emb", emit: ipd
	
	script:
	"""
	obscure_id.py ${embl}	
	"""
}

process SPLIT_EMBL_FILE {
	
	// Here, we use a formatting-agnostic method (i.e., not Biopython)
	// to split the EMBL file containing all novel sequences into one 
	// file per sequence, as is required by the IPD submission process.
	
	publishDir params.indiv_embl, mode: 'copy', overwrite: true
	
	cpus 1
	
	input:
	path ipd_emb
	
	output:
	path "*.embl"
	
	when:
	ipd_emb.exists()
	
	script:
	"""
	ipd_embl_splitter.R ${ipd_emb} ${params.experiment_number}
	"""
}

process TAR_EMBL_FILES {
	
	// Finally, we now compress all individual EMBL files into one
	// tarball, as is required by the IPD submission process.
	
	publishDir params.results, mode: 'copy', overwrite: true
	
	cpus 1
	
	input:
	path embl_files
	
	output:
	path "*.tar.gz"
	
	script:
	"""
	tar -czf ${params.experiment_number}.tar.gz *.embl
	"""
}

// --------------------------------------------------------------- //
