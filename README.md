# Nextflow pipeline for converting novel allele sequences to Immuno Polymorphism Database format

## Overview
This converter was meant to solve a routine problem for the [AVRL](https://dholk.primate.wisc.edu/wiki/home/page.view?name=home_index) Host Genomics Team: Converting putative novel allele sequences, discovered in macaque MHC genes, to the EMBL-like format required by the [Immuno Polymorphism Database](https://www.ebi.ac.uk/ipd/) (IPD). This is a crucial step in documenting new MHC sequences that are relevant to biomedical research around the world. If we can automate this step, IPD's uploading process can run with fewer errors, and  novel sequences can thus become available to the research community more quickly.

This converter is part of a broader allele discovery and curation process, described below. The bold text shows where this converter is used.

1. Generate PacBio amplicon data for novel allele assembly.
2. Run novel allele discovery analysis, like what we have created and documented for [discovering novel MHC Class II alleles](https://github.com/dholab/MHC-II-allele-discovery).
3. Conduct expert, manual review of these alleles as well as intron/exon annotations in Geneious or other genetic data visualization software.
4. Export novel allele sequences, with annotations, into GenBank flat file format. At the same time, export the sequences _without annotations_ in FASTA format.
5. **Run IPD converter on the GenBank and FASTA file. Instructions below.**
6. Review these converted files and add animal or isolate names to each review.
7. **Rerun IPD converter on the EMBL file with animal names inserted, if desired. Note that this step is not necessary—it simply splits the full EMBL into individual files and compresses them in a tarball.
8. [Submit to IPD](https://www.ebi.ac.uk/ipd/submission/).

## Quick start
Please not that this converter requires that you have Docker, Git, Java, and NextFlow installed on your system. If you do not have these dependencies installed, proceed to the **Detailed Instructions** section below.

Additionally, this converter has only been tested on x86 systems and was developed with NextFlow version 22.10.0.5826. Mileage may vary on other compute architectures or Nextflow versions.

With those disclaimers out of the way, here is an example invocation of the converter that does not require you to `git clone` this repo:

```
nextflow dholab/IPD-converter -latest \
--input_data data/ \
--species "Macaca mulatta" \
--experiment_number 1
```

The experiment number gives users the opportunity to use numbers from their own experiment and file tracking systems. However, if the user is not working with such a system, simply use this argument to name the final folder where the converted files will be placed.
