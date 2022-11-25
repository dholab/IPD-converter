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
7. **Rerun IPD converter on the EMBL file with animal names inserted, if desired. Note that this step is not necessary—it simply splits the full EMBL into individual files and compresses them in a tarball**.
8. [Submit to IPD](https://www.ebi.ac.uk/ipd/submission/).

## Quick start
Please not that this converter requires that you have Docker, Git, Java, and NextFlow installed on your system. If you do not have these dependencies installed, proceed to the **Detailed Setup Instructions** section below.

Additionally, this converter has only been tested on x86 systems and was developed with NextFlow version 22.10.0.5826. Mileage may vary on other compute architectures or Nextflow versions.

With those disclaimers out of the way, here is an example invocation of the converter that does not require you to `git clone` this repo:

```
nextflow dholab/IPD-converter -latest \
--input_data data/ \
--species "Macaca mulatta" \
--experiment_number 1
```
By default, the converter looks for a _single_ FASTA and a single GenBank in a subdirectory of where you launched it called `data`. However, you may place the FASTA and GenBank files anywhere as long as you supply the absolute path after the `--input_data` flag.

The experiment number gives users the opportunity to use numbers from their own experiment and file tracking systems. However, if the user is not working with such a system, simply use this argument to name the final folder where the converted files will be placed.

## Detailed Setup Instructions

First, make sure you are using an POSIXct-compatible system (e.g. MacOS, Linux, [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install)) with Git installed.

Then, download the converter files by running `git clone` in a directory of your choice:

```
git clone https://github.com/dholab/AVRL-pango-updator.git .
```

Next, the Docker engine must be installed if it isn't already. To do so, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/) to find instructions for your system.

Most systems already java installed; run `java --version` in a terminal to check if yours does as well. If java is not installed, visit [Oracle's Java installation page](https://www.oracle.com/java/technologies/downloads/).

### Nextflow Installation

This pipeline was built with the [NextFlow](https://www.nextflow.io/) pipeline manager (v22.10.0.5826). We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. If you do not have install conda installed, go to the following link to download Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line, which will make conda usage much faster:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

If you would prefer not to use conda at all, you can also install NextFlow with curl:

1. Run the following line in a directory where you'd like to install NextFlow:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, NextFlow has successfully installed.

Once all the above dependencies are installed, you are ready to run the converter.

### Running the pipeline

The pipeline is invoked from the working directory where you downloaded it with a shell command, which looks something like this:

```
nextflow run main.nf \
--input_data data/ \
--species "Macaca mulatta" \
--experiment_number 27909
```

If you are re-running the converter to split and tar a manually edited EMBL file, simply do so with:

```
nextflow run main.nf \
--curated_embl data/curated.embl \
--input_data data/ \
--species "Macaca mulatta" \
--experiment_number 27909
```

You may change `curated.embl` to whatever filename you've used and place it anywhere on your system as long as you provide the absolute file path.

