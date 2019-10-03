![Python 3.5](https://img.shields.io/badge/python-3.5-blue.svg)
![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)

# VarAnno
### Variant annotation tool

This is a small VCF annotator script written in Python 3 that will query ExAC, providing severity predictions and a couple other things (e.g. ExAC allele freq). The script itself is capable of taking
in a user provided VCF and outputting an annotated VCF. The annotation step involves reading in the variant calls,
followed by decomposition of multiallelic variant calls, and annotation with data from ExAC. The output is in the form
of a tab-delimited text file that can be opened in a text editor.

Requirements: requests>=2.18.1 and pathlib2>=2.1.0

Requirements can be installed with:
    `pip install -r requirements.txt`

Usage:
    `python varanno.py`

Input: User-provided VCF or defaults to a local Challenge_data.vcf
Output: Tab-delimited file with decomposed and annotated variants to a user defined location or current working directory

Output file format: 10 columns
1. chromosome: chromosome that the variant was called on

2. position: chromosomal coordinate of the variant

3. ref_allele: the reference allele

4. alt_allele: the alternative allele called

5. seq_depth: sequencing depth/number of reads covering that loci (DP VCF tag)

6. reads_supporting_variant: the number of reads supporting the alternative allele/variant (AO VCF tag)

7. percent_alt_relative_to_total_support: percentage of the total reads that support the variant/alt allele (AO / DP)

8. variant_type: the type of variant that was called (e.g. snp, ins, del, etc.)

9. variant_effect: the effect the variant has within its genomic context - the effect reported is the most severe as
estimated by Ensembl based on the sequence ontology terms (if the ExAC variant consequence term doesn't exist in the
sequence ontology list or is not reported, the string 'No_matching_severity_term' is returned)

10. ExAC_allele_freq: the allele frequency found in the ExAC database (if known) for the called variant - if the allele
frequency is not reported, the string 'Unknown_in_ExAC' is returned

Author: Ben K. Johnson 2018

TODO: add unit tests for test suite and use something like `tox`
