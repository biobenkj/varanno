#!/usr/bin/env python3

"""
This is a small python VCF annotator script that will query ExAC, providing severity predictions and a couple other things (e.g. ExAC allele freq). The script itself is capable of taking
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
chromosome: chromosome that the variant was called on
position: chromosomal coordinate of the variant
ref_allele: the reference allele
alt_allele: the alternative allele called
seq_depth: sequencing depth/number of reads covering that loci (DP VCF tag)
reads_supporting_variant: the number of reads supporting the alternative allele/variant (AO VCF tag)
percent_alt_relative_to_total_support: percentage of the total reads that support the variant/alt allele (AO / DP)
variant_type: the type of variant that was called (e.g. snp, ins, del, etc.)
variant_effect: the effect the variant has within its genomic context - the effect reported is the most severe as
estimated by Ensembl based on the sequence ontology terms (if the ExAC variant consequence term doesn't exist in the
sequence ontology list or is not reported, the string 'No_matching_severity_term' is returned)
ExAC_allele_freq: the allele frequency found in the ExAC database (if known) for the called variant - if the allele
frequency is not reported, the string 'Unknown_in_ExAC' is returned

Author: Ben K. Johnson 2018
"""

import os
import argparse
import logging
import csv
import requests
import json

# Import utilities from utils.py
from utils import *

import version

# Set the ExAC bulk variant query API URL
EXAC_API_URL = 'http://exac.hms.harvard.edu/rest/bulk/variant'

# Set the sequence ontology ranking as a dictionary (e.g. 1 is most severe and 35 is least severe)
# Source of ranking is estimated from Ensembl using terms from the Sequence Ontology consortium
# http://www.ensembl.org/info/genome/variation/predicted_data.html
SEVERITY_RANKING = {
    'transcript_ablation':1,
    'splice_acceptor_variant':2,
    'splice_donor_variant':3,
    'stop_gained':4,
    'frameshift_variant':5,
    'stop_lost':6,
    'start_lost':7,
    'transcript_amplification':8,
    'inframe_insertion':9,
    'inframe_deletion':10,
    'missense_variant':11,
    'protein_altering_variant':12,
    'splice_region_variant':13,
    'incomplete_terminal_codon_variant':14,
    'stop_retained_variant':15,
    'synonymous_variant':16,
    'coding_sequence_variant':17,
    'mature_miRNA_variant':18,
    '5_prime_UTR_variant':19,
    '3_prime_UTR_variant':20,
    'non_coding_transcript_exon_variant':21,
    'intron_variant':22,
    'NMD_transcript_variant':23,
    'non_coding_transcript_variant':24,
    'upstream_gene_variant':25,
    'downstream_gene_variant':26,
    'TFBS_ablation':27,
    'TFBS_amplification':28,
    'TF_binding_site_variant':29,
    'regulatory_region_ablation':30,
    'regulatory_region_amplification':31,
    'feature_elongation':32,
    'regulatory_region_variant':33,
    'feature_truncation':34,
    'intergenic_variant':35
}

# Set up the logging
logger = logging.getLogger('varanno')
logger.setLevel(logging.INFO)


def parse_args():

    """Parse command line arguments for the varanno.py script."""

    parser = argparse.ArgumentParser(description='Variant annotation script for the technical challenge. Requires python 3.5 or greater',
                                     usage='python varanno.py [options]',
                                     epilog='Written by Ben K Johnson, PhD\nVARI Bioinformatics and Biostatistics Core\n2018',
                                     prog='varanno')

    parser.add_argument("-v", "--version", help="Installed VarAnno version",
                        action="version",
                        version="%(prog)s " + str(version.__version__))

    parser.add_argument('-o', '--output',
                        help='Output directory within which to place the final results. (default=current directory)',
                        required=False, default=os.getcwd())

    parser.add_argument('-i', '--inputvcf',
                        help='Alternative VCF input. Defaults to local Challenge_data.vcf from Git repo.',
                        required=False, default='Challenge_data.vcf')
    # Currently non-functional due to time restrictions
    #parser.add_argument('-t', '--tests',
    #                    help='Turn on unit testing suite. This will run tests and then exit. (default=False)',
    #                    required=False, type=bool, default=False, choices=[True, False])

    return parser.parse_args()


class VcfParser(object):
    """
    A class for reading and parsing VCFs
    """
    def __init__(self, vcf):
        self.reader = self.vcfreader(vcf)  # Call the reader constructor

    @classmethod  # Build a class method that can be called with the dot operator
    def vcfreader(cls, inputvcfpath):
        """
        A function that will read in a VCF file from inputvcfpath and will not parse metadata. The returned object is
        a generator with which we can iterate over each variant call by row.
        :param inputvcfpath: path to a VCF file
        :return: generator to iterate over variant calls by row
        """
        try:
            with open(inputvcfpath, "r") as vcfreader:  # Open in read mode
                for row in csv.reader(vcfreader, dialect="excel-tab"):  # Treat this as a tab separated file
                    if len(row):  # This will address whether there is an empty VCF
                        if row[0].startswith('##') or row[0].startswith('#'):  # Skip the header metadata
                            continue  # Skip the vcf header information
                        else:
                            yield row  # Grab each row and return as an iterator
                    else:
                        logger.warning("This appears to be either an empty VCF or an empty row!")

        except IOError:
            logger.error("Couldn't open the vcf.")  # If the VCf cannot be opened, exit with status 1
            logger.error("Exiting")
            exit(1)

class ExacQuery(object):
    """
    Class to query the ExAC API with a bulk payload instead of hitting it for each variant.
    """
    def __init__(self, json_variant_info):
        # Query the ExAC db with a bulk payload in the form of a json array
        variantanno = requests.post(EXAC_API_URL, data=json.dumps(json_variant_info))

        # Check to see if the request was successful
        if variantanno.status_code != 200:
            # If the status code is not 200 (success), raise an error in the logging and quit
            logger.error("The ExAC API appears to be down. Try again later.")
            exit(1)
        else:
            # If the status code indicates success, continue with annotations
            self.variantannojson = variantanno.json()  # Convert the output to a list of dictionaries

        super().__init__()  # Remove need to ref base class

    def varianteffect(self, variant):
        """
        A function to parse the consequence field returned for a given variant from ExAC. It attempts to catch as many
        possible issues that I could think of if either a value doesn't exist in the current SO terms or no known
        consequence.
        :param variant: an ExAC dictionary
        :return: the most severe SO term (if multiple) or place holder string if no term is found
        """
        consequence = variant['consequence']  # Get the consequence field, if exists from the ExAC variant query id

        if consequence and not consequence == None:  # Check if dict is empty or a NoneType
            so_val = [SEVERITY_RANKING.get(severity) for severity in consequence]  # Get value by SO key
            if so_val == [None]:  # Check if the key does not exist in the SO list
                so_val = -1  # Set the SO value to a number that does not exist in the SO list
            else:
                so_val = min([seval for seval in so_val if seval is not None])  # Get the most severe SO value
        else:
            so_val = -1  # Capture the variants with no known consequence in ExAC or not in SO list

        if not so_val in range(1, 35):  # Check if the SO value exists in the SO list above
            so_term = "No_matching_severity_term"  # If not exist, annotate with 'No_matching_severity_term'
        else:
            # A way to get keys from values
            so_term = list(SEVERITY_RANKING.keys())[list(SEVERITY_RANKING.values()).index(int(so_val))]
        return so_term

    def allelefreq(self, variant):
        """
        A function to return the ExAC allele frequency for a given variant and return a place holder string if no value
        is found.
        :param variant: an ExAC dictionary
        :return: the ExAC allele frequency for a given variant if it exists
        """
        return variant['variant'].get('allele_freq', 'Unknown_in_ExAC')


def main():
    #--Launch the tool--#

    logger.info("Welcome to VarAnno.")
    logger.info("Parsing command line arguments")

    args = parse_args()  # Parse command line arguments

    #--Locate and read input--#

    logger.info("Locating input...")

    ivcf = getvcfpath(args.inputvcf)  # Find the location of the input VCF

    logger.info("Reading the vcf...")

    vcf = VcfParser(ivcf)  # set up the VCF reader generator

    logger.info("Done reading vcf")

    #--Check output location--#

    logger.info("Checking output location...")

    outpath = getoutputpath(args.output)  # Check if the output location exists. If not, fall back to PWD

    #--Build up the JSON array for ExAC query--#

    logger.info("Generating the json array, decomposing multiallelic variants, and gathering variant info...")

    exac_query = []  # Set up an empty list to be populated by variant values to query the ExAC API

    variant_anno_info = []  # Set up an empty list to be populated for gathering additional variant information/stats

    for row in vcf.reader:
        chromosome = row[0]  # CHROM field
        position = row[1]  # POS field
        refallele = row[3]  # REF field
        if len(row[4].split(',')) > 1:  # Check if the ALT field is multiallelic
            altallele = row[4].split(',')  # Decompose ALT field
            for aa in range(len(altallele)):  # Decompose the called multiple ALT alleles at the genomic location
                exac_query.append('-'.join((chromosome, position, refallele, altallele[aa])))  # ExAC query specific format
                # Get additional pertinent information as we go and append to another list
                seqdepth = infodict(row[7])['DP']  # Get the total sequencing depth at that location
                altcounts = infodict(row[7])['AO'].split(',')[aa]  # Get the ALT counts for each allele call
                if int(altcounts) != int(seqdepth):  # Prevent division by 0
                    pctalttotoal = 100 * (int(altcounts) / (int(seqdepth)))  # ALT read count support as pct of total
                else:
                    pctalttotoal = 100  # This checks for division by 0 and if alt == ref counts, all reads support alt
                vartype = infodict(row[7])['TYPE'].split(',')[aa]  # Get the variant call type (e.g. snp, del, ins, etc)
                variant_anno_info.append([chromosome,
                                          position,
                                          refallele,
                                          altallele[aa],
                                          seqdepth,
                                          altcounts,
                                          pctalttotoal,
                                          vartype])  # Append to the list
        else:
            altallele = row[4]  # If biallelic
            exac_query.append('-'.join((chromosome, position, refallele, altallele[0])))  # ExAC query specific format
            # Get additional pertinent information as we go and append to another list
            seqdepth = infodict(row[7])['DP']  # Get the total sequencing depth at that location
            altcounts = infodict(row[7])['AO']  # Get the ALT counts
            if int(altcounts) != int(seqdepth):  # Prevent division by 0
                pctalttotoal = 100 * (int(altcounts) / (int(seqdepth)))  # ALT read count support as pct of total
            else:
                pctalttotoal = 100  # This is a check for division by 0 and if alt == ref counts, all reads support alt
            vartype = infodict(row[7])['TYPE']  # Get the variant call type (e.g. snp, del, ins, etc)
            variant_anno_info.append([chromosome,
                                      position,
                                      refallele,
                                      altallele,
                                      seqdepth,
                                      altcounts,
                                      pctalttotoal,
                                      vartype])  # Append to the list

    logger.info("Done generating json array, decomposing multiallelic variants, and gathering variant info.")

    #--Query and parse ExAC data--#

    logger.info("Querying ExAC...")

    exac_data = ExacQuery(exac_query)  # Query ExAC API with JSON array

    logger.info("Received response and data from ExAC.")

    logger.info("Annotating variants with ExAC data...")

    exac_annotations = []  # Initialize an empty list to be populated with parsed annotations from ExAC

    for variant in exac_query:
        exac_info = exac_data.variantannojson[str(variant)]  # Find the corresponding variant in the ExAC list of dicts
        exac_consequence = exac_data.varianteffect(exac_info)  # Get the most severe SO term
        exac_af = exac_data.allelefreq(exac_info)  # Get the ExAC allele frequency if it exists
        exac_annotations.append([variant, exac_consequence, exac_af])  # Append to the list

    logger.info("Finished annotating variants.")

    logger.info("Writing data out to {}".format(os.path.join(outpath, "Challenge_data.annotated.tsv")))

    # Write out the data
    with open(os.path.join(outpath, "Challenge_data.annotated.tsv"), "w") as tsv:
        # Write data in a tab-delimited format
        writer = csv.writer(tsv, dialect="excel-tab")
        # Write the header row for the annotated variants
        writer.writerow(['chromosome',
                         'position',
                         'ref_allele',
                         'alt_allele',
                         'seq_depth',
                         'reads_supporting_variant',
                         'percent_alt_relative_to_total_support',
                         'variant_type',
                         'variant_effect',
                         'ExAC_allele_freq'])
        try:
            # This is imperative to check since we add annotations from one list to another
            # and need to make sure they line up
            assert len(variant_anno_info) == len(exac_annotations)

        except ValueError:
            # Exit if lists are not of equal length. Something strange must have happened.
            logger.error("Variant annotations are not consistent.")
            logger.error("Something strange happened during the run.")
            logger.error("Try running in test mode to narrow down the error - <python varanno.py -t>")
            logger.error("Exiting")
            exit(1)

        for annotated_variants in range(len(variant_anno_info)):
            # Add ExAC annotation data to variant info/stats
            # Need to append each value one at a time since a slice returns a list
            variant_anno_info[annotated_variants].append(exac_annotations[annotated_variants][1])
            variant_anno_info[annotated_variants].append(exac_annotations[annotated_variants][2])
            # Write out the data by row - this is in a list format
            writer.writerow(variant_anno_info[annotated_variants])


    logger.info("Finished!")
    return

if __name__ == '__main__':
    # Launch the prototype tool
    main()
