#!/usr/bin/env python3

import os
import logging
import re
from pathlib import PurePath as pp


# Set up logging module for basic log level
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(module)s - %(message)s', datefmt='%Y-%m-%d %H:%M')
logger = logging.getLogger('varanno')

def converttofilename(sample_name):
    """
    Convert to a valid filename.
    Removes leading/trailing whitespace, converts internal spaces to underscores.
    Allows only alphanumeric, dashes, underscores, unicode.
    :param sample_name: A sample file name to remove offending characters
    :return: A stripped, valid file name
    """
    return re.sub(r'(?u)[^-\w]', '', sample_name.strip().replace(' ', '_'))

def getcurrentlocation():
    """
    Set the analysis file path location regardless of where the python script is located.
    """
    anpath = os.path.abspath(__file__)
    parpath = pp(anpath)
    srcdir = str(parpath.parents[1])
    return srcdir

def getvcfpath(inputvcf):
    """
    The goal of this function is to determine whether we have a properly formatted vcf and the vcf version
    :param inputvcf:
    :return: vcf version, path to the vcf
    """
    try:
        if os.path.isfile(inputvcf):
            vcfpath = os.path.abspath(inputvcf)
            logger.info("Input VCF location: {}".format(vcfpath))
            return vcfpath
        elif os.path.exists(os.path.join(getcurrentlocation(), inputvcf)):
            vcfpath = os.path.join(getcurrentlocation(), inputvcf)
            logger.info("Input VCF location: {}".format(vcfpath))
            return vcfpath
        else:
            logger.error("Could not determine the location of the vcf file")

    except FileNotFoundError:
        logger.error("Exiting because the input vcf could not be found")
        exit(1)

def getoutputpath(outputloc):
    """
    The goal of this function is to determine whether the output location exists
    :param: outputloc is the output location as a path
    """
    try:
        if os.path.exists(outputloc):
            outpath = os.path.abspath(outputloc)
            logger.info("Using output location: {}".format(outpath))
            return outpath
        else:
            outpath = os.getcwd()
            logger.info("Using output location: {}".format(outpath))
            return outpath

    except NotADirectoryError:
        logger.error("Could not determine an output location or identify the current working directory")
        logger.error("Exiting.")
        exit(1)


def getvcfversion(inputvcf):
    """
    The goal of this function is to determine whether we have a properly formatted vcf and the vcf version
    :param inputvcf:
    :return: vcf version, path to the vcf
    """
    try:
        with open(getvcfpath(inputvcf), "r") as vcfreader:
            vcfheader = vcfreader.readline()
            if vcfheader.startswith('##'):
                vcfver = re.sub('[#\n]', '', vcfheader)
                logging.info("VCF formatted as {}".format(vcfver))
                return vcfver
            else:
                logging.error("This file does not appear to be a properly formatted vcf")

    except FileNotFoundError:
        logging.error("Exiting because the input vcf could not be found")
        exit(1)


def infodict(infofield):
    """
    Function to parse and convert the INFO field for a given variant call to a dictionary
    :param infofield: INFO string to be parsed
    :return: INFO field as a dictionary
    """
    infofield = infofield.split(';')
    infodict = dict(map(lambda s: s.split('='), infofield))
    return infodict