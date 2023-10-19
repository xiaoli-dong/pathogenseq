#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/pathogen samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,long_fastq,basecaller_mode
    SAMEA6451102,read_1.fastq.gz,read_2.fastq.gz,longread.fastq.gz,r1041_e82_400bps_hac_v4.2.0

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """
    medaka_current_models = [
       
         # r1041 e82 (kit14) consensus
        'r1041_e82_400bps_hac_v4.2.0',
        'r1041_e82_400bps_sup_v4.2.0',
    ]
    medaka_archived_models = [
        # r9 consensus
        'r941_sup_plant_g610',
        'r941_min_fast_g507', 'r941_prom_fast_g507',
        'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
        'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
        'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
        'r941_prom_high_g344', 'r941_prom_high_g360', 'r941_prom_high_g4011',
        # r10 consensus
        'r10_min_high_g303', 'r10_min_high_g340',
        'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
        'r103_fast_g507', 'r103_hac_g507', 'r103_sup_g507',
        # r104 e81 consensus
        'r104_e81_fast_g5015', 'r104_e81_sup_g5015', 'r104_e81_hac_g5015',
        'r104_e81_sup_g610',
       
        # r1041 e82 consensus
        'r1041_e82_400bps_hac_g615',  'r1041_e82_400bps_fast_g615',
        'r1041_e82_400bps_fast_g632', 'r1041_e82_260bps_fast_g632',
        'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_g615',
        'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
        'r1041_e82_400bps_hac_v4.0.0', 'r1041_e82_400bps_sup_v4.0.0',
        'r1041_e82_260bps_hac_v4.0.0', 'r1041_e82_260bps_sup_v4.0.0',
        'r1041_e82_260bps_hac_v4.1.0', 'r1041_e82_260bps_sup_v4.1.0',
        'r1041_e82_400bps_hac_v4.1.0', 'r1041_e82_400bps_sup_v4.1.0',
        
        # rle consensus
        'r941_min_high_g340_rle',
        # r9 consensus
        'r941_min_hac_g507', 'r941_min_sup_g507',
        'r941_prom_hac_g507', 'r941_prom_sup_g507',
        
        # r941 e81 consensus
        'r941_e81_fast_g514', 'r941_e81_hac_g514', 'r941_e81_sup_g514'
       
    ]
    
    medaka_allowed_models = sorted(medaka_current_models + medaka_archived_models)

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        import re
        regex=re.compile('^#')
        

        ## Check header
        MIN_COLS = 2
        # TODO nf-core: Update the column names for the input samplesheet
        HEADER = ["sample", "fastq_1", "fastq_2", "long_fastq", "basecaller_mode"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)
        
        ## Check sample entries
        for line in fin:
            if re.match(regex, line):
                continue
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )
            
            ## Check sample name entries
            sample, fastq_1, fastq_2, long_fastq, basecaller_mode = lspl[: len(HEADER)]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2, long_fastq]:
                
                if fastq:
                    if fastq.upper() == 'NA':
                        continue
                    elif fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    elif not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Check basecalling mode
            if basecaller_mode:
                #if basecaller_mode.upper() != 'NA' and basecaller_mode.lower() not in ["fast", "hac", "sup"]:
                if basecaller_mode.upper() != 'NA' and basecaller_mode.lower() not in medaka_allowed_models:
                    modelStr = ' '.join(medaka_allowed_models)
                    print_error(
                        f"Long read bascaling mode is not valid, it can only be one of: {modelStr}",
                        "Line",
                        line,
                    )
            
            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["False", fastq_1, fastq_2, long_fastq, basecaller_mode]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["True", fastq_1, fastq_2, long_fastq, basecaller_mode]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2", "long_fastq", "basecaller_mode"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join(["{}_T{}".format(sample, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
