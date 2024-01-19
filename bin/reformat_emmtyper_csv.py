#!/usr/bin/env python

import argparse
import csv

def main():

    description = "add header to emmtyper output and get rid of .tmp from the sampleid"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"emmtyper csv output'\n",
    )
    parser.add_argument(
        "-o", 
        "--output", 
        required=True, 
        default="", 
        help=f"Output file name\n"
    )
    args = parser.parse_args()

    header = [
        "sampleid", 
        "num_of_blast_hits",
        "num_of_clusters",
        "emm-type",
        "emm-type-positions",
        "emm-like",
        "emm-like-positions",
        "EMM-cluster"
    ]
    with open(args.input, "r", encoding="utf8") as f_input:
        with open(args.output, 'w', newline='') as f_output:
            reader = csv.reader(f_input, delimiter='\t')
            writer = csv.writer(f_output, delimiter=',')
            writer.writerow(header)
            for row in reader:
                row[0] = row[0].rstrip(".tmp")
                writer.writerow(row)
       
        
if __name__ == "__main__":
    main()