#!/usr/bin/env python

import argparse
import csv

def main():

    description = "add header to emmtyper output and get rid of .tmp from the isolate name"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"Comma seperated json file list, for example: 'f1.json,f2.json,f3.json'\n",
    )
    parser.add_argument(
        "-o", 
        "--output", 
        required=True, 
        default="", 
        help=f"Output file name\n"
    )
    parser.add_argument(
        "-s", 
        "--sname", 
        required=True, 
        default="", 
        help=f"delimiter\n"
    )
    parser.add_argument(
        "-d", 
        "--delimiter", 
        default="\t", 
        help=f"input and output delimiter\n"
    )
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf8") as f_input:
            with open(args.output, "w") as f_output:
                csvreader = csv.DictReader(f_input, delimiter=args.delimiter)
                header = csvreader.fieldnames
                header[0] = "sampleid"
            
                rows = []
                for row in csvreader:
                    row[header[0]] = args.sname
                    rows.append(row)

                
                writer = csv.DictWriter(f_output, fieldnames=header, delimiter=args.delimiter)
                writer.writeheader()
                writer.writerows(rows)
       
if __name__ == "__main__":
    main()