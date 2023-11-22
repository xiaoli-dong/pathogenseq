#!/usr/bin/env python
import click
import os
import csv

@click.command()
@click.option("-n", "--sample-name", default="", help="Sample Name.")
@click.option("-o", "--output-csv", default="", help="output csv file")
@click.option("-s", "--sep", default="\t", help="csv file separator")
@click.argument("input_csv", nargs=-1)
def reformat_csv(sample_name, output_csv, input_csv, sep):
    
    for csvfile in input_csv:
        with open(csvfile, "r", encoding="utf8") as in_file:
            csvreader = csv.DictReader(in_file, delimiter=sep)
            header = csvreader.fieldnames
            header[0] = "sampleid"
           
            rows = []
            for row in csvreader:
                row[header[0]] = sample_name
                rows.append(row)

            with open(output_csv, "w") as file:
                writer = csv.DictWriter(file, fieldnames=header, delimiter=sep)
                writer.writeheader()
                writer.writerows(rows)
        
if __name__ == "__main__":
    reformat_csv()
