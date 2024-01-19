#!/usr/bin/env python

import argparse
import json

def main():

    description = "Combine multiple json files into a single json file"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"Comma seperated json file list, for example: 'f1.json,f2.json,f3.json'\n",
    )
    parser.add_argument("-o", "--output", required=True, default="combined.json", help=f"Output file name\n")
    
    args = parser.parse_args()
    json_files = args.input.split(sep=',')

    # Create a list of all the JSON files that you want to combine.
    #json_files = ["file1.json", "file2.json", "file3.json"]

    # Create an empty list to store the Python objects.
    python_objects = []

    # Load each JSON file into a Python object.
    for json_file in json_files:
        print(json_file)
        with open(json_file, "r") as fin:
            data = json.load(fin)
            #python_objects.append(json.load(fin, strict=False))
            python_objects.append(data)
    # Dump all the Python objects into a single JSON file.
    with open(args.output, "w") as fout:
        json.dump(python_objects, fout, indent=4)
           
    fout.close()
   

if __name__ == "__main__":
    main()
