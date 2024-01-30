#!/usr/bin/env python

import argparse
import xml.etree.ElementTree as ET
import sys
def main():

    description = "Combine multiple xml files into a single xml file"
    parser = argparse.ArgumentParser(description=description)

    # help=f"Space seperated xml file list, for example: 'f1.xml f2.xml f3.xml'\n",
    parser.add_argument('-i', "--input", required=True, help=f"space seperated xml file name list\n")
    parser.add_argument("-o", "--output", required=True, default="combined.xml", help=f"Output file name\n")
    
    args = parser.parse_args()

    xml_files = args.input.split()

    with open(args.output, "a+") as fout:
        # Load each JSON file into a Python object.
        xml_element_tree = None
        for xml_file in xml_files:
            data = ET.tostring(ET.parse(xml_file).getroot()).decode("utf-8")
            fout.write(data)
            fout.write('\n')    
    fout.close()
   

if __name__ == "__main__":
    main()
