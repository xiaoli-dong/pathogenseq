#!/usr/bin/env python

import csv
import argparse
import xml.etree.ElementTree as ET

def xml_to_csv(element, csv_writer): 
    # Extract element name and text 
    name = element.tag 
    text = element.text 
 
    # Extract element attributes 
    attrib = element.attrib 
 
    # Write element name and text as well as attributes to CSV file 
    row = [name, text] + list(attrib.values()) 
    csv_writer.writerow(row) 
 
    # Recursively process child elements 
    for child in element: 
        xml_to_csv(child, csv_writer) 

def main():
    
    description = "Combine multiple xml files into a single xml file"
    parser = argparse.ArgumentParser(description=description)

    # help=f"Space seperated xml file list, for example: 'f1.xml f2.xml f3.xml'\n",
    parser.add_argument('-i', "--input", required=True, help=f"space seperated xml file name list\n")
    parser.add_argument("-o", "--output", required=True, default="combined.xml", help=f"Output file name\n")

    args = parser.parse_args()

    cols = ["name", "phone", "email", "date", "country"] 
    rows = [] 

    tree = ET.parse(args.input)
    root = tree.getroot()
    
    for result in root[1]:
        sample_data = []
        #print(result.tag)
        #print(result.attrib)
        for detail in result:
            #print(detail.tag)
            #print(detail.attrib)
            print(detail.attrib.get("type"))
            print(detail.attrib.get("value"))
            
    
    #Open CSV file for writing 
    with open(args.output, "w", newline="") as csv_file: 
        # Create CSV writer 
        csv_writer = csv.writer(csv_file) 
    
        # Convert XML to CSV 
        xml_to_csv(root, csv_writer) 
    
if __name__ == "__main__":
    main()

