#!/usr/bin/env python
import sys
import gzip

#!/usr/bin/env python3

import sys
import gzip

def open_file(path, mode='rt'):
    """
    Opens a file normally or with gzip based on file extension.
    mode = 'rt' for reading text, 'wt' for writing text.
    """
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def parse_fastq(filepath):
    """
    Parses a FASTQ file and returns a dictionary mapping read IDs to full headers.
    """
    headers = {}
    with open_file(filepath, 'rt') as f:
        while True:
            header = f.readline()
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not header:
                break
            header = header.strip()
            read_id = header.split()[0].lstrip('@')
            headers[read_id] = header
    return headers

def replace_headers(source_fastq, target_fastq, output_fastq):
    """
    Replaces headers in target FASTQ using headers from source FASTQ, matching by read ID.
    Writes output to output_fastq.
    """
    source_headers = parse_fastq(source_fastq)
    with open_file(target_fastq, 'rt') as fin, open_file(output_fastq, 'wt') as fout:
        while True:
            header = fin.readline()
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not header:
                break
            read_id = header.strip().split()[0].lstrip('@')
            new_header = source_headers.get(read_id, header.strip())
            #fout.write(f"@{new_header if not new_header.startswith('@') else new_header}\n")
            fout.write(f"{new_header}\n")
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)

def main():
    if len(sys.argv) != 4:
        print("Usage: python replace_fastq_headers.py source.fastq[.gz] target.fastq[.gz] output.fastq[.gz]")
        sys.exit(1)
    
    source_fastq = sys.argv[1]
    target_fastq = sys.argv[2]
    output_fastq = sys.argv[3]

    replace_headers(source_fastq, target_fastq, output_fastq)

if __name__ == "__main__":
    main()
