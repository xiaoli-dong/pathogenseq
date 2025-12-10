import csv
import sys

def extract_columns(file_name, columns):
    with open(file_name, mode='r') as file:
        csv_reader = csv.reader(file)
        header = next(csv_reader)  # Skip header
        
        # Convert column indices to integers
        columns = [int(col) - 1 for col in columns]  # Convert to 0-based index
        
        # Extract the specified columns
        for row in csv_reader:
            selected_columns = [row[col] for col in columns]
            
            # if  ''.join(selected_columns) == '':
            #     print(row)
            print("*****".join(selected_columns))

if __name__ == "__main__":
    # Expect file name as the first argument and column indices as the rest
    if len(sys.argv) < 3:
        print("Usage: python extract_columns.py <csv_file> <column_indices>")
        sys.exit(1)

    file_name = sys.argv[1]  # The CSV file name
    columns = sys.argv[2:]  # Column indices (1-based)

    extract_columns(file_name, columns)
