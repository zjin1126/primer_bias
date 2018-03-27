"""
Merge result from usearch global or local
"""
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+', help='Output file from USEARCH')
    args = parser.parse_args()
    
    partial_match = {}
    for input_file in args.input_files:
        with open(input_file) as fh:
            for line in fh:
                col = line.rstrip().split('\t')
                col[0] = col[0].split()[0]
                col[1] = col[1].split()[0]
                if float(col[6]) < 90:
                    continue
                if col[0] not in partial_match:
                    partial_match[col[0]] = col
                elif float(partial_match[col[0]][2]) < float(col[2]):
                        partial_match[col[0]] = col
    for _, match in partial_match.items():
        print('\t'.join(match))

if __name__ == '__main__':
    main()