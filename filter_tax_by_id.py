"""
filter SILVA taxonomy file by IDs
"""
import sys


def main():
    # read ids
    ids = set()
    with open(sys.argv[1]) as fh:
        for line in fh:
            ids.add(line.rstrip())
    # read tax file
    with open(sys.argv[2]) as fh:
        for line in fh:
            silva_id = line.rstrip().split('\t')[0]
            if silva_id in ids:
                print(line, end='')


if __name__ == '__main__':
    main()
