"""
Wrapper of USEARCH
"""
import argparse
import os
import shlex
import subprocess


def do_search_oligodb(refs: str, oligodb: str, output_file: str, diff: str) -> None:
    command = 'usearch -search_oligodb {} -db {} -userout {} -strand plus ' + \
        '-userfields query+target+qstrand+diffs+tlor+thir+qlor+qhir+trowdots -maxdiffs {}'
    command = command.format(refs, oligodb, output_file, diff)
    args = shlex.split(command)
    subprocess.run(args)


def revcomp(input_file: str, output_file: str) -> None:
    command = 'usearch -fastx_revcomp {} -fastaout {}'
    command = command.format(input_file, output_file)
    args = shlex.split(command)
    subprocess.run(args)


def merge_file(input_files: list, output_file: str) -> None:
    fhout = open(output_file, mode='w+')
    for file in input_files:
        with open(file) as fh:
            for line in fh:
                fhout.write(line)
    fhout.close()


def detect_reverse(input_file: str) -> None:
    fhout_fwd = open('.tmp_fwd.fa', mode='w+')
    fhout_rev = open('.tmp.fa', mode='w+')
    with open(input_file) as fhin:
        for line in fhin:
            if line.startswith('>'):
                temp = line.rstrip().split('-')
                if temp[5] == 'A':  # antisense
                    print(line, end='', file=fhout_rev)
                    print(fhin.readline(), end='', file=fhout_rev)
                else:
                    print(line, end='', file=fhout_fwd)
                    print(fhin.readline(), end='', file=fhout_fwd)
    fhout_fwd.close()
    fhout_rev.close()
    if os.path.getsize('.tmp.fa') > 0:
        revcomp('.tmp.fa', '.tmp_rev.fa')
        merge_file(['.tmp_fwd.fa', '.tmp_rev.fa'], '.tmp_resolved.fa')
    else:
        os.rename('.tmp_fwd.fa', '.tmp_resolved.fa')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'reference', help='reference sequences in fasta format')
    parser.add_argument('primers', help='primersin fasta format')
    parser.add_argument('output', help='search_oligodb result')
    args = parser.parse_args()
    detect_reverse(args.primers)
    do_search_oligodb(args.reference, '.tmp_resolved.fa', args.output, 1)


if __name__ == '__main__':
    main()
