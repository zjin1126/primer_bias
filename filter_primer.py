"""
Filter primer pair by converage
"""
import argparse
import itertools

from Bio.SeqUtils import MeltingTemp as mt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'primer_seq', help='primer sequences in fasta format'
    )
    parser.add_argument(
        'coverage', help='coverage result from calculate_coverage.py'
    )
    parser.add_argument(
        'binding_site', help='binding site from binding_site.py'
    )
    parser.add_argument('--tax', nargs=1, default=['Bacteria,Archaea,Eukaryota'])
    args = parser.parse_args()
    
    tax_set = set(args.tax[0].split(','))

    # read primer seq
    primer_seq = {}
    with open(args.primer_seq) as fh:
        for line in fh:
            if line.startswith('>'):
                primer_name = line.rstrip()[1:]
                seq = fh.readline()
                primer_seq[primer_name] = seq.rstrip()

    # read primer coverage
    primer_cov = {}
    with open(args.coverage) as fh:
        for line in fh:
            if line.startswith('>'):
                primer_name = line.rstrip()[1:]
                bact = arch = euk = 0
                if 'Bacteria' in tax_set:
                    bact = float(fh.readline().rstrip().split()[2])
                if 'Archaea' in tax_set:
                    arch = float(fh.readline().rstrip().split()[2])
                if 'Eukaryota' in tax_set:
                    euk = float(fh.readline().rstrip().split()[2])
                fh.readline()  # discard last line
                primer_cov[primer_name] = {
                    'Bacteria': bact, 'Archaea': arch, 'Eukaryota': euk}

    # read primer binding sites
    sites = {}
    with open(args.binding_site) as fh:
        for line in fh:
            primer_name, start, end = line.rstrip().split()
            start = int(start)
            end = int(end)
            sites[primer_name] = (start, end)

    # filter primer by converage
    filtered_primer = []
    for primer in primer_cov:
        if primer_cov[primer]['Bacteria'] >= 0.85 or primer_cov[primer]['Archaea'] >= 0.85:
            filtered_primer.append(primer)

    # filter primer pair by delta Tm
    for primer_pair in itertools.combinations(filtered_primer, 2):
        pp = sorted(primer_pair, key=lambda x: sites[x][0])
        if sites[pp[1]][1] - sites[pp[0]][0] <= 900:  # filter short primer pair
            continue
        if abs(mt.Tm_NN(primer_seq[pp[0]], strict=False) - mt.Tm_NN(primer_seq[pp[1]], strict=False)) > 5:
            continue
        print('\t'.join(pp))


if __name__ == '__main__':
    main()
