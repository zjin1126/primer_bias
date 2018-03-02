"""
Filter primer pair by converage
"""
import argparse
import itertools


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'coverage', help='coverage result from calculate_coverage.py')
    parser.add_argument(
        'binding_site', help='binding site from binding_site.py')
    args = parser.parse_args()

    primer_cov = {}
    with open(args.coverage) as fh:
        for line in fh:
            if line.startswith('>'):
                primer_name = line.rstrip()[1:]
                bact = float(fh.readline().rstrip().split()[2])
                arch = float(fh.readline().rstrip().split()[2])
                euk = float(fh.readline().rstrip().split()[2])
                fh.readline()  # discard last line
                primer_cov[primer_name] = {
                    'Bacteria': bact, 'Archaea': arch, 'Eukaryota': euk}
    sites = {}
    with open(args.binding_site) as fh:
        for line in fh:
            primer_name, start, end = line.rstrip().split()
            start = int(start)
            end = int(end)
            sites[primer_name] = (start, end)

    filtered_primer = []
    for primer in primer_cov:
        if primer_cov[primer]['Bacteria'] >= 0.75 or primer_cov[primer]['Archaea'] >= 0.75:
            filtered_primer.append(primer)

    for primer_pair in itertools.combinations(filtered_primer, 2):
        pp = sorted(primer_pair, key=lambda x: sites[x][0])
        if sites[pp[1]][1] - sites[pp[0]][0] >= 900:  # only long primer pair
            print('\t'.join(pp))


if __name__ == '__main__':
    main()
