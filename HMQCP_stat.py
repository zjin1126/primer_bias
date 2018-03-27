"""
HMQCP converage calculator
"""

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('otu_silva', help='search HMQCP against SILVA databsae (in blast6 format)')
    parser.add_argument('usearch_output', help='Position converted USEARCH output')
    parser.add_argument('primer_pair', help='primer pairs in tsv format')
    parser.add_argument('ref_site', help='Reference binding site')
    args = parser.parse_args()

    # read HMP OTU information (search against SILVA)
    total_otu = 0
    silva_otu = {}
    with open(args.otu_silva) as fh:
        for line in fh:
            otu, silva_target = line.rstrip().split()[:2]
            if silva_target in silva_otu:
                silva_otu[silva_target] += 1
            else:
                silva_otu[silva_target] = 1
            total_otu += 1

    # read primer pair
    primer_pairs = []
    with open(args.primer_pair) as fh:
        for line in fh:
            fwd, rev = line.rstrip().split()
            primer_pairs.append((fwd, rev))

    # read primer binding site
    ref_site = {}
    with open(args.ref_site) as fh:
        for line in fh:
            primer_name, start, end = line.rstrip().split('\t')
            ref_site[primer_name] = (int(start), int(end))

    # read USEARCH result
    usearch_result = {}
    with open(args.usearch_output) as fh:
        for line in fh:
            silva_target, primer, start, end = line.rstrip().split()
            if primer not in usearch_result:
                usearch_result[primer] = set()
            start = int(start)
            end = int(end)
            if (abs(ref_site[primer][0] - start) <= 1 and
                abs(ref_site[primer][1] - end) <= 1):
                usearch_result[primer].add(silva_target)
    
    # calculate coverage
    for primer_pair in primer_pairs:
        count = 0
        for otu in silva_otu:
            if otu in usearch_result[primer_pair[0]] and otu in usearch_result[primer_pair[1]]:
                count += silva_otu[otu]
        print('>' + primer_pair[0] + '+' + primer_pair[1])
        print('Total\t{}\t{}'.format(count, count / total_otu))


if __name__ == '__main__':
    main()