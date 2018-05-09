#!/usr/bin/env python
"""
Calculate taxonomy coverage from converted USEARCH search_oligodb
"""
import argparse
import pickle
import sys

from tax_node import TaxNode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tax_tree', help='Taxonomy tree')
    parser.add_argument(
        'tax_tid', help='Taxonomy TID (species acc to tree node id)')
    parser.add_argument('ref_site', help='Reference binding site')
    parser.add_argument(
        'usearch_result', help='Position converted USEARCH output')
    parser.add_argument('--deepth', default=1, help='Deepth of lineage')
    args = parser.parse_args()

    # read binding site
    print('Read binding site from {}...'.format(args.ref_site), file=sys.stderr)
    ref_site = {}
    with open(args.ref_site) as fh:
        for line in fh:
            primer_name, start, end = line.rstrip().split('\t')
            ref_site[primer_name] = (int(start), int(end))
    print('Done!', file=sys.stderr)

    # print('Read taxonomy from {}...'.format(args.tax_tree), file=sys.stderr)
    with open(args.tax_tid, mode='rb') as fh:
        acc_tid = pickle.load(fh)
    # print('Done!', file=sys.stderr)
    
    print('Read USEARCH result from {}...'.format(
        args.usearch_result), file=sys.stderr)

    hit_result = {}
    # hit = set()
    with open(args.usearch_result) as fh:
        for line in fh:
            acc, primer_name, start, end = line.rstrip().split('\t')
            if primer_name not in hit_result:
                hit_result[primer_name] = set()
            start = int(start)
            end = int(end)
            if (abs(ref_site[primer_name][0] - start) <= 1 and
                    abs(ref_site[primer_name][1] - end) <= 1):
                # hit.add(acc_tid[acc])
                hit_result[primer_name].add(acc_tid[acc])

    for primer in hit_result:
        with open(args.tax_tree, mode='rb') as fh:
            root_node = pickle.load(fh)

        stack = [root_node]
        while stack:
            pre = stack[-1]
            if 'child_done' in pre.data:  # child_done, so cal
                pre.data['count'] = sum([x.data['count'] for x in pre.child])
                pre.data['cov'] = sum([x.data['cov'] for x in pre.child])
                pre.data['done'] = True
                stack.pop()
                continue
            count = 0
            for curr in pre.child:
                if 'done' not in curr.data:
                    if not curr.child:  # child empty
                        if (curr.data['start'] <= ref_site[primer][0] and
                            curr.data['end'] >= ref_site[primer][1]):
                            curr.data['count'] = 1
                            if curr.tid in hit_result[primer]:
                                curr.data['cov'] = 1
                            else:
                                curr.data['cov'] = 0
                        else:
                            curr.data['count'] = 0
                            curr.data['cov'] = 0
                        curr.data['done'] = True
                        continue  # cal next child
                    stack.append(curr)
                    break
                count += 1
            if count == len(pre.child):
                pre.data['child_done'] = True
        print('Done!', file=sys.stderr)
        print('>' + primer)
        for c in root_node.child:
            print('{}\t{}\t{:.4f}'.format(
                c.data['name'], c.data['cov'], c.data['cov'] / c.data['count']))
        not_match = []
        for child in root_node.child[0].child:
            if child.data['cov'] == 0:
                not_match.append(child.data['name'])
        print('\t'.join(not_match))


if __name__ == '__main__':
    main()
