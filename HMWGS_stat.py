"""
HMP MetaPhlan2 stat
"""
import argparse

from tax_node import TaxNode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'species_list', help='HMP MetaPhlan2 species list in genus level')
    parser.add_argument('silva_tax', help='SILVA taxonamy')
    args = parser.parse_args()

    genus_list = set()
    with open(args.species_list) as fh:
        for line in fh:
            rank = line.rstrip().split('|')
            if rank[5].endswith('_unclassified') or rank[5].endswith('_noname'):
                continue
            genus = rank[5][3:].lower().replace('_', ' ')
            genus_list.add(genus)

    gene_list = []
    with open(args.silva_tax) as fh:
        for line in fh:
            silva_id, tax = line.rstrip().split('\t')
            tax = tax.split(';')
            if len(tax) < 6:
                continue
            genus = tax[5].lower()
            if genus in genus_list:
                gene_list.append(silva_id)
            else:
                if genus == 'escherichia-shigella':
                    gene_list.append(silva_id)
                elif genus == 'burkholderia-paraburkholderia':
                    gene_list.append(silva_id)
                elif genus.startswith('coprococcus'):
                    gene_list.append(silva_id)
                elif genus.startswith('clostridium sensu stricto'):
                    gene_list.append(silva_id)
                elif genus.startswith('ruminococcus'):
                    gene_list.append(silva_id)
    print('\n'.join(gene_list))


if __name__ == '__main__':
    main()
