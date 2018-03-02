"""
Get binding site from converted USEARCH result
"""
import argparse


def read_usearch_result(input_file: str) -> dict:
    usearch_result = {}
    with open(input_file) as fh:
        for line in fh:
            target, primer, start, end = line.rstrip().split()
            if primer not in usearch_result:
                usearch_result[primer] = {target: (start, end)}
            else:
                usearch_result[primer][target] = (start, end)
    return usearch_result


def binding_site_stat(sites: dict) -> tuple:
    site_stat = {}
    for _, site in sites.items():
        if site not in site_stat:
            site_stat[site] = 1
        else:
            site_stat[site] += 1
    site_stat = sorted(site_stat.items(), key=lambda x: x[1], reverse=True)
    return site_stat[0][0]


def get_binding_site(result: dict, reference_id: str) -> None:
    for primer in result:
        if reference_id in result[primer]:
            print(primer, *result[primer][reference_id], sep='\t')
        else:
            site = binding_site_stat(result[primer])
            print(primer, *site, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('usearch_result', help='USEARCH resault (converted)')
    parser.add_argument('reference_id', help='reference sequence id')
    args = parser.parse_args()
    result = read_usearch_result(args.usearch_result)
    get_binding_site(result, args.reference_id)


if __name__ == '__main__':
    main()
