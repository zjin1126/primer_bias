#!/usr/bin/env python
import sys

import requests
from bs4 import BeautifulSoup


def get_primer(url):
    response = requests.get(url)
    if response.status_code != 200:
        print('Fail to connect!', file=sys.stderr)
        return
    soup = BeautifulSoup(response.text, 'lxml')
    primer_data = soup.select('table > tr')
    # get primer name
    primer_name_tag = primer_data[1].select('td')[1]
    primer_seq_tag = primer_data[9].select('td')[1]
    primer_seq = ''.join(primer_seq_tag.text.split()[1:-1])
    print(primer_name_tag.text, primer_seq, sep='\t')


def main():
    url = 'http://probebase.csb.univie.ac.at/pb_results/categories/26'
    response = requests.get(url)
    if response.status_code != 200:
        print('Fail to connect!', file=sys.stderr)
        return
    soup = BeautifulSoup(response.text, 'lxml')
    primer_url = soup.select('#myTable > tr > td > a')
    for primer_tag in primer_url:
        get_primer(primer_tag['href'])


if __name__ == '__main__':
    main()
