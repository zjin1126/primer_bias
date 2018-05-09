# download SILVA database
wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/ARB_files/SSURef_NR99_128_SILVA_07_09_16_opt.arb.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/taxmap_slv_ssu_ref_nr_128.txt.gz
gzip -d *.gz
# use arb software to extract aligned sequences (SILVA_128_aligned.fasta)

# download primers from probebase
python scripts/get_primers.py > scripts/primers.txt
cat primers.txt | awk '{ printf ">%s\n%s\n", $1, $2 }' > primers.fa  # convert to fasta
# make fai index
samtools faidx SILVA_128_aligned.fasta
# extract reference sequence
python scripts/get_reference.py SILVA_128_aligned.fasta U00096.223771.225312 > ecoli_ref.fa
# extract position information
python scripts/extract_posision.py ecoli_ref.fa SILVA_128_aligned.fasta pos_info_new.pickle
# build taxonomy tree
python scripts/tax_node.py pos_info_new.pickle silva_tax.txt tax_tree.pickle tax_tid.pickle
# run usearch search_oligodb
python scripts/usearch.py SILVA_128_SSURef_Nr99_tax_silva.fasta scripts/primers.fa temp_result.txt
# convert unaligned position to aligned position
python scripts/convert_position.py pos_info_new.pickle temp_result.txt > temp_result_conv.txt
# get reference binding sites
python scripts/binding_site.py temp_result_conv.txt U00096.223771.225312 > temp_binding_site.txt
# calculate primer coverage
python scripts/calculate_coverage.py tax_tree.pickle tax_tid.pickle temp_binding_site.txt temp_result_conv.txt > result.txt
# filter primer (coverage > 85%)
python scripts/filter_primer.py scripts/primers.fa result.txt temp_binding_site.txt > primer_pair.txt
# calculate primer pair coverage
python scripts/primer_pair.py tax_tree.pickle tax_tid.pickle temp_binding_site.txt primer_pair.txt temp_result_conv.txt > pair_result.txt

# DEGEPRIME
# prepare
MakeSilvaTaxonomy.pl -i SILVA_128_aligned.fasta -o degeprime_tax.txt
TrimAlignment.pl -i SILVA_128_aligned.fasta -o all.fa -ref U00096.223771.225312
# filter (baceria)
fgrep "Bacteria;" degeprime_tax.txt > bacteria_tax.txt
python -c "
with open('bacteria_tax.txt') as fh:
    tax = set()
    for line in fh:
        name, des = line.split('\t')
        tax.add('>'+name+'\n')
with open('all.fa') as fh, open('bacteria.fa', mode='w') as fho:
    for line in fh:
        if line.startswith('>'):
            if line in tax:
                print(line, end='', file=fho)
                print(fh.readline(), end='', file=fho)
"
# run DegePrime.pl (Length: 16-25)
seq 16 25 | xargs -P 6 -I@ DegePrime.pl -i bacteria.fa -l @ -d 8 -o bacteria_primer_@.txt -taxfile bacteria_tax.txt -taxlevel 2
# select top 100 primers
tail -n +2 bacteria_primer_16.txt | cut -f 1-7 | awk '{print $1, $6/$2, $7}' | sort -nk2 | tail -n 100 | cut -f 1,3 -d ' ' | tr ' ' \t > bact_top100.txt
# calculate primer coverage (DEGEPRIME)
python scripts/usearch_custom.py SILVA_128_SSURef_Nr99_tax_silva.fasta primer_design/bact_top100.fa bact_top100_result.txt
python scripts/convert_position.py pos_info_new.pickle bact_top100_result.txt > bact_top100_result_conv.txt
python scripts/binding_site.py bact_top100_result_conv.txt U00096.223771.225312 > bact_top100_site.txt
python scripts/calculate_coverage.py tax_tree.pickle tax_tid.pickle bact_top100_site.txt bact_top100_result_conv.txt > bact_top100_coverage.txt
python scripts/filter_primer.py bact_top100.fa bact_top100_coverage.txt bact_top100_site.txt > bact_top100_pair.txt
python scripts/primer_pair.py tax_tree.pickle tax_tid.pickle bact_top100_site.txt bact_top100_pair.txt bact_top100_result_conv.txt > bact_top100_pair_coverage.txt

# HMQCP (16S)
wget http://downloads.ihmpdcc.org/data/HMQCP/rep_set_v35.fna.gz
# split database into 2 files (RAM limitation)
head -n 6000001 SILVA_128_SSURef_Nr99_tax_silva.fasta > test.1.fa
tail -n +6000002 SILVA_128_SSURef_Nr99_tax_silva.fasta > test.2.fa
# make udb
usearch -makeudb_usearch test.1.fa -output test.1.udb
usearch -makeudb_usearch test.2.fa -output test.2.udb
# usearch against SILVA SSURef database
usearch -usearch_local rep_set_v35.fna -db test.1.udb -strand both -id 0.9 -userfields query+target+id+alnlen+mism+opens+qcov+evalue+bits -userout HMP_usearch.1.txt
usearch -usearch_local rep_set_v35.fna -db test.2.udb -strand both -id 0.9 -userfields query+target+id+alnlen+mism+opens+qcov+evalue+bits -userout HMP_usearch.2.txt
python scripts/merge_usearch.py HMP_usearch.*.txt > HMP_usearch.txt
# calculate coverage
python scripts/HMQCP_stat.py HMP_usearch.txt temp_result_conv.txt primer_pair.txt temp_binding_site.txt > HMQCP_converage.txt

# HMSMCP2 (shotgun)
wget http://downloads.hmpdacc.org/data2/hhs/genome/microbiome/wgs/analysis/hmsmcp/v2/hmp1-II_metaphlan2-mtd-qcd.pcl.bz2
# filter of genus level data
grep 'g__' hmp1-II_metaphlan2-mtd-qcd.pcl | grep -v 's__' | egrep -v '(k__Viruses|k__Eukaryota)' | cut -f 1 > hmwgs_list.txt
# prepare database
# select related sequence in SILVA SSURef database
python scripts/HMWGS_stat.py hmwgs_list.txt silva_tax.txt > hmwgs_silva_list.txt
cat hmwgs_silva_list.txt | xargs samtools faidx SILVA_128_SSURef_Nr99_tax_silva.fasta > hmwgs_silva_subset.fa
cat hmwgs_silva_list.txt | xargs samtools faidx SILVA_128_aligned.fasta > hmwgs_silva_subset_aligned.fa
# build tax tree & extract position information
python scripts/extract_posision.py ecoli_ref.fa hmwgs_silva_subset_aligned.fa hmwgs_pos_info.pickle
python scripts/filter_tax_by_id.py hmwgs_silva_list.txt silva_tax.txt > hmwgs_tax.txt
python scripts/tax_node.py hmwgs_pos_info.pickle hmwgs_tax.txt hmwgs_tree.pickle hmwgs_tid.pickle
# calculate coverage
python scripts/usearch.py hmwgs_silva_subset.fa scripts/primers.fa hmwgs_result.txt
python scripts/convert_position.py hmwgs_pos_info.pickle hmwgs_result.txt > hmwgs_result_conv.txt
python scripts/binding_site.py hmwgs_result_conv.txt U00096.223771.225312 > hmwgs_site.txt
python scripts/calculate_coverage.py hmwgs_tree.pickle hmwgs_tid.pickle hmwgs_site.txt hmwgs_result_conv.txt > hmwgs_coverage.txt
python scripts/filter_primer.py --tax 'Bacteria,Archaea' scripts/primers.fa hmwgs_coverage.txt hmwgs_site.txt > hmwgs_pair.txt
python scripts/primer_pair.py hmwgs_tree.pickle hmwgs_tid.pickle hmwgs_site.txt hmwgs_pair.txt hmwgs_result_conv.txt > hmwgs_pair_converage.txt
