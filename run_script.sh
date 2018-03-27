# download primers from probebase
python scripts/get_primers.py > scripts/primers.txt
# make fai index
samtools faidx SILVA_128_aligned.fasta
# extract reference sequence
python scripts/get_reference.py SILVA_128_aligned.fasta U00096.223771.225312 > ecoli_ref.fa
# extract position information
python scripts/extract_posision.py ecoli_ref.fa SILVA_128_aligned.fasta pos_info_new.pickle
# run usearch search_oligodb
python scripts/usearch.py SILVA_128_SSURef_Nr99_tax_silva.fasta scripts/primers.fa temp_result.txt
# convert unaligned position to aligned position
python scripts/convert_position.py pos_info_new.pickle temp_result.txt > temp_result_conv.txt
# get reference binding sites
python scripts/binding_site.py temp_result_conv.txt U00096.223771.225312 > temp_binding_site.txt
# build taxonomy tree
python scripts/tax_node.py pos_info_new.pickle silva_tax.txt tax_tree.pickle tax_tid.pickle
# calculate primer coverage
python scripts/calculate_coverage.py tax_tree.pickle tax_tid.pickle temp_binding_site.txt temp_result_conv.txt > result.txt
# filter primer (coverage > 75%)
python scripts/filter_primer.py scripts/primers.fa result.txt temp_binding_site.txt > primer_pair.txt
# calculate primer pair coverage
python scripts/primer_pair.py tax_tree.pickle tax_tid.pickle temp_binding_site.txt primer_pair.txt temp_result_conv.txt > pair_result.txt
# calculate primer coverage (DEGEPRIME)
python scripts/filter_primer.py bact_top100.fa bact_top100_coverage.txt bact_top100_site.txt > bact_top100_pair.txt
python scripts/primer_pair.py tax_tree.pickle tax_tid.pickle bact_top100_site.txt bact_top100.fa bact_top100_result_conv.txt > bact_top100_pair_coverage.txt

usearch -usearch_global rep_set_v35.fna -db test.1.udb -id 0.9 -strand both -blast6out HMP_usearch.1.b6
usearch -usearch_global rep_set_v35.fna -db test.2.udb -id 0.9 -strand both -blast6out HMP_usearch.2.b6
python scripts/merge_usearch.py HMP_usearch.*.txt > HMP_usearch.txt
python scripts/HMQCP_stat.py HMP_usearch.txt temp_result_conv.txt primer_pair.txt temp_binding_site.txt > HMQCP_converage.txt


grep 'g__' hmp1-II_metaphlan2-mtd-qcd.pcl | grep -v 's__' | egrep -v '(k__Viruses|k__Eukaryota)' | cut -f 1 > hmpwgs_list.txt