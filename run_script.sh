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
python scripts/filter_primer.py result.txt temp_binding_site.txt > primer_pair.txt
# calculate primer pair coverage
