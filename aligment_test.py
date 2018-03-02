# coding: utf-8
from Bio import SeqIO
ecoli = SeqIO.read('ecoli_ref.fa', 'fasta')
ecoli
ecoli.seq
[b for b in ecoli.seq]
ecoli.seq[13132]
ecoli.seq[13132 + 1]
ecoli.seq[13132:13132+10]
ecoli.seq[13132:13132+100]
ecoli
ecoli.seq.ungap()
ecoli.seq.ungap('.')
ecoli.seq.ungap('.').ungap('-')
len(ecoli.seq.ungap('.').ungap('-'))
ecoli_ungap = ecoli.seq.ungap('.').ungap('-')
ecoli_ungap[521]
ecoli_ungap[520]
ecoli_ungap[520:520+20]
ecoli.seq.find('GCAGCCGCGGUAAUACGGAG')
ecoli_ungap.find('GCAGCCGCGGUAAUACGGAG')
s = 13132
x = 0
while s < 50000:
    if ecoli[s] == 'GCAGCCGCGGUAAUACGGAG'[x]:
        x += 1
        print(s)
    s += 1
    
x = 0
while s < 50000 and x < 20:
    if ecoli[s] == 'GCAGCCGCGGUAAUACGGAG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 13132
x = 0
while s < 50000 and x < 20:
    if ecoli[s] == 'GCAGCCGCGGUAAUACGGAG'[x]:
        x += 1
        print(s)
    s += 1
    
ecoli[13132:13878]
ecoli[13132:13879]
ecoli.seq[13132:13879]
ecoli.seq[13132:13879].ungap('-')
ecoli.seq[13132:13879].ungap('-').__len__()
ecoli.seq[13132:13879]
print(ecoli.seq[13132:13879])
from pyfaidx import Fasta
gene = Fasta('./SILVA_128_aligned.fa')
gene['JQ467910.1.1405']
gene['JQ467910.1.1405'].seq
gene['JQ467910.1.1405'][:]
gene['JQ467910.1.1405'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
gene['JQ467910.1.1405'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
gene['JQ467910.1.1405'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
from Bio.Seq import Seq
JQ467910 = Seq(gene['JQ467910.1.1405'][13132:13879])
JQ467910 = Seq(gene['JQ467910.1.1405'][13132:13879].seq)
JQ467910
JQ467910.ungap('-')
len(JQ467910.ungap('-'))
JQ467910 = Seq(gene['JQ467910.1.1405'][:].seq)
len(JQ467910.ungap('-'))
JQ467910
len(JQ467910.ungap('-').ungap('.'))
JQ467910 = JQ467910.ungap('-').ungap('.')
JQ467910
JQ467910[519:539]
len(JQ467910[519:539])
gene['AY592675.1.1286'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
gene['DQ159175.1.1521'][13132:13879]
gene['JQ460811.1.1389'][13132:13879]
gene['JQ460811.1.1389'][13132:13879].replace('-', '')
gene['JQ460811.1.1389'][13132:13879]
gene['JQ460811.1.1389'][13132:13879]
gene['JQ193480.1.1365'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
gene['U00096.223771.225312'][13132:13879].seq
gene['U00096.223771.225312'][13132:13879].seq == gene['JQ193480.1.1365'][13132:13879].seq
test = Seq(gene['JQ193480.1.1365'][:].seq).ungap('.').ungap('-')
len(test)
test[497:517]
len(test[497:517])
gene['JF428926.1.1482'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
Seq(gene['U00096.223771.225312'][13132:13879].seq).ungap('-')
Seq(gene['JF428926.1.1482'][13132:13879].seq).ungap('-')
Seq(gene['JF428926.1.1482'][13132:13879].seq).ungap('-').__len__
Seq(gene['JF428926.1.1482'][13132:13879].seq).ungap('-').__len__()
s = 13132
x = 0
while s < 50000 and x < 20:
    if Seq(gene['JF428926.1.1482'][:].seq) == 'GCAGCCGCGGUAAUACGGAG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 13132
x = 0
temp = Seq(gene['JF428926.1.1482'][:].seq)
while s < 50000 and x < 20:
    if temp == 'GCAGCCGCGGUAAUACGGAG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 13132
x = 0
temp = Seq(gene['JF428926.1.1482'][:].seq)
while s < 50000 and x < 20:
    if temp == 'GCAGCCGCGGUAAUACGAG'[x]:
        x += 1
        print(s)
    s += 1
    
Seq(gene['JF428926.1.1482'][:].seq)
len(Seq(gene['JF428926.1.1482'][:].seq))
Seq(gene['JF428926.1.1482'][13132:13879].seq).ungap('-')
Seq(gene['JF428926.1.1482'][13132:13879].seq)
gene['JF428926.1.1482'][13132:13879].seq
gene['U00096.223771.225312'][13132:13879].seq
gene['U00096.223771.225312'][13132:13879]
gene['U00096.223771.225312'][13132:13879+2]
gene['U00096.223771.225312'][13132:13879]
gene['JF428926.1.1482'][13132:13879+1]
gene['JF428926.1.1482'][13132:13879+2]
gene['U00096.223771.225312'][13132:13879]
gene['U00096.223771.225312'][13132:13879+2]
gene['JF428926.1.1482'][13132:13879+2]
Seq(gene['JF428926.1.1482'][13132:13879+2].seq).ungap('-')
Seq(gene['JF428926.1.1482'][13132:13879+2].seq).ungap('-').__len__()
Seq(gene['U00096.223771.225312'][13132:13879].seq).ungap('-')
Seq(gene['HG975445.39219252.39220730'][13132:13879].seq).ungap('-')
gene['HG975445.39219252.39220730'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
gene['U00096.223771.225312'][:]
gene['U00096.223771.225312'][:].seq
str(gene['U00096.223771.225312'][:])
str(gene['U00096.223771.225312'][:]).replace('.', '')
str(gene['U00096.223771.225312'][:]).replace('.', '').replace('-', '')
ecoli_ungap
ecoli_ungap.find('GUCGUCAGCUCGUGUUG')
s = 0
x = 0
temp = ecoli
while s < 50000 and x < 20:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 0
x = 0
temp = ecoli.seq
while s < 50000 and x < 20:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
ecoli.seq
ecoli.seq[3]
s = 0
x = 0
temp = ecoli.seq
while s < 50000 and x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
ecoli[3]
ecoli[50000]
'GUCGUCAGCUCGUGUUG'[2]
s = 0
x = 0
temp = ecoli.seq
while s < 50000 or x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 0
x = 0
temp = ecoli
while s < 50000 or x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 0
x = 0
temp = str(ecoli.seq)
while s < 50000 or x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 0
x = 0
temp = str(ecoli.seq)
while s < 50000 and x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
s = 30000
x = 0
temp = str(ecoli.seq)
while s < 50000 and x < 17:
    if temp == 'GUCGUCAGCUCGUGUUG'[x]:
        x += 1
        print(s)
    s += 1
    
ref_seq = str(SeqIO.read('ecoli_ref.fa', 'fasta').seq)
pos = -1
ref_pos = []
for base in ref_seq:
    if base != '.' and base != '-':
        pos += 1
    ref_pos.append(pos)
ref_pos
ref_pos[:-10]
ref_pos[1000]
ref_pos[:-100]
pos
ref_pos[-1]
gene['JQ467910.1.1405'][13132:13879]
gene['U00096.223771.225312'][13132:13879]
