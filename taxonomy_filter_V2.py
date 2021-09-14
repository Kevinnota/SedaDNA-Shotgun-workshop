import pandas as pd
import re
from tqdm import tqdm
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import textwrap
parser = argparse.ArgumentParser(description='Input and output')

parser.add_argument('-in', '--input', 
    help="Path to the input file", type=str)
parser.add_argument('-o', '--output', 
    help="Path and basename of the output files", type=str)
parser.add_argument('-t', '--taxonomy', type=str,
    help=textwrap.dedent('''
    Path to the ncbi taxonomy dump files: get by
    |$
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    |$
    tar -zxvf taxdump.tar.gz '''))
parser.add_argument('-r', '--reads', 
    help="Path to the input reads that need to be filtered", type=str)
parser.add_argument('-ids', '--taxids', 
    help="Taxids of the highest taxonomic level to keep | such as for plants use 33090 (Viridiplantae)", type=str)
  
args = parser.parse_args()
print(args.input)
print(args.output)

print("\n")

PIA_tabel = pd.read_table(args.input,  header=None, delimiter=",")
ncbi_taxonomy_nodes = pd.read_table(args.taxonomy+"/nodes.dmp",  header=None, delimiter="\t")
#ncbi_taxonomy_merged = pd.read_table("merged.dmp",  header=None, delimiter="\t")
ncbi_taxonomy_names = pd.read_table("names.dmp",  header=None, delimiter="\t")
filtered_reads = pd.DataFrame()

print("\n")
print("Overview of which taxa will be kept")
taxids=(re.split(",", args.taxids))
taxids=[int(x) for x in taxids]
names_to_keep = []
for i in range(len(taxids)) :
    names = ncbi_taxonomy_names.loc[taxids[i]==ncbi_taxonomy_names.loc[:,0]]
    scientific_name =  names.loc["scientific name"==names.loc[:,6],2].item()
    print("Keep ",scientific_name)
    names_to_keep.append(scientific_name)

print("\n")

#Viridiplantae 33090
#Actinopterygii 7898

print("Matching taxonomy to the taxa of interest")
for i in tqdm(range(len(PIA_tabel))): #there might be a problem here when the taxid has been merged with another taxa. In this case the loop needs to be modified so that if the taxid is not found it goes through the ncbi_taxonomy_merged 
    parent_id=PIA_tabel.loc[i, 1]
    level_reached=False
    if(parent_id in taxids):
        level_reached=True    
    while not (level_reached==True):
        parent_id=ncbi_taxonomy_nodes.loc[(parent_id==ncbi_taxonomy_nodes.loc[:,0]),2].item()
        if(parent_id in taxids):
            level_reached=True
            for x in range(len(taxids)):
                if(taxids[x]==parent_id):
                    taxon_name=names_to_keep[x]
            data = {'read':[PIA_tabel.loc[i,0]],
                'original_taxid':[PIA_tabel.loc[i,1]],
                'taxon_group':[taxon_name]}
            filtered_reads=filtered_reads.append(pd.DataFrame(data), ignore_index=True)
            
        if(parent_id==1):
            level_reached=True

print("\n\nWriting fasta file with only sequences of interest")

reads = list(SeqIO.parse(args.reads, "fasta"))
new_fasta = []
for i in tqdm(range(len(filtered_reads))):
    r=0
    found_read=False
    while (found_read!=True):
        if(filtered_reads.loc[i,"read"]==reads[r].id):
            new_fasta.append(reads[r])
            found_read=True
        r=r+1
    
print()

SeqIO.write(new_fasta, args.output+".fasta", "fasta")
filtered_reads.to_csv(args.output+".csv", header=True, index=False, sep=',')
            
print("Done")