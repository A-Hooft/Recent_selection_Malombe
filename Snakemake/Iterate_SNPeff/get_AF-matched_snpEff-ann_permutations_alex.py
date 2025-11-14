import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys 
sys.path.append("/data/antwerpen/grp/asvardal/hs_tools")
from pypopgen3.modules import vcfpandas as vp
import timeit

# Parameters - Define before running
ana_dir = '/scratch/antwerpen/grp/asvardal/projects/cichlid/analyses/Copadichromis_FIE_Alex'

allsnp_alex_fn = snakemake.input[0]
outlierspersnp_alex_fn = snakemake.input[1]
anno_vcf_fn_gz = snakemake.input[2]

snps_regions_fn = snakemake.input[3]

annotations = ['3_prime_UTR_variant',
 '5_prime_UTR_variant',
 'downstream_gene_variant',
 'intergenic_region',
 'intron_variant',
 'missense_variant',
 'missense_variant&splice_region_variant',
 'splice_acceptor_variant&intron_variant',
 'splice_region_variant&intron_variant',
 'splice_region_variant&synonymous_variant',
 'stop_gained',
 'synonymous_variant',
 'upstream_gene_variant']

n_iterations = 10


############ Read in data ###############
# MAF data
outlierspersnp_alex = pd.read_csv(outlierspersnp_alex_fn,
                          sep = '\t',
                          names=['chr','ps','maf'])


allsnp_alex = pd.read_csv(allsnp_alex_fn,
                     sep = '\t',
                     names = ['chr','ps','maf'])    # All SNP present in the vcf used for the GWA
############ Functions ###############

# Bin MAF (0.05 bins)
# This way, MAF are grouped according to their value into different categories
def get_maf_bin(maf):
    return round(maf/0.05)

########### Run ####################


# Apply maf_bin function to datasets (select)
allsnp_alex['mafbin'] = allsnp_alex['maf'].apply(get_maf_bin)
outlierspersnp_alex['mafbin'] = outlierspersnp_alex['maf'].apply(get_maf_bin)

# Create a groupby object
mafgroup = allsnp_alex.groupby('mafbin')

# Get number of SNP per bin in test SNP dataset (select)
# 0.01%:
snp_per_bin = outlierspersnp_alex.groupby('mafbin').apply(len)

########## Run ################
start = timeit.default_timer()

ann_per_pos = {}
for annotat in annotations:
    # Lists of effects:
    ann_per_pos[annotat] = []

for i in range(n_iterations):
    # Initialize a dictionary to store the counts
    annotation_counts = {annotat: 0 for annotat in annotations}
    print(f'iteration {i}')
    
    for mafbin,nsnps in snp_per_bin.items():
        df = mafgroup.get_group(mafbin)      
        target_snps = df.sample(nsnps)
        
        for chrom,pos in target_snps[['chr','ps']].itertuples(index=False):
            header, info_dic = vp.parse_vcf_header(anno_vcf_fn_gz.format(chrom))
            a = vp.get_vcf_df(anno_vcf_fn_gz.format(chrom), chrom=str(chrom), start=pos ,end=pos , header=header)
            # Split the 'INFO' string
            info_string = a.iloc[0]['INFO'].split('ANN')[1]
            # Iterate over the annotations and count occurrences
            for annotat in annotations:
                if annotat in info_string:
                    annotation_counts[annotat] += 1            
    
    for annotat in annotations:
        # Lists of effects:
        ann_per_pos[annotat].append(annotation_counts[annotat])
    
# wf = open(f"/scratch/antwerpen/grp/asvardal/projects/cichlid/analyses/Copadichromis_FIE_Alex/_data/Genomic_features/eff-per-pos_AF-matched_SNPs_{n_iterations}_iterations.txt","w")
wf = open(snakemake.output[0],"w")
wf.write(str(ann_per_pos))
wf.close()

stop = timeit.default_timer()
print('Time: ', stop - start) 
