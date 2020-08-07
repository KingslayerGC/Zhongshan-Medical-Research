# %%
args = (r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\示例代码及结果\iRIGS_code\SNP_file\SCZ_108_loci",
        1000000, r"C:\Users\Mac\Desktop\result", "SCZ")

import numpy as np
import pandas as pd

## file loading function
def csv_read(path):
    file = pd.read_csv(path, header=None)
    file = file[0].str.split(pat=None, n=0, expand=True)
    file.columns= file.loc[0]
    file = file.drop(axis=0, index=[0]).set_index(np.arange(file.shape[0]-1))
    return file

## match and return the index
def match(a, b):
    index=[]
    for g in b:
        try:
            ind = list(b).index(g)
        except ValueError:
            ind = np.nan
        index.append(ind)
    return index

# %%
# =============================================================================
# extract_candidate_genes
# =============================================================================
# %%
## load the data
def extract_candidate_genes(args):
    gwas_file = args[0]
    flanking = args[1]
    gwas = csv_read(gwas_file)
    
    ## check column names and flanking
    col = pd.Index(['SNP','Chr','Pos_hg19'])
    if not col.isin(gwas.columns).all():
        raise('Wrong column names!')
        if flanking <0:
            raise('Please verify the flank region number!')
            
    ## add prefix "chr" in column "Chr"
    gwas['Chr'] = gwas['Chr'].apply(lambda x : 'chr'+x)        
    gene_path = r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\示例代码及结果\iRIGS_code\supporting_files\All_human_genes"
    gene = csv_read(gene_path)
    
    ## 筛出备选基因，并指定到对应SNP
    gwas['Pos_hg19'] = gwas['Pos_hg19'].astype(int)
    gene[['start_hg19', 'end_hg19']] = gene[['start_hg19', 'end_hg19']].astype(int, copy=False)
    output = pd.DataFrame([])
    for i in range(gwas.shape[0]):
        start = max(gwas.loc[i,'Pos_hg19']-flanking, 0)
        end = gwas.loc[i,'Pos_hg19'] + flanking
        temp = gene.loc[gene['chrom']==gwas.loc[i,'Chr']]
        index1 = (temp['start_hg19']<end) & (temp['start_hg19']>start)
        index2 = (temp['end_hg19']<end) & (temp['end_hg19']>start)
        temp = temp.loc[index1|index2, :]
        if temp.shape[0] >0:
            temp['SNP'] = gwas.loc[i, 'SNP']
            temp['SNP_chr'] = gwas.loc[i, 'Chr']
            temp['SNP_pos_hg19'] = gwas.loc[i, 'Pos_hg19']
            output = pd.concat([output, temp], axis=0)

    return output

# %%
# =============================================================================
# extra_evi
# =============================================================================
def extra_evi(gene):
    gene['Name'] = gene['Name'].apply(lambda x:x[0:15])
    
    ## distance to TSS
    def compute_dist(gene):
        index = gene['strand']=='+'
        tss_dist = pd.Series(np.arange(gene.shape[0]), index=gene.index)
        tss_dist[index] = abs(gene.loc[index]['start_hg19']-gene.loc[index]['SNP_pos_hg19'])
        tss_dist[~index] = abs(gene.loc[~index]['end_hg19']-gene.loc[~index]['SNP_pos_hg19'])
        gene['dist'] = tss_dist
        return gene
    gene = compute_dist(gene)
    
    ## adding capture Hi-C
    print("Collecting extra evidence: capture Hi-C...")
    def adding_caphic(gene):
        path = r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\示例代码及结果\iRIGS_code\supporting_files/capHiC/"
        cap4 = path + "GM12878_DRE_number"
        cap4 = csv_read(cap4)
        index = match(gene['Name'],cap4['Name'])
        gene['cap4_enhancer_no'] = np.array(cap4.loc[index]['cap4_enhancer_no'])
        return gene
    gene = adding_caphic(gene)

### still got something in this function


# %%
# =============================================================================
# sampling_procedure
# =============================================================================

## maximum sampling round allowed at burn and post burn in step
burnin_round = 3000
after_burnin_round = 3000
exclude_samegene = True

## specify different network
print("Loading propagation probabilities...")
pro_p = pd.read_pickle(r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\示例代码及结果\iRIGS_code\supporting_files\All_human_genes.pkl",)
pro_p.drop(columns=['Unnamed: 0'], axis=1, inplace=True)
pro_p.index = pro_p.columns

nodes = list(pro_p.columns)
gene = extract_candidate_genes(args)
gene = gene.loc[~(gene['official_name'] == 'NA')]
gene = gene.loc[gene['official_name'].isin(nodes)]
region = []
for snp in gene['SNP'].unique():
    region.append(gene.loc[gene['SNP']==snp].set_index(['SNP'])['official_name'])
print(len(gene['official_name'].unique()), "genes from",
      len(region), "loci were found with propagation probability")

pro_p = pro_p.loc[:, pro_p.columns.isin(gene['official_name'].unique())]

## burn in step
thres = 0.01
pickup = 0
num_region = len(region)
circle = -1
chosen = None
remaining = pd.concat([loci.sample(n=1, replace=True) for loci in region])
num0 = np.zeros(gene.shape[0])

dif = thres + 1
dif_record = None

while (dif>thres) & (circle<(burnin_round+1)):
    pickup = pickup %num_region +1
    if pickup == 1:
        if not chosen == None:
            num1 = None
            for j in range(len(rigion)):
                num1 = 
    
    pickup_p = pro_p.loc[: , pro_p.columns.isin(
        remaining.drop(remaining.index[[pickup-1]])
        )]
    pickup_p = pickup_p.loc[pickup_p.index.isin(region[pickup-1])]
    
    if (exclude_samegene) & (sum(pickup_p.shape)!=0):
        pickup_p = pickup_p.loc[~pickup_p.columns.isin(pickup_p.index)]
        if (sum(pickup_p.shape)!=0) & (pickup_p.shape[1]==0)
            raise("Error: no conditional genes!")
    if sum(pickup_p.shape)==0:
        pickup_p = 1
        pickup_p.columns = region[pickup-1]
    else:
        pickup_p = pickup_p.sum(axis=1)
        pickup_p = extra_weight.loc[
            match(pickup_p.columns, extra_weight['gene'])
            ]['extra_weight'] *pickup_p

    if sum(pickup_p)==0:
        for i in range(len(pickup_p)):
            pickup_p.iloc[i] = 1 / len(pickup_p)
    
    remaining.iloc[pickup-1] = pickup_p.columns.sample(
        n=1, replace=True, weights=pickup_p)
    chosen = pd.concat([chosen,remaining], axis=1)








