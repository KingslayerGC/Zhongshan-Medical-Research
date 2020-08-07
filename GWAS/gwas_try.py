# %%
import numpy as np
import pandas as pd
dataset = pd.read_excel(r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\GWAS\gwas-bipolar disorder.xlsx")

# %%
data = dataset.copy()
data = data[data['P-VALUE']<10**-8]

# %%
gene_df = data['SNP_GENE_IDS'].str.split(pat=',', expand=True)
gene_df = pd.concat([gene_df[col] for col in list(gene_df.columns)])
gene_df = gene_df[~gene_df.isnull()].to_frame()
gene_df['SNPS'] = data.loc[gene_df.index, 'SNPS']
gene_df[0] = [item.lstrip() for item in list(gene_df[0])]
#gene_df['count'] = 1
#gene_df = gene_df.groupby(by=[0,'SNP']).count()

# %%
nearby_gene_df = data.dropna(subset=['UPSTREAM_GENE_ID'])
nearby_gene_df = nearby_gene_df[nearby_gene_df['UPSTREAM_GENE_DISTANCE']<2000]

nearby_gene_df = data.dropna(subset=['DOWNSTREAM_GENE_ID'])
nearby_gene_df = nearby_gene_df[nearby_gene_df['DOWNSTREAM_GENE_DISTANCE']<2000]
nearby_gene_df = nearby_gene_df[['DOWNSTREAM_GENE_ID', 'SNPS']]














