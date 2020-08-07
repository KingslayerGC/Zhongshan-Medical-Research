# %%
import os
import pandas as pd
import numpy as np
os.chdir(r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\GWAS")
df_list=[]
for parents, dirname, filename in os.walk(r"C:\Users\Mac\Desktop\过程\项目\医统\GWAS\GWAS"):
    for file in filename:
        df = pd.read_excel(os.path.join(parents, file))
        df_list.append(df)

def clean(table):
    for i in range(len(table)):
        try:
            float(table.iloc[i])
        except ValueError:
            table.iloc[i] = np.nan

df = df_list[0]
df = df.rename(columns={'CHR_ID':'Chr', 'CHR_POS':'Pos_hg19', 'SNPS':'SNP'})
clean(df['Pos_hg19'])

df = df[['Chr','Pos_hg19','SNP']]
df = df.dropna(subset=df.columns)

df = df.iloc[0:50, :]

df['SNP'].is_unique
df.info()

df.to_csv(r"C:\Users\Mac\Desktop\new_snp.csv", sep='\t', index=False)
