# 用卡方检验和t检验筛选变量，导入temp1已经生成的矩阵（含有所有哑变量和用相关性筛选的其他变量）
# %%
import numpy as np
import pandas as pd
import sklearn.feature_selection as sf
import scipy.stats as ss
# %%
ori_observation_df=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\data for ehr\training\observation.csv',header=0,index_col=0,sep=',',low_memory=False)
ori_drug_df=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\data for ehr\training\drug_exposure.csv',header=0,index_col=0,sep=',',low_memory=False)
ori_measurement_df=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\data for ehr\training\measurement.csv',header=0,index_col=0,sep=',',low_memory=False)
# 以下两个文件过大，故采用迭代读法，每次读取100万行，储存在列表中
procedure_reader=list(pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\data for ehr\training\procedure_occurrence.csv',header=None,sep=',',low_memory=False,chunksize=1000000))
condition_reader=list(pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\data for ehr\training\condition_occurrence.csv',header=None,sep=',',low_memory=False,chunksize=1000000))
# %%
total_matrix=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\total_matrix_2.csv',\
                         header=0,index_col=0,sep=',',low_memory=False)
death_df=total_matrix['death']
person_list=list(death_df.index)
# %%
dummy_variable_df=total_matrix[['gender_concept_id', 'race_concept_id', 'ethnicity_concept_id']]
select = sf.SelectKBest(sf.chi2,k='all')
select.fit_transform(dummy_variable_df,death_df)
#输出所有哑变量的pvalue
select.pvalues_
# %%
#total_pvalue_df=total_pvalue_df.fillna(0.9999)
#total_pvalue_df.to_csv(r'C:\Users\Mac\Desktop\过程\项目\total_pvalue.csv')
total_pvalue_df=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\total_pvalue.csv',header=0,index_col=0,sep=',',low_memory=False)
# %%
import math
loglist=[]
string_list=['observation','drug','procedure','condition','measurement']
total_pvalue_df=total_pvalue_df.sort_values(axis=1,by='p_value',ascending=True)
for i in range(len(list(total_pvalue_df.columns))):
    loglist.append(-math.log(total_pvalue_df.iloc[0,i],10))
    if total_pvalue_df.iloc[0,i] > 0.001:
        break
varient_list=list(total_pvalue_df.columns)[0:i]
num_dict={'observation':0,'drug':0,'procedure':0,'condition':0,'measurement':0}
for term in varient_list:
    if term[0]=='o':
        num_dict['observation']+=1
    if term[0]=='d':
        num_dict['drug']+=1
    if term[0]=='p':
        num_dict['procedure']+=1 
    if term[0]=='c':
        num_dict['condition']+=1
    if term[0]=='m':
        num_dict['measurement']+=1
# %%
# 处理数量不大的几种变量
total_pvalue_df=pd.DataFrame(index=['p_value'])
other_variable_df=total_matrix[['visit_id_0', 'visit_id_9201', 'visit_id_9202','age','death']]
'''
# 以下计算pvalue
X1=other_variable_df[other_variable_df['death']==0].drop(['death'],axis=1)
X2=other_variable_df[other_variable_df['death']==1].drop(['death'],axis=1)
[s,pvalue] = ss.ttest_ind(X1,X2)
total_pvalue_df=pd.DataFrame(pvalue,columns=['p_value'],index=X1.columns).T
'''
other_variable_df=other_variable_df[[v for v in list(other_variable_df.columns)\
                                     if v in varient_list]]
# %%
observation_df=ori_observation_df.copy()
# 将person_id一列作为行索引
observation_df.set_index(['person_id'], inplace=True)
# 仅保留有用变量
observation_df['count']=1
observation_df=observation_df[['observation_concept_id','count']]
# 将特征量拆分成多个子特征量，对每个样本每个子量单独计数作为该子量的值
observation_df['count']=observation_df['count'].astype('float')
# 依据drugid—personid进行分类
grouped=observation_df.groupby(['observation_concept_id','person_id']).sum()
'''
# 以下生成p_value的dataframe
for concept_id in grouped.index.levels[0]:
    df=grouped.loc[concept_id].join(death_df, how='right')
    df=df.fillna(0.0)
    X1=df[df['death']==0].drop(['death'],axis=1)
    X2=df[df['death']==1].drop(['death'],axis=1)
    [s,pv]=ss.ttest_ind(X1,X2)
    total_pvalue_df['observation_id_'+str(concept_id)]=pv
'''
df=pd.DataFrame(columns=[],index=person_list)
for concept_id in grouped.index.levels[0]:
    if 'observation_id_'+str(concept_id) not in varient_list:
        continue
    # 生成每个子特征量
    df1=grouped.loc[concept_id]
    df=df1.join(df, how='right')
    # 重命名子特征量
    df.rename(columns={'count':'observation_id_'+str(concept_id)},inplace=True)
observation_df=df
# %%
drug_df=ori_drug_df.copy()
# 将person_id一列作为行索引
drug_df.set_index(['person_id'], inplace=True)
# 仅保留有用变量
drug_df=drug_df[['drug_concept_id','quantity']]
# 将特征量拆分成多个子特征量，对每个样本每个子量单独计数作为该子量的值
drug_df['quantity']=drug_df['quantity'].astype('float')
# 依据drugid—personid进行分类
grouped=drug_df.groupby(['drug_concept_id','person_id']).sum()
'''
# 以下生成p_value的dataframe
for concept_id in grouped.index.levels[0]:
    df=grouped.loc[concept_id].join(death_df, how='right')
    df=df.fillna(0.0)
    X1=df[df['death']==0].drop(['death'],axis=1)
    X2=df[df['death']==1].drop(['death'],axis=1)
    [s,pv]=ss.ttest_ind(X1,X2)
    total_pvalue_df['drug_id_'+str(concept_id)]=pv
'''
df=pd.DataFrame(columns=[],index=person_list)
for concept_id in grouped.index.levels[0]:
    if 'drug_id_'+str(concept_id) not in varient_list:
        continue
    # 生成每个子特征量
    df1=grouped.loc[concept_id]
    df=df1.join(df, how='right')
    # 重命名子特征量
    df.rename(columns={'quantity':'drug_id_'+str(concept_id)},inplace=True)
drug_df=df
# %%
measurement_df=ori_measurement_df.copy()
# 将person_id一列作为行索引
measurement_df.set_index(['person_id'], inplace=True)
# 仅保留有用变量
measurement_df['count']=1
measurement_df=measurement_df[['measurement_concept_id','count']]
# 将特征量拆分成多个子特征量，对每个样本每个子量单独计数作为该子量的值
measurement_df['count']=measurement_df['count'].astype('float')
# 依据drugid—personid进行分类
grouped=measurement_df.groupby(['measurement_concept_id','person_id']).sum()
'''
# 以下生成p_value的dataframe
for concept_id in grouped.index.levels[0]:
    df=grouped.loc[concept_id].join(death_df, how='right')
    df=df.fillna(0.0)
    X1=df[df['death']==0].drop(['death'],axis=1)
    X2=df[df['death']==1].drop(['death'],axis=1)
    [s,pv]=ss.ttest_ind(X1,X2)
    total_pvalue_df['measurement_id_'+str(concept_id)]=pv
'''
df=pd.DataFrame(columns=[],index=person_list)
for concept_id in grouped.index.levels[0]:
    if 'measurement_id_'+str(concept_id) not in varient_list:
        continue
    # 生成每个子特征量
    df1=grouped.loc[concept_id]
    df=df1.join(df, how='right')
    # 重命名子特征量
    df.rename(columns={'count':'measurement_id_'+str(concept_id)},inplace=True)
measurement_df=df
# %%
k=0
flag=1
condition_df_dict={}
# 循环地处理列表中的每个df
for chunk in condition_reader:
    condition_df=chunk.copy()
    # 将第一个表的第一行作为列索引
    if flag==1:
        cond_index=list(condition_df.iloc[0,:])
        condition_df.drop(condition_df.index[0],inplace=True)
        flag=0
    condition_df.columns=cond_index
    # 将person_id一列作为行索引
    condition_df.set_index(['person_id'], inplace=True)
    # 仅保留有用变量
    condition_df['count']=1
    condition_df=condition_df[['condition_concept_id','count']]
    # 将特征量拆分成多个子特征量，对每个样本每个子量单独计数作为该子量的值
    grouped=condition_df.groupby(['condition_concept_id','person_id']).sum()
    for concept_id in grouped.index.levels[0]:
        if 'condition_id_'+str(concept_id) not in varient_list:
            continue
        df1=grouped.loc[concept_id]
        df1.rename(columns={'count':'condition_id_'+str(concept_id)},inplace=True)
        if 'condition_id_'+str(concept_id) in list(condition_df_dict.keys()):
            condition_df_dict['condition_id_'+str(concept_id)]=condition_df_dict\
            ['condition_id_'+str(concept_id)].add(df1,fill_value=0)
        else:
            condition_df_dict['condition_id_'+str(concept_id)]=df1
    k=k+1
    print(k)
condition_df=pd.concat(list(condition_df_dict.values()),axis=1,join='outer')
condition_df=condition_df.join(death_df, how='right').drop(['death'],axis=1)
'''
    # 以下生成相关系数dataframe
    for concept_id in grouped.index.levels[0]:
        df1=grouped.loc[concept_id]
        df1.rename(columns={'count':'condition_id_'+str(concept_id)},inplace=True)
        if 'condition_id_'+str(concept_id) in list(condition_df_dict.keys()):
            condition_df_dict['condition_id_'+str(concept_id)]=condition_df_dict\
            ['condition_id_'+str(concept_id)].add(df1,fill_value=0)
        else:
            condition_df_dict['condition_id_'+str(concept_id)]=df1
    flag=0
    k+=1
    print(k)
i=0
for (k,v) in condition_df_dict.items():
    df=v.join(death_df, how='right')
    df=df.fillna(0.0)
    X1=df[df['death']==0].drop(['death'],axis=1)
    X2=df[df['death']==1].drop(['death'],axis=1)
    [s,pv]=ss.ttest_ind(X1,X2)
    total_pvalue_df[k]=pv
    i+=1
    if i%500==0:
        print(i)
'''
# %%
k=0
flag=1
procedure_df_dict={}
# 循环地处理列表中的每个df
for chunk in procedure_reader:
    procedure_df=chunk.copy()
    # 将第一个表的第一行作为列索引
    if flag==1:
        cond_index=list(procedure_df.iloc[0,:])
        procedure_df.drop(procedure_df.index[0],inplace=True)
        flag=0
    procedure_df.columns=cond_index
    # 将person_id一列作为行索引
    procedure_df.set_index(['person_id'], inplace=True)
    # 仅保留有用变量
    procedure_df['count']=1
    procedure_df=procedure_df[['procedure_concept_id','count']]
    # 将特征量拆分成多个子特征量，对每个样本每个子量单独计数作为该子量的值
    grouped=procedure_df.groupby(['procedure_concept_id','person_id']).sum()
    for concept_id in grouped.index.levels[0]:
        if 'procedure_id_'+str(concept_id) not in varient_list:
            continue
        df1=grouped.loc[concept_id]
        df1.rename(columns={'count':'procedure_id_'+str(concept_id)},inplace=True)
        if 'procedure_id_'+str(concept_id) in list(procedure_df_dict.keys()):
            procedure_df_dict['procedure_id_'+str(concept_id)]=procedure_df_dict\
            ['procedure_id_'+str(concept_id)].add(df1,fill_value=0)
        else:
            procedure_df_dict['procedure_id_'+str(concept_id)]=df1
    k=k+1
    print(k)
procedure_df=pd.concat(list(procedure_df_dict.values()),axis=1,join='outer')
procedure_df=procedure_df.join(death_df, how='right').drop(['death'],axis=1)
'''
    # 以下生成相关系数dataframe
    for concept_id in grouped.index.levels[0]:
        df1=grouped.loc[concept_id]
        df1.rename(columns={'count':'procedure_id_'+str(concept_id)},inplace=True)
        if 'procedure_id_'+str(concept_id) in list(procedure_df_dict.keys()):
            procedure_df_dict['procedure_id_'+str(concept_id)]=procedure_df_dict\
            ['procedure_id_'+str(concept_id)].add(df1,fill_value=0)
        else:
            procedure_df_dict['procedure_id_'+str(concept_id)]=df1
    flag=0
    k+=1
    print(k)
i=0
for (k,v) in procedure_df_dict.items():
    df=v.join(death_df, how='right')
    df=df.fillna(0.0)
    X1=df[df['death']==0].drop(['death'],axis=1)
    X2=df[df['death']==1].drop(['death'],axis=1)
    [s,pv]=ss.ttest_ind(X1,X2)
    total_pvalue_df[k]=pv
    i+=1
    if i%500==0:
        print(i)
'''        
#%%
total_df=pd.concat([dummy_variable_df,death_df,other_variable_df,observation_df,drug_df,\
                        measurement_df,condition_df,procedure_df],axis=1,join='outer')
total_df=total_df.fillna(0).astype('float')
total_df.index=total_df.index.astype('int')
total_df.sort_index(inplace=True)
#total_df.to_csv(r'C:\Users\Mac\Desktop\过程\项目\total_matrix3.csv')
