# %%
# =============================================================================
#  基础模板与参数设置
# =============================================================================
import os
import re
import pandas as pd
import numpy as np

# 设置pandas输出格式
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# 设置文件路径
dir_path = r"C:\Users\Mac\Desktop\original data\training"
useful_path = r"C:\Users\Mac\Desktop\original data\OMOP_usefule_columns.csv"
pkl_path = r"C:\Users\Mac\Desktop\过程\项目\医统\EHR challenge\临时数据-pickle"

# 仅考虑最后一次visit前m个月，后n个月内(下称目标周期)的数据
m = 6
n = 6

# 特征至少需要有 死亡人数 × percent 个非空值才视为有参考价值
percent = 0.05

# =============================================================================
# 定义操作函数
# =============================================================================
# 定义csv读取函数 
def csvread(path, usecols=None):
    return pd.read_csv(path, sep=',', header=0, index_col=0, usecols=usecols)

# 定义文件查找和读取函数
def loadcsv(name):
    df = csvread(useful_path)
    useful_col = df['ColNam'].loc[(df.primary==True) | (df.Useful==True)]
    for parents, _, filename in os.walk(dir_path):
        for file in filename:
            v = re.match(r'(\w*)\.\w*', file).group(1)
            # 若全部读取则不输入参数，建议逐个读取
            if (v==name) | (name==None):
                # 仅提取标记为primary或useful的列
                usecols = useful_col[useful_col.index.isin([v])]
                globals()[v + '_df'] = csvread(os.path.join(parents, file),
                                                  usecols=usecols)

# 滤除目标周期外的数据
from datetime import timedelta
def timefilter(df, column):   
    # 计算间隔
    df[column] = pd.to_datetime(df[column])
    last_visit_date = last_visit_df.loc[df['person_id']]
    last_visit_date.index = df.index
    gap_df = df[column] - last_visit_date
    # 计算目标间隔
    before = timedelta(30 * m)
    after = timedelta(30 * n)
    # 仅保留目标周期内的数据
    return df.loc[(-before < gap_df) & (gap_df < after)]

# 数据透视，丢弃过稀疏特征
def pivot(df, columns, prefix, values=None, drop=True):
    # 若没有计数列，则生成之
    if values==None:
        df['count'] = 1
        values = 'count'
    # 生成数据透视表
    pivot_df = pd.pivot_table(df, values=values, index='person_id',
                              columns=columns, aggfunc=np.sum)
    # 判断是否要进行丢弃过稀疏特征的步骤
    if drop==True:
        pivot_df = pivot_df.dropna(axis=1, thresh=nan_thresh)
    # 填充缺失值，转为int类型节省空间，然后给列名加上前缀
    return pivot_df.fillna(0).astype(int).add_prefix(prefix)

# 输出表格到pickle文件
def topkl(name):
    path = pkl_path + "\\" + name + ".pkl"
    globals()[name + '_df'].to_pickle(path)
    # 释放空间
    del globals()[name + '_df']

# %%
# =============================================================================
#  首先处理visit_occurrence_df
# =============================================================================
loadcsv('visit_occurrence')

# 查看数据类型和数据缺失情况
visit_occurrence_df.info(null_counts=True)

# 查看数据描述
visit_occurrence_df.astype(str).describe()

# visit_type_concept_id只有一种分类，visit_start_date明显不相关，放弃这些特征
visit_occurrence_df.drop(columns=['visit_type_concept_id', 'visit_start_date'],
                         axis=1, inplace=True)

# 将visit_end_date转为时间戳对象
visit_occurrence_df['visit_end_date'] = pd.to_datetime(
    visit_occurrence_df['visit_end_date'])

# 计算最后一次的日期（下面称目标周期)
last_visit_df = visit_occurrence_df\
    .sort_values(by=['visit_end_date'], ascending=False)\
        .drop_duplicates(subset=['person_id'], keep='first')\
            .set_index(keys='person_id').visit_end_date

# 滤除目标周期外的数据
visit_occurrence_df = timefilter(visit_occurrence_df, 'visit_end_date')

# 数据透视，丢弃过稀疏特征
visit_occurrence_df = pivot(visit_occurrence_df, 'visit_concept_id',
                            'visit_', drop=False)

# 合并
visit_occurrence_df = pd.concat([visit_occurrence_df, last_visit_df], axis=1)

# 储存为pickle格式
#topkl('visit_occurrence')

# %%
# =============================================================================
#  处理death_df
# =============================================================================
loadcsv('death')

# 查看数据类型和数据缺失情况
death_df.info()

# 查看数据描述
death_df.astype(str).describe()

# death_type_concept_id只有一种分类，放弃这个特征
death_df.drop(columns=['death_type_concept_id'], axis=1, inplace=True)

# 将死亡时间转为时间戳对象
death_df['death_date'] = pd.to_datetime(death_df['death_date'])

# 设置一个阈值，非零样本少于此值的特征将被认为没有参考意义
nan_thresh = len(death_df) * percent

# 储存为pickle格式
#topkl('death')

# %%
# =============================================================================
#  处理person_df
# =============================================================================
loadcsv('person')

# 查看数据类型和数据缺失情况
person_df.info()

# 查看数据描述
person_df = person_df.astype(str)
person_df.describe()

# time_of_death是空的，放弃这个特征
person_df.drop(columns=['time_of_birth'], axis=1, inplace=True)

# 用时间戳对象取代出生年月日
person_df['birth_date'] = person_df['year_of_birth'].map(str)\
    + '-' + person_df['month_of_birth'].map(str)\
        + '-' + person_df['day_of_birth'].map(str)
person_df.drop(columns=['year_of_birth', 'month_of_birth', 'day_of_birth'],
               axis=1, inplace=True)
person_df['birth_date'] = pd.to_datetime(person_df['birth_date'])

# 对分类变量单热编码
person_df = pd.get_dummies(person_df, {'gender_concept_id':'gender',
                                       'race_concept_id':'race',
                                       'ethnicity_concept_id':'ethnicity'},
                           drop_first=True, dtype=int)

# 储存为pickle格式
#topkl('person')

# %%
# =============================================================================
#  处理observation_period_df
# =============================================================================
loadcsv('observation_period')

# 查看数据类型和数据缺失情况
observation_period_df.info()

# 查看数据描述
observation_period_df.astype(str).describe()

# period_type_concept_id只有一种分类，放弃这个特征
observation_period_df.drop(columns=['period_type_concept_id'],
                           axis=1, inplace=True)

# 重设索引
observation_period_df.set_index(keys='person_id', inplace=True)

# 储存为pickle格式
#topkl('observation_period')

# %%
# =============================================================================
#  处理observation_df
# =============================================================================
loadcsv('observation')

# 查看数据类型和数据缺失情况
observation_df.info(null_counts=True)

# 查看数据描述
observation_df.astype(str).describe()

# 滤除目标周期外的数据
#observation_df = timefilter(observation_df, 'observation_date')

# value_as_string，unit_concept_id为空
# observation_type_concept_id，value_as_concept_id只有一种分类
# observation_date，visit_occurrence_id明显不相关，放弃这些特征
observation_df = observation_df[['person_id', 'observation_concept_id']]

# 数据透视，丢弃过稀疏特征
observation_df = pivot(observation_df,
                       'observation_concept_id', 'observation_')

# 储存为pickle格式
#topkl('observation')

# %%
# =============================================================================
#  处理drug_exposure_df
# =============================================================================
loadcsv('drug_exposure')

# 查看数据类型和数据缺失情况
drug_exposure_df.info(null_counts=True)

# 查看数据描述
drug_exposure_df.astype(str).describe()

# drug_type_concept_id只有一种分类，visit_occurrence_id明显不相关，放弃这些变量
drug_exposure_df.drop(columns=['drug_type_concept_id',
                               'visit_occurrence_id'], axis=1, inplace=True)

# 用对应药种的quantity中位数去填充缺失值
grouped = drug_exposure_df.groupby(by='drug_concept_id').median().fillna(1)
nan_ind = drug_exposure_df['quantity'].isnull()
nan_type = drug_exposure_df['drug_concept_id'][nan_ind]
drug_exposure_df['quantity'][nan_ind] = grouped.loc[nan_type, 'quantity']
del grouped, nan_ind, nan_type

# 滤除目标周期外的数据
drug_exposure_df = timefilter(drug_exposure_df, 'drug_exposure_start_date')

# 数据透视，丢弃过稀疏特征
drug_exposure_df = pivot(drug_exposure_df, 'drug_concept_id',
                         'drug_', values='quantity')

# 储存为pickle格式
#topkl('drug_exposure')

# %%
loadcsv('measurement')

# 查看数据类型和数据缺失情况
measurement_df.info(null_counts=True)

# 查看数据描述
measurement_df.astype(str).describe()

# measurement_type_concept_id只有一种分类，visit_occurrence_id明显不相关
# measurement_time,value_as_number等为空，放弃这些特征
measurement_df = measurement_df[['person_id', 'measurement_concept_id',
                                 'measurement_date']]

# 滤除目标周期外的数据
measurement_df = timefilter(measurement_df, 'measurement_date')

# 数据透视，丢弃过稀疏特征
measurement_df = pivot(measurement_df,
                       'measurement_concept_id', 'measurement_')

# 储存为pickle格式
#topkl('measurement')

# %%
loadcsv('condition_occurrence')

# 查看数据类型和数据缺失情况
condition_occurrence_df.info(null_counts=True)

# 查看数据描述
condition_occurrence_df.astype(str).describe()

# visit_occurrence_id明显不相关, 放弃这个特征
condition_occurrence_df.drop('visit_occurrence_id', axis=1, inplace=True)

# 滤除目标周期外的数据
condition_occurrence_df = timefilter(condition_occurrence_df,
                                     'condition_start_date')

# 数据透视，丢弃过稀疏特征
condition_occurrence_df = pd.concat([
    pivot(condition_occurrence_df,
          'condition_type_concept_id', 'condition_type_'),
    pivot(condition_occurrence_df,
          'condition_concept_id', 'condition_')
    ], axis=1)

# 储存为pickle格式
#topkl('condition_occurrence')

# %%
loadcsv('procedure_occurrence')

# 查看数据类型和数据缺失情况
procedure_occurrence_df.info(null_counts=True)

# 查看数据描述
procedure_occurrence_df.astype(str).describe()

# visit_occurrence_id明显不相关, 放弃这个特征
procedure_occurrence_df.drop('visit_occurrence_id', axis=1, inplace=True)

# 滤除目标周期外的数据
procedure_occurrence_df = timefilter(procedure_occurrence_df,
                                     'procedure_date')

# 对两列数据进行数据透视，并丢弃过稀疏特征
procedure_occurrence_df = pd.concat([
    pivot(procedure_occurrence_df,
          'procedure_type_concept_id', 'procedure_type_'),
    pivot(procedure_occurrence_df,
          'procedure_concept_id', 'procedure_')
    ], axis=1)

# 储存为pickle格式
#topkl('procedure_occurrence')
