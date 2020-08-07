# %%
# =============================================================================
#  基础模板和参数设置
# =============================================================================
import pandas as pd
import os
import re

pkl_path = r"C:\Users\Mac\Desktop\过程\项目\医统\EHR challenge\临时数据-pickle"

# =============================================================================
#  文件读取和前期处理
# =============================================================================
# 读取pkl文件，若不存在则合成
for parents, _, filename in os.walk(pkl_path):
    path = parents + "\\total.pkl"
    if "total.pkl" in filename:
        total_df = pd.read_pickle(path)
    else:
        total_df = None
        for file in filename:
            v = re.match(r'(\w*)\.\w*', file).group(1)
            total_df = pd.concat([total_df, pd.read_pickle(
                os.path.join(parents, file))], join='outer', axis=1)

        # 丢弃没有visit记录的样本,并进行零填充
        total_df = total_df.loc[total_df['visit_end_date'].notnull()]

        # 在最后一次visit后六个月内死亡的样本视为真阳性
        from datetime import timedelta
        gap_df = total_df['death_date'] - total_df['visit_end_date']
        total_df['death'] = 0
        total_df['death'].loc[gap_df < timedelta(6 * 30)] = 1


        # 计算最后一次visit时的年龄
        total_df['age'] = ((total_df['visit_end_date'] - total_df['birth_date']
                    ).apply(lambda x : x.days) / 365).astype(int)

        # 丢弃明显不相关变量
        total_df = total_df.drop(['death_date', 'birth_date', 'visit_end_date',
                                  'observation_period_start_date',
                                  'observation_period_end_date'], axis=1).fillna(0)
        
        # 将所有特征合并储存
        total_df.to_pickle(path)
        
# 采取欠抽样平衡死亡人数偏斜情况
death_df = total_df[total_df['death'] == 1]
sample_df = pd.concat(
    [total_df[total_df['death'] == 1],
     total_df[total_df['death'] == 0].sample(
         death_df.shape[0], random_state=42)], axis=0)

# 按照8：2的比例切割训练集和验证集
from sklearn.model_selection import StratifiedShuffleSplit
split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
for train_index, test_index in split.split(sample_df, sample_df['death']):
    X_train = sample_df.drop('death', axis=1).iloc[train_index]
    y_train = sample_df['death'].iloc[train_index]
    X_test = sample_df.drop('death', axis=1).iloc[test_index]
    y_test = sample_df['death'].iloc[test_index]

# %%
# =============================================================================
#  特征工程
# =============================================================================
# PCA降维
from sklearn.decomposition import PCA
pca = PCA(n_components=0.9)
pca.fit(X_train)

X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)

 # %%
# =============================================================================
#  模型与效果
# =============================================================================
from sklearn.metrics import accuracy_score
from prettytable import PrettyTable

## 一个展示准确度的表格
table = PrettyTable([
    'Decomposition','Method','Train Accuracy','Test Accuracy'])

# 输出准确度函数
def accuracy_summary(clf_name, clf, X, X_test, decomposition):
    clf.fit(X, y_train)
    y_pred_train = clf.predict(X)
    y_pred_test = clf.predict(X_test)
    table.add_row([decomposition, clf_name,
                   str(accuracy_score(y_train, y_pred_train)*100)[:5] + '%',
                   str(accuracy_score(y_test, y_pred_test)*100)[:5] + '%'])

## 训练一些分类器查看效果
# 线性svm
from sklearn.svm import LinearSVC
svm_clf = LinearSVC(C=0.1)

# 逻辑回归
from sklearn.linear_model import LogisticRegression
log_clf = LogisticRegression()

# 随机森林
from sklearn.ensemble import RandomForestClassifier
forest_clf = RandomForestClassifier(
    n_estimators=100, max_depth=30, random_state=42)

# 极端随机树
from sklearn.ensemble import ExtraTreesClassifier
extree_clf = ExtraTreesClassifier(
    n_estimators=100, max_depth=30, random_state=42)

# 训练并记录结果
accuracy_summary('Linear SVM', svm_clf, X_train_pca, X_test_pca, 'PCA')
accuracy_summary('Logistic Regression', log_clf,
                 X_train_pca, X_test_pca, 'PCA')
accuracy_summary('Random Forest', forest_clf, X_train, X_test, 'NONE')
accuracy_summary('Extra Trees', extree_clf, X_train, X_test, 'NONE')

# 查看模型效果
print(table)


