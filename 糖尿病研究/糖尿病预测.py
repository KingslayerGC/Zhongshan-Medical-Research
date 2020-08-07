# %%
## 导入数据
import numpy as np
import pandas as pd
GDM_data = pd.read_excel(r"C:\Users\Mac\Desktop\过程\项目\医统\糖尿病研究\基因分析数据199例.xls",
              sheet_name='基因与GDM的关系资料', header=1)
Bloodsugar_data = pd.read_excel(r"C:\Users\Mac\Desktop\过程\项目\医统\糖尿病研究\基因分析数据199例.xls",
              sheet_name='基因与产后血糖转归的调查资料', header=1)

# %%
## 筛选数据
# 先筛掉一些显然无关或样本量严重不足的特征，得到所有备选特征和标签
X = GDM_data.drop(columns=['编号','姓名','备注（剖宫产原因等）',
                           '孕期合并症','Unnamed: 79'])
X.replace(to_replace={'/':np.nan}, inplace=True)
y = Bloodsugar_data[['产后血糖异常（有=1，无=0）']]
y.replace(to_replace={'/':np.nan}, inplace=True)

# %%
## 填充缺失值
X_necessary = X[['rs2297508（GG=0，GC=1，CC=2）','rs11868035(TT=0，TC=1，CC=2）','组别（妊娠期糖尿病GDM=1，正常=0）']]
cat_list = ['糖尿病家族史','一级亲属','二级亲属','父亲','母亲','父系','母系',
            '孕次（次）','产次（次）','新生儿性别（男=1，女=2）',
            '胎膜早破（无=0，有=1）','早产（无=0，有=1）','羊水过多（无=0，有=1）',
            '妊娠期高血压（无=0，有=1）','产后出血（无=0，有=1）',
            '胎膜早剥（无=0，有=1）','羊水过少（无=0，有=1）','流产（无=0，有=1）',
            '孕期并发症（无=0，有=1）','胎儿宫内生长受限/发育迟缓（无=0，有=1）',
            '巨大儿（无=0，有=1）','胎儿宫内窘迫（无=0，有=1）',
            '新生儿窒息（无=0，有=1）','新生儿黄疸/高胆红素血症（无=0，有=1）',
            '低体重儿或小于胎龄儿（无=0，有=1）','先天畸形（无=0，有=1）',
            '新生儿低血糖（无=0，有=1）','新生儿合并症（无=0，有=1）',
            '营养咨询或治疗（GDM患者进行营养治疗，NGT接受孕妇学校讲座。咨询或治疗=1，无=0）',
            'GDM孕妇用胰岛素治疗（是=1，否=0）']
X_cat = X[cat_list]     
X_num = X.drop(columns=cat_list)

# 定性变量用众位数填充
from sklearn.impute import SimpleImputer
impute = SimpleImputer(missing_values = np.nan, strategy='most_frequent')
y_tr = pd.DataFrame(impute.fit_transform(y), columns=y.columns)
X_necessary_tr = pd.DataFrame(impute.fit_transform(X_necessary),
                              columns=X_necessary.columns)
X_cat_tr = pd.DataFrame(impute.fit_transform(X_cat),
                        columns=X_cat.columns)

# 定量变量用随机森林填充
X_num_tr = X_num.copy()
df = pd.concat([X_necessary_tr, y_tr],axis=1)
from sklearn.ensemble import RandomForestRegressor
forest_reg = RandomForestRegressor(n_estimators=100)
forest_reg.get_params()
for col in list(X_num_tr.columns):
    if X_num_tr[col].isna().sum() == 0:
        continue
    fill = X_num_tr[col]
    Ytrain = fill[fill.notnull()]
    Ytest = fill[fill.isnull()]
    Xtrain = df.iloc[Ytrain.index,:]
    Xtest = df.iloc[Ytest.index,:]
    forest_reg.fit(Xtrain, Ytrain)
    fill[fill.isnull()] = forest_reg.predict(Xtest)
    X_num_tr[col] = fill

X_before_select = pd.concat([X_num_tr, X_cat_tr],axis=1)

# %%
from sklearn.ensemble import RandomForestClassifier
forest_clf = RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=42)
forest_clf.fit(X_before_select, y_tr)
# 输出所有特征的重要度
feature_importance = pd.DataFrame(forest_clf.feature_importances_.reshape(1,-1),
                                  columns=list(X_before_select.columns))
feature_importance.sort_values(by=0, axis=1, inplace=True, ascending=False)
# 选择排名靠前的特征
select_list = list(feature_importance.columns)[0:7] 
X_select = X_before_select[select_list]
X_tr = pd.concat([X_necessary, X_select],axis=1)

X = X_tr.values
y = y_tr.values.ravel()

# %%
from sklearn.linear_model import LogisticRegression
log_clf = LogisticRegression()

from sklearn.svm import SVC
svm_clf = SVC(kernel='rbf')

from sklearn.neighbors import KNeighborsClassifier
knn_clf = KNeighborsClassifier(algorithm='auto')

# %%
## 绘制ROC曲线
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_predict

log_y_scores = cross_val_predict(log_clf, X, y, cv=3, method="decision_function")

svm_y_scores = cross_val_predict(svm_clf, X, y, cv=3, method="decision_function")

knn_y_probas = cross_val_predict(knn_clf, X, y, cv=3, method="predict_proba")
knn_y_scores = knn_y_probas[:, 1]
    
# 绘制ROC曲线
def plot_roc_curve(fpr, tpr, linetype='b', label=None):
    plt.plot(fpr, tpr, linetype, linewidth=2, label=label)
    plt.plot([0, 1], [0, 1],'r')
    plt.axis([0, 1, 0, 1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

from sklearn.metrics import roc_curve

log_fpr, log_tpr, thresholds = roc_curve(y, log_y_scores)
plot_roc_curve(log_fpr, log_tpr, 'b:', label="logit")

svm_fpr, svm_tpr, thresholds = roc_curve(y, svm_y_scores)
plt.plot(svm_fpr, svm_tpr, label="SVM")

knn_fpr, knn_tpr, thresholds = roc_curve(y, knn_y_scores)
plt.plot(knn_fpr, knn_tpr, label="KNN")

plt.legend(loc="best")

plt.show()

# %%
# 计算AUC、f1-score和准确率
from sklearn.metrics import roc_auc_score
auc_log = roc_auc_score(y, log_y_scores)
auc_svm = roc_auc_score(y, svm_y_scores)
auc_knn = roc_auc_score(y, knn_y_scores)
# 计算准确率
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import cross_val_predict
for clf in[knn_clf, svm_clf, log_clf]:    
    y_pred = cross_val_predict(clf, X, y, cv=5)
    accuracy = accuracy_score(y, y_pred) 
    f1 = f1_score(y, y_pred)
    print(accuracy, f1)
auc_knn
auc_svm
auc_log

# %%
before_test = pd.concat([y_tr, X_select], axis=1)

from scipy.stats.mstats import ttest_ind
for col in list(before_test.columns):
    if col == '产后血糖异常（有=1，无=0）':
        continue
    n1 = before_test[col][before_test['产后血糖异常（有=1，无=0）'].isin([0])]
    n2 = before_test[col][before_test['产后血糖异常（有=1，无=0）'].isin([1])]
    stat, p = ttest_ind(n1, n2)
    print(col, p)
    
before_test['产后血糖异常（有=1，无=0）'].value_counts()


