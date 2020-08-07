# %%
import pandas as pd
GDM_data = pd.read_excel(r"C:\Users\Mac\Desktop\过程\项目\医统\糖尿病研究\基因分析数据199例.xls",
              sheet_name='基因与GDM的关系资料', header=1)
Bloodsugar_data = pd.read_excel(r"C:\Users\Mac\Desktop\过程\项目\医统\糖尿病研究\基因分析数据199例.xls",
              sheet_name='基因与产后血糖转归的调查资料', header=1)

# %%
df = pd.concat([GDM_data[['组别（妊娠期糖尿病GDM=1，正常=0）',
                          'rs11868035(TT=0，TC=1，CC=2）',
                          'rs2297508（GG=0，GC=1，CC=2）']],
          Bloodsugar_data[['产后血糖异常（有=1，无=0）']]], axis=1)
df.rename(columns={'组别（妊娠期糖尿病GDM=1，正常=0）':'GDM',
                   'rs11868035(TT=0，TC=1，CC=2）':'rs1',
                   'rs2297508（GG=0，GC=1，CC=2）':'rs2',
                   '产后血糖异常（有=1，无=0）':'bloodsugar'}, inplace=True)
df = df[~df['bloodsugar'].isin(['/'])].astype(float)

X = df.drop(columns=['bloodsugar']).values
y = df[['bloodsugar']].values.ravel()

# %%
from sklearn.linear_model import LogisticRegression
log_clf = LogisticRegression()

from sklearn.svm import SVC
svm_clf = SVC(kernel='rbf')

from sklearn.neighbors import KNeighborsClassifier
knn_clf = KNeighborsClassifier(algorithm='auto')

from sklearn.model_selection import cross_val_predict

log_y_scores = cross_val_predict(log_clf, X, y, cv=3, method="decision_function")

svm_y_scores = cross_val_predict(svm_clf, X, y, cv=3, method="decision_function")

knn_y_probas = cross_val_predict(knn_clf, X, y, cv=3, method="predict_proba")
knn_y_scores = knn_y_probas[:, 1]



# %%
# 计算AUC
import matplotlib.pyplot as plt

from sklearn.metrics import roc_auc_score
auc_log = roc_auc_score(y, log_y_scores)
auc_svm = roc_auc_score(y, svm_y_scores)
auc_knn = roc_auc_score(y, knn_y_scores)

from sklearn.metrics import accuracy_score, f1_score

score_pd = pd.DataFrame(index=range(2), columns=['逻辑回归', 'SVM', 'KNN'])


svm_clf.fit(X, y)
y_pred = log_clf.predict(X)
accuracy = accuracy_score(y, y_pred) 
f1 = f1_score(y, y_pred)
print(accuracy, f1)

    
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
df['rs1'][df['rs1'].isin([2])] = 1
df['rs2'][df['rs2'].isin([2])] = 1

# %%
df1 = df[df['bloodsugar'].isin([0])]
df2 = df[df['bloodsugar'].isin([1])]
from scipy.stats.mstats import ttest_ind
for str in ['GDM','rs1','rs2']:
    stat, p = ttest_ind(df1[str], df2[str])
    print(str,'   ',p)


df['bloodsugar'].value_counts()

# %%
df['combine'] = 10
'''
change='GDM'
control='rs'
flag=1

if flag==0:
    df = df[df[control].isin([0])]
else:
    df = df[~df[control].isin([0])]

'''
df['combine'][df['GDM'].isin([0]) & df['rs'].isin([0])] = 0
df['combine'][~df['GDM'].isin([0]) & ~df['rs'].isin([0])] = 1
df = df[~df['combine'].isin([10])]

X=df[['combine']]

y=df[['bloodsugar']]

from sklearn.linear_model import LogisticRegression
log_clf = LogisticRegression()
log_clf.fit(X,y)
coef=list(log_clf.coef_)
import math
OR = math.exp(log_clf.coef_[0][0])
#import scipy.stats as ss
#e = math.exp((1/61+0.5+1/12+0.1)**0.5*ss.norm.ppf(0.975,0,1))
contingency_table = pd.crosstab(df['combine'], df['bloodsugar'])
from scipy.stats import chi2_contingency
chi2, p, dof, ex = chi2_contingency(contingency_table)
print(OR,'\n',p)

# %%
from statsmodels.stats.contingency_tables import Table
table = Table(contingency_table)
table.test_nominal_association()
print(table.test_nominal_association())
before_contingency

df.to_excel(r'C:\Users\Mac\Desktop\数据.xlsx')
'''