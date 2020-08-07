# %%
import pandas as pd
#import numpy as np
#total_matrix=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\total_matrix_2.csv',\
#                         header=0,index_col=0,sep=',',low_memory=False)
total_matrix=pd.read_csv(r'C:\Users\Mac\Desktop\过程\项目\total_matrix_wil0.1.csv',\
                         header=0,index_col=0,sep=',',low_memory=False)

# %%
X=total_matrix.drop(['death'],axis=1)
X=(X-X.mean())/(X.std())
Y=total_matrix['death']

# %%
import random
alive_list=list(total_matrix[total_matrix['death']==0].index)
death_list=list(total_matrix[total_matrix['death']==1].index)
sample_list=random.sample(alive_list,3000)+death_list
sample_matrix=total_matrix[total_matrix.index.isin(sample_list)]
X=sample_matrix.drop(['death'],axis=1)
X=(X-X.mean())/(X.std())
Y=sample_matrix['death']

#%%
from sklearn.decomposition import PCA
pca = PCA()
pca.fit(X)
contribute_rate_list=pca.explained_variance_ratio_.tolist()
k=0
for i in range(len(contribute_rate_list)):
    k=k+contribute_rate_list[i]
    if k>0.9:
        print(i)
        break
'''
pca = PCA(1000)
pca.fit(X)  #降低维度
X=pd.DataFrame(pca.transform(X))
'''

# %%
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
import sklearn.metrics as sm
x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=3)
model = LogisticRegression(C=10)
model.fit(x_train, y_train)
train_score = model.score(x_train, y_train)
y_test_pre = model.predict(x_test)
y_train_pre = model.predict(x_train)
# %%
log_train_score = model.score(x_train, y_train)
print("LOG Train_Accuracy = %f"%log_train_score)
log_accuracy = sm.accuracy_score(y_test,y_test_pre) 
print("LOG CV_Accuracy_Score = %f"%log_accuracy)  
log_precision = sm.precision_score(y_test,y_test_pre, average=None)
print("LOG Precision = ",log_precision)
log_recall = sm.recall_score(y_test,y_test_pre, average=None)
print("LOG Recall = ",log_recall) 
log_f1_score = sm.f1_score(y_test,y_test_pre)
print("LOG F1-Score  = %f"%log_f1_score)
# %%
from sklearn.svm import SVC
import pandas as pd
from sklearn.model_selection import train_test_split
import sklearn.metrics as sm
x_train2,x_test2,y_train2,y_test2 = train_test_split(X,Y,test_size=0.3, random_state=3)
clf = SVC (C=10,kernel='rbf',gamma='auto')
clf.fit(x_train2,y_train2)
y_test_pre2 = clf.predict(x_test2)
# %%
svm_train_score = clf.score(x_train2, y_train2)
print("SVM Train_Accuracy = %f"%svm_train_score)
svm_accuracy = sm.accuracy_score(y_test2,y_test_pre2) 
print("SVM CV_Accuracy_Score = %f"%svm_accuracy)  
svm_precision = sm.precision_score(y_test2,y_test_pre2, average=None)
print("SVM Precision = ",svm_precision)
svm_recall = sm.recall_score(y_test2,y_test_pre2, average=None)
print("SVM Recall = ",svm_recall) 
svm_f1_score = sm.f1_score(y_test2,y_test_pre2)
print("SVM F1-Score  = %f"%svm_f1_score)

# %%






# %%
import numpy as np
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import pandas as pd
total_matrix=pd.read_csv(r'C:\Users\Mac\Desktop\过程\编程\项目\EHR challenge\生成表格\total_matrix_wil0.1.csv',\
                         header=0,index_col=0,sep=',',low_memory=False)
# 定义选取一批次数据的函数
def shuffle_batch(X, y, batch_size):
    rnd_idx = np.random.permutation(len(X))
    n_batches = len(X) // batch_size
    for batch_idx in np.array_split(rnd_idx, n_batches):
        X_batch, y_batch = X[batch_idx], y[batch_idx]
        yield X_batch, y_batch
# 定义初始化函数
def reset_graph(seed=42):
    tf.reset_default_graph()
    tf.set_random_seed(seed)
    np.random.seed(seed)

# %%
import random
alive_list=list(total_matrix[total_matrix['death']==0].index)
death_list=list(total_matrix[total_matrix['death']==1].index)
sample_list=random.sample(alive_list,3000)+death_list

sample_matrix=total_matrix[total_matrix.index.isin(sample_list)]
X_train=sample_matrix.drop(['death'],axis=1).values
y_train=np.array(sample_matrix['death'])

valid_matrix=total_matrix.loc[total_matrix.index.drop(sample_list)]
X_valid=valid_matrix.drop(['death'],axis=1).values
y_valid=np.array(valid_matrix['death'])

# %%
reset_graph()
from functools import partial
n_inputs = 319
n_hidden1 = 100
n_hidden2 = 30
n_outputs = 2
learning_rate = 0.01
n_epochs = 40
batch_size = 100
dropout_rate=0.5

X = tf.placeholder(tf.float32, shape=(None, n_inputs), name='X')
y = tf.placeholder(tf.int32, shape=(None), name='y')
training = tf.placeholder_with_default(False, shape=(), name='training')

X_drop = tf.layers.dropout(X, dropout_rate, training=training)
training = tf.placeholder_with_default(False, shape=(), name='training')
he_init=tf.variance_scaling_initializer()
my_dense_layer = partial(tf.layers.dense, kernel_initializer=he_init)
my_batch_norm_layer = partial(tf.layers.batch_normalization,
                              training=training, momentum=0.9)

with tf.name_scope('dnn'):
    hidden1 = my_dense_layer(X_drop, n_hidden1, name='hidden1')
    bn1 = tf.nn.elu(my_batch_norm_layer(hidden1))
    hidden1_drop = tf.layers.dropout(bn1, dropout_rate, training=training)
    
    hidden2 = my_dense_layer(hidden1_drop, n_hidden2, name='hidden2')
    bn2 = tf.nn.elu(my_batch_norm_layer(hidden2))
    hidden2_drop = tf.layers.dropout(bn2, dropout_rate, training=training)
    logits_before_bn = my_dense_layer(hidden2_drop, n_outputs, name='outputs')
    logits = my_batch_norm_layer(logits_before_bn)

with tf.name_scope('loss'):
    xentropy = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=y,
                                                              logits=logits)
    loss = tf.reduce_mean(xentropy, name='loss')
with tf.name_scope('train'):
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate)
    training_op = optimizer.minimize(loss)
with tf.name_scope('eval'):
    correct = tf.nn.in_top_k(logits, y, 1)
    accuracy = tf.reduce_mean(tf.cast(correct, tf.float32), name="accuracy")
init = tf.global_variables_initializer()
saver = tf.train.Saver()

extra_update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
with tf.Session() as sess:
    init.run()
    for epoch in range(n_epochs):
        for X_batch, y_batch in shuffle_batch(X_train, y_train, batch_size):
            sess.run([training_op, extra_update_ops],
                     feed_dict={training:True, X:X_batch, y:y_batch})
        accuracy_train = accuracy.eval(feed_dict={X:X_train, y:y_train})
        accuracy_valid = accuracy.eval(feed_dict={X:X_valid, y:y_valid})
        print('Training score', accuracy_train, 'Valid score', accuracy_valid)











