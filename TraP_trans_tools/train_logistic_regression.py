import generic_tools
import numpy as np
import random
from scipy import optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties

def shuffle_datasets(data):
    shuffled=[]
    val_list=range(len(data))
    random.shuffle(val_list)
    for row in range(len(data)):
        shuffled.append(data[val_list[row]])
    shuffled=np.array(shuffled)
    return shuffled

def create_datasets(data, n, m):
    shuffle_datasets(data)
    train=data[:n,:]
    valid=data[n:m,:]
    test=data[m:,:]
    return train, valid, test

def create_X_y_arrays(data):
    i=data.shape[1]-1
    X=np.matrix(data[:,:i])
    X = np.c_[np.ones(len(X)), X]
    y=np.matrix(data[:,i])
    return X, y

def sigmoid(z):
    g = 1/(1+np.exp(-z))
    return g

def reg_cost_func(theta, X, y, lda):
    m=y.shape[1]
    J=0
    sig=sigmoid(X * np.c_[theta])
    J=(1/float(m))*(-y.dot(np.log(sig)) - (1-y).dot(np.log(1-sig)))
    temp=np.copy(theta)
    temp[0]=0.0
    J=J+(lda/(2*m))*np.sum(np.multiply(temp,temp))
    Jnum=np.array(J)[0][0]
    return Jnum

def quadratic_features(X):
    X_quad=np.matrix(np.zeros((X.shape[0],X.shape[1]*2)))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            X_quad[i,j]=X[i,j]
            X_quad[i,j+X.shape[1]]=X[i,j]**2
    return X_quad

def learning_curve(X, y, Xval, yval, lda, options):
    m=X.shape[0]
    n=Xval.shape[0]
    error_train=np.zeros((m))
    error_val=np.zeros((m))
    for i in range(1,m):
        theta=np.zeros((X.shape[1]))
        initial_theta=np.zeros((X.shape[1]))
        theta, cost, _, _, _ = optimize.fmin(lambda t: reg_cost_func(t,X[:i,:],y[:,:i],lda), initial_theta, **options)
        error_train[i]=check_error(X,y,theta)
        error_val[i]=check_error(Xval,yval,theta)
    return error_train, error_val, theta

def check_error(X,y,theta):
    error = (1./(2.*float(X.shape[0])))*np.sum(np.power((sigmoid(X * np.c_[theta])-y.T),2))
    return error

def validation_curve(X, y, Xval, yval,options):
    lambda_vec = np.array([10.**a for a in np.arange(-5,5,0.1)])
    error_train=np.zeros((lambda_vec.shape[0]))
    error_val=np.zeros((lambda_vec.shape[0]))
    theta_rec=np.zeros((lambda_vec.shape[0],X.shape[1]))
    m=X.shape[0]
    n=Xval.shape[0]
    for i in range(0,lambda_vec.shape[0]):
        theta=np.zeros((X.shape[1]))
        initial_theta=np.zeros((X.shape[1]))
        lda=lambda_vec[i]
        theta, cost, _, _, _ = optimize.fmin(lambda t: reg_cost_func(t,X,y,lda), initial_theta, **options)
        theta_rec[i]=theta
        error_train[i]=check_error(X,y,theta)
        error_val[i]=check_error(Xval,yval,theta)
        print theta
    min_err_val = min(error_val)
    for i in range(0,lambda_vec.shape[0]):
        if error_val[i] == min_err_val:
            lda=lambda_vec[i]
    return error_train, error_val, lambda_vec, lda

def plotLC(num, error_train, error_val, fname, xlog, ylog, xlabel):
    plt.figure()
    plt.plot(num, error_train, 'b-')
    plt.plot(num, error_val, 'g-')
    if ylog:
        plt.yscale('log')
    if xlog:
        plt.xscale('log')
    plt.xlabel(xlabel)
    plt.ylabel('Error')
    plt.axis([min(num)*0.8, max(num)*1.2, 1e-4,2e-2])
    plt.legend(['training', 'validation'], loc=4)
    plt.savefig('LR_'+fname+'_curve.png')
    plt.close()
    return

def classify_data(X,y,theta):
    tp=0
    fp=0
    fn=0
    tn=0
    classified_data=[]
    predictions=predict(X,theta)
    y=y.T
    for i in range(predictions.shape[0]):
        if predictions[i] > 0.5 and y[i] == 1:
            tp=tp+1
            classified_data.append([X[i,1],X[i,2],X[i,3],X[i,4],1])
        elif predictions[i] > 0.5 and y[i] == 0:
            fp=fp+1
            classified_data.append([X[i,1],X[i,2],X[i,3],X[i,4],2])
        elif predictions[i] < 0.5 and y[i] == 1:
            fn=fn+1
            classified_data.append([X[i,1],X[i,2],X[i,3],X[i,4],3])
        elif predictions[i] < 0.5 and y[i] == 0:
            tn=tn+1
            classified_data.append([X[i,1],X[i,2],X[i,3],X[i,4],4])
    return tp, fp, fn, tn, classified_data

def predict(X,theta):
    predictions=sigmoid(X * np.c_[theta])
    return predictions
