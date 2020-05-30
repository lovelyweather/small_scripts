import xarray as xr
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from metpy import calc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV

'''calculate temperature C from perturbed potential temperature (K)'''
rp0 = 1.E-5   # p0=100000 的倒数
Rd = 287.04
Cp = 7.*Rd/2.
RCP = Rd/Cp
# Set random seed to ensure reproducible runs
RSEED = 50

'''read data and flatten as an 1D array'''
wrf=xr.open_dataset('wrfout_d01_2017-05-16_22_00_00')
refdata= xr.open_dataset('grid_ref_200.nc')
zdrfile= xr.open_dataset('grid_zdr_200.nc')
kdpfile= xr.open_dataset('grid_zdp_200.nc')
rhvfile= xr.open_dataset('grid_rhv_200.nc')
zdr=zdrfile['grid_zdr'][1:41,1:501,1:501].values.flatten()  #last no not included
kdp=kdpfile['grid_zdp'][1:41,1:501,1:501].values.flatten()
rhv=rhvfile['grid_rhv'][1:41,1:501,1:501].values.flatten()
ref=refdata['grid_ref'][1:41,1:501,1:501].values.flatten()
print(type(ref))

qv3d=wrf['QVAPOR'][0,:,:,:]
TEMP=wrf['T'][0,:,:,:]
Pres=wrf['P'][0,:,:,:]+wrf['PB'][0,:,:,:]
T= (TEMP + 300.) * ( Pres * rp0 )**RCP # get the temperature (K) of from WRF TEMP(perturbed potential temperature)
print('qv3d coords: {}, \n qv3d attrs: {}'.format(qv3d.coords,qv3d.attrs['units']))
qv3d.attrs['units']='dimensionless'
T.attrs={"units":'K'}
Pres.attrs={"units":'Pa'}
Pres.assign_attrs({"units": 'pressure'})
print(qv3d.attrs['units'],T.attrs,Pres.attrs)
rh=calc.relative_humidity_from_mixing_ratio(qv3d,T,Pres).flatten() # calculate relative humidity & transfer to 1D array
T=T.values.flatten()
P=Pres.values.flatten()

'''create the DataFrame for easier use'''
dataset = pd.DataFrame({'Reflectivity':ref,'ZDR':zdr,'KDP':kdp,'rhv':rhv,'Temp':T,'Pres':P,'RH':rh})
df=dataset[dataset['Reflectivity']>5]
#df['ZDRCol_Exist']=df['ZDR'].apply(lambda x:int(x>1)) & df['Temp'].apply(lambda x: int(x<273.55))

#print(df['ZDR,KDP,rhv'.split(',')].describe())
#print(df.isnull().any())
#print(df.isna().any())
#print(df.columns)
df=df.sample(frac=0.3,random_state=100,replace=True)  # sample 30% grids
X=df['Reflectivity,ZDR,KDP,rhv,Temp,Pres'.split(',')] ###,ZDRCol_Exist
Y=df['RH']
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.20, random_state = 42) #, stratify=qv)

rfr=RandomForestRegressor()
rfr.fit(X_train, y_train)
#rfr.feature_importances_
#rfr.get_params()

'''visualize the importance of features'''
feats = {}
for feature, importance in zip(X.columns, rfr.feature_importances_):
    feats[feature] = importance
    importances = pd.DataFrame.from_dict(feats, orient='index').rename(
        columns={0: 'MSE-Importance'})
    importances = importances.sort_values(by='MSE-Importance', ascending=False)
    importances = importances.reset_index()
    importances = importances.rename(columns={'index': 'Features'})
sns.barplot(x=importances['MSE-Importance'], y=importances['Features'], data=importances, color='skyblue')
plt.savefig('importance.jpg')

'''predict'''
y_pred = rfr.predict(X_test)
y_Ptrain = rfr.predict(X_train)
print('the train score is {}'.format(rfr.score(X_train, y_train)))
print('the predict score is {}'.format(rfr.score(X_test, y_test)))
#print(result.sample(frac=0.001,random_state=100,replace=True))
result=pd.DataFrame({'reflectivity':X_test['Reflectivity'],'ZDR':X_test['ZDR'],'KDP':X_test['KDP'],\
                     'rhv':X_test['rhv'],'Temp':X_test['Temp'],'Pres':X_test['Pres'],\
                     'real':y_test,'predict':y_pred})
train=pd.DataFrame({'reflectivity':X_train['Reflectivity'],'ZDR':X_train['ZDR'],'KDP':X_train['KDP'],\
                     'rhv':X_train['rhv'],'Temp':X_train['Temp'],'Pres':X_train['Pres'],\
                     'real':y_train,'predict':y_Ptrain})

#result.corr()
#result.describe()
#sns.distplot(result['Temp'],kde=True)
#result[abs(result['real']-result['predict'])<0.2].describe()
plt.clf()
plt.figure(figsize=(8,8))
#cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
sns.scatterplot('real','predict',size='reflectivity',data=result,sizes=(5, 100),palette="coolwarm",color='black',linewidth=0,alpha=0.3) #,edgecolors='red')
plt.xlabel('Real',fontsize=15)
plt.ylabel('Predict',fontsize=15)
plt.savefig('test.jpg')

'''below is for tuning parameters, to be continued'''
'''
n_estimators = [int(x) for x in np.linspace(start = 100, stop = 500, num = 5)]
max_features = ['log2', 'sqrt']
max_depth = [int(x) for x in np.linspace(start = 3, stop = 15, num = 3)]
min_samples_split = [int(x) for x in np.linspace(start = 2, stop = 50, num = 5)]
min_samples_leaf = [int(x) for x in np.linspace(start = 2, stop = 50, num = 5)]
bootstrap = [True, False]
param_dist = {'n_estimators': n_estimators, 'max_features': max_features, 'max_depth': max_depth, \
              'min_samples_split': min_samples_split, 'min_samples_leaf': min_samples_leaf, \
              'bootstrap': bootstrap}
rfr=RandomForestRegressor()
rs = RandomizedSearchCV(rfr,  param_dist, n_iter = 100,  cv = 3,  verbose = 1,  n_jobs=-1,  random_state=0)
rs.fit(X_train, y_train)
print(rs.best_params_)

#output is as follows:
# {'n_estimators': 200, 'min_samples_split': 14, 'min_samples_leaf': 2,
# \'max_features': 'log2', 'max_depth': 15, 'bootstrap': False}
'''