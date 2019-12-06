#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt

f1, f2 = open('dict_sasa.pkl','rb'), open('st_dict_sasa.pkl','rb')
sasa = pickle.load(f1)
st_sasa = pickle.load(f2)

#plt.plot(sasa['PHE66'])
#plt.plot(st_sasa['PHE66'])
plt.hist(st_sasa['PHE66'], bins=[i for i in np.arange(0, 1.0, 0.02)], histtype='step', linewidth = 2)
plt.hist(st_sasa['PHE58'], bins=[i for i in np.arange(0, 1.0, 0.02)], histtype='step', linewidth = 2)
plt.hist(st_sasa['TRP98'], bins=[i for i in np.arange(0, 1.0, 0.02)], histtype='step', linewidth = 2)
plt.hist(st_sasa['HIS138'], bins=[i for i in np.arange(0, 1.0, 0.02)], histtype='step', linewidth = 2)
plt.show()
