#!/Users/siida/anaconda3/bin/python
import pickle

f1, f2 = open('dict_sasa.pkl','rb'), open('st_dict_sasa.pkl','rb')
sasa = pickle.load(f1)
st_sasa = pickle.load(f2)

print(sasa['PHE66'][0], st_sasa['PHE66'][0])
