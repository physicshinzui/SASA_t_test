#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
import pickle
from scipy import stats

def two_sample_t_test(A, B, target='PHE66', alpha=0.05):
    # inpus a and b are dict 

    counter = 0
    n_reject= 0
    rejected_residues     = []
    non_rejected_residues = []
    fout = open(f'{target}_to_others.out','w')
    fout.write('# pair, t-statisitcs, p-value\n')
    for key in A.keys():
        header = key[0:3]

        if key == target:
            continue

        elif header == "PHE" or header == "TYR" or header == "TRP" or header == "HIS":

            t2, p2 = stats.ttest_ind(A[target], B[key])

            fout.write(f'{target}-{key}, {t2:>12.8}, {p2:<12.8}\n')
            if p2 < alpha:
                print(f'{target}-{key} can be different. (i.e., null hypo. was rejected.)')
                rejected_residues.append(key)
                n_reject += 1

            else:
                print(f'{target}-{key} might be identical.')
                non_rejected_residues.append(key)
                
            print(f'    {t2:<10.5}, {p2:<10.5}')
            counter += 1
    fout.close()
    print(f'{n_reject} {counter} {(n_reject / counter) * 100.0}%')
    print(rejected_residues)
    print(non_rejected_residues)

def main():
    fin = open('st_dict_sasa.pkl', 'rb')
    #fin = open('dict_sasa.pkl', 'rb')
    sasa = pickle.load(fin)

    two_sample_t_test(sasa, sasa) 
#    t2, p2 = stats.ttest_ind(sasa['PHE66'], sasa['HIS74'])
#    print(t2, p2)

if __name__ == "__main__":
    main()
