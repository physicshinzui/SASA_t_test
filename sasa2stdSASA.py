#!/Users/siida/anaconda3/bin/python
import mdtraj as md

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--sasa' , required=True, help='sasa.pkl')
    p.add_argument('-b', '--base', required=True, help='baseline.inp')
    args = p.parse_args()
    return args

def read_baseline(filename):
    with open(filename, 'r') as fin:
        dict_sasa = {}
        for line in fin:
            line = line.split(':')
            resn, sasa = line[0], line[1]
            dict_sasa[resn] = float(sasa)
    return dict_sasa

def standardized_sasa(dict_sasa, baseline):
    # The values were computed by mdtraj's sasa class.
    # unit, nm
    sasa_values_of_naked_aa = baseline

    st_dict_sasa = {}
    print('Note: Standardized SASA of non-aromatic residues is not computed.')
    for key in dict_sasa.keys():

        if key[0:3] == 'PHE':
            print(key)
            st_dict_sasa[key] = dict_sasa[key] / sasa_values_of_naked_aa['PHE']

        elif key[0:3] == 'TYR':
            print(key)
            st_dict_sasa[key] = dict_sasa[key] / sasa_values_of_naked_aa['TYR']

        elif key[0:3] == 'TRP':
            print(key)
            st_dict_sasa[key] = dict_sasa[key] / sasa_values_of_naked_aa['TRP']

        elif key[0:3] == 'HIS':
            print(key)
            st_dict_sasa[key] = dict_sasa[key] / sasa_values_of_naked_aa['HIS']

        elif key[0:3] == 'ARG':
            print(key)
            st_dict_sasa[key] = dict_sasa[key] / sasa_values_of_naked_aa['ARG']

        else:
            st_dict_sasa[key] = dict_sasa[key]

    return st_dict_sasa

def main():
    import pickle
    import sys
    args = parser()
    inp, baseline = args.sasa, args.base
    print(inp, baseline)

    fsasa     = open(inp,'rb')
    sasa_dict = pickle.load(fsasa)
    base      = read_baseline(baseline)

    print('--- baseline ---')
    print(base)

    st_sasa_dict = standardized_sasa(sasa_dict, base)
    fout = open('st_dict_sasa.pkl', 'wb')
    pickle.dump(st_sasa_dict, fout)

if __name__ == "__main__":
    main()
