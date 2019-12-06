#!/Users/siida/anaconda3/bin/python
import numpy as np
import mdtraj as md
from MDAnalysis import Universe
import sys
import pickle

def get_residue_tag(topology):
    u = Universe(topology)
    residue_tags = []
    for resn, resi in zip(u.residues.resnames, u.residues.resids):
        residue_tags.append(resn+str(resi))
    return residue_tags

def standardized_sasa(dict_sasa):
    # the sasa [nm^2] were calculated for the atoms without hydrogens
    sasa_values_of_naked_aa = {'PHE': 3.12518,
                               'TYR': 3.37539,
                               'TRP': 3.70723,
                               'HIS': 3.17074}
    st_dict_sasa = {}
    print('Note: I do not standardize sasa of non-aromatic residues cuz they are out of my scope.')
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

        else:
            st_dict_sasa[key] = dict_sasa[key]
    
    return st_dict_sasa

def main():
    trj_file = 'input_traj/cat_npt_prod_skip100.xtc'
    #trj_file = 'test/test.xtc'
    #trj_file = '../../../data/bcl_xl@ambient/data/traj/cat_npt_prod.xtc'
    top = 'template.pdb'
    traj = md.load(trj_file, top='template.pdb')
    print(traj)

    keys = get_residue_tag(top)
    print(keys)
    trj_sasa = md.shrake_rupley(traj, 
                                probe_radius=0.14, 
                                n_sphere_points=100, 
                                mode='residue', 
                                change_radii=None, 
                                get_mapping=False)
   
    print('sasa data shape', trj_sasa.shape)
    print(len(trj_sasa[0,:]))

    i = 0
    sasa_dict = {}
    for key in keys:
        sasa_dict[key] = trj_sasa[:,i]
        i += 1

    st_sasa_dict = standardized_sasa(sasa_dict)
#    print(sasa_dict['TRP98'][0:3], st_sasa_dict['TRP98'][0:3])

    fout1 = open('dict_sasa.pkl', 'wb')
    fout2 = open('st_dict_sasa.pkl', 'wb')
    pickle.dump(sasa_dict, fout1)
    pickle.dump(st_sasa_dict, fout2)

if __name__ == "__main__":
    main()

