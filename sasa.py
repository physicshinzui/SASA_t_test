#!/Users/siida/anaconda3/bin/python
import numpy as np
import mdtraj as md
import sys
import pickle

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--traj', required=True, help='.pdb/.gro/.xtc/.trr ')
    p.add_argument('-s', '--top' , required=True, help='.pdb/.gro')
    args = p.parse_args()
    return args

def get_residue_tag(topology):
    from MDAnalysis import Universe
    u = Universe(topology)
    residue_tags = []
    for resn, resi in zip(u.residues.resnames, u.residues.resids):
        residue_tags.append(resn+str(resi))
    return residue_tags

def standardized_sasa(dict_sasa):
    # the sasa [nm^2] were calculated for the atoms without hydrogens via
    # "compute -> surfae area -> solvent accessible" in PyMol
    sasa_values_of_naked_aa = {'PHE': 3.280,
                               'TYR': 3.433,
                               'TRP': 3.737,
                               'HIS': 3.136,
                               'ARG': 3.675}
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
    args = parser()
    trj_file, top = args.traj, args.top
    traj = md.load(trj_file, top=top)
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

    fout1 = open('dict_sasa.pkl'   , 'wb')
    fout2 = open('st_dict_sasa.pkl', 'wb')
    pickle.dump(sasa_dict   , fout1)
    pickle.dump(st_sasa_dict, fout2)

if __name__ == "__main__":
    main()
