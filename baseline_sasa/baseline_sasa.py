#!/Users/siida/anaconda3/bin/python
import mdtraj as md

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--file', nargs='*', required=True, help='.pdb/.gro')
    args = p.parse_args()
    return args

def get_residue_tag(topology):
    from MDAnalysis import Universe
    u = Universe(topology)
    residue_tags = []
    for resn, resi in zip(u.residues.resnames, u.residues.resids):
        residue_tags.append(resn+str(resi))
    return residue_tags

def main():
    import numpy as np
    args = parser()
    file_names = args.file
    for file_name in file_names:
        traj = md.load(file_name, top=file_name)

        keys = get_residue_tag(file_name)
        print(keys)
        trj_sasa = md.shrake_rupley(traj,
                                    probe_radius=0.14,
                                    n_sphere_points=100, # the default value of gmx sasa is 24 
                                    mode='residue')

        print('sasa data shape', trj_sasa.shape)

        i = 0
        sasa_dict = {}
        for key in keys:
            sasa_dict[key] = trj_sasa[:,i]
            i += 1

        name = file_name.split('.')[0].upper()
        np.savetxt(f'sasa_{name}.dat', trj_sasa, fmt='%.2f',header=f'SASA in isolation of {name}' )

if __name__ == "__main__":
    main()
