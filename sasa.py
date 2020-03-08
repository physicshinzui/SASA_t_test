#!/Users/siida/anaconda3/bin/python

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--traj', required=True, help='.pdb/.gro/.xtc/.trr ')
    p.add_argument('-s', '--top' , required=True, help='.pdb/.gro')
    p.add_argument('-rp', '--probe_radius'   ,type=float, default = 0.14, help='Default: 0.14')
    p.add_argument('-np', '--n_sphere_points' ,type=float, default = 100 , help='Default: 100')
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
    import mdtraj as md
    import pickle
    args = parser()
    trj_file, top = args.traj, args.top
    probe_radius, n_sphere_points = args.probe_radius, args.n_sphere_points
    traj = md.load(trj_file, top=top)
    print(traj)

    keys = get_residue_tag(top)
    print(keys)
    trj_sasa = md.shrake_rupley(traj,
                                probe_radius=probe_radius,
                                n_sphere_points=n_sphere_points, # the default value of gmx sasa is 24 
                                mode='residue')

    print('sasa data shape', trj_sasa.shape)
    print(len(trj_sasa[0,:]))

    i = 0
    sasa_dict = {}
    for key in keys:
        sasa_dict[key] = trj_sasa[:,i]
        i += 1
    fout1 = open('dict_sasa.pkl'   , 'wb')
    pickle.dump(sasa_dict, fout1)

if __name__ == "__main__":
    main()
