#!/Users/siida/anaconda3/bin/python
import numpy as np
import sys
from MDAnalysis import Universe
import scipy.spatial.distance as distance

def classifyResiduesIntoTwo(apo_pdb, holo_pdb, ligname, cutoff=4.0):

    u_holo, u_apo = Universe(holo_pdb), Universe(apo_pdb)
    ligand = u_holo.select_atoms(f'resname {ligname}')
    holo   = u_holo.select_atoms(f'not resname {ligname}')
    apo    = u_apo.select_atoms(f'protein')

    resids = []
    S_cryptic, S_not_cryptic = [], []
    # -- calculate distances from atoms of a ligand to those of residues in an apo state
    # -- the aim is to detect residues in a cryptic site.
    # -- if the distance is less than a threshold (i.e., CRASHED!), then the aromatic residue is considered as cryptic one.
    for iatom in ligand:
        #for jatom in holo:
        for jatom in apo:
            distance = np.linalg.norm(iatom.position - jatom.position)

            if distance <= cutoff and jatom.resname in ['PHE','TRP','HIS','TYR']:
                #print(f'{iatom.name}-{iatom.resname}, {jatom.name}-{jatom.resname}{jatom.resid}, {distance}')
                resids.append(jatom.resid)
                S_cryptic.append(f'{jatom.resname}{jatom.resid}')

# -- a set of aromatic residue's names are generated here. note that this is specialised for aromatic residues
    S_not_cryptic = [ f'{residue.resname}{residue.resid}' for residue in holo.residues
                        if f'{residue.resname}{residue.resid}' not in S_cryptic and f'{residue.resname}' in ['PHE','TYR','TRP','HIS'] ]

# -- degug
    #print(set(S_not_cryptic).intersection(set(S_cryptic)), 'must be nothing.')

    print('This is for PyMol: sele resi ' + '+'.join([f'{i}' for i in set(resids)]))

    return set(S_cryptic), set(S_not_cryptic)

def main():
    S_cryptic, S_not = classifyResiduesIntoTwo(sys.argv[1], 'X0B')

if __name__ == "__main__":
    main()
