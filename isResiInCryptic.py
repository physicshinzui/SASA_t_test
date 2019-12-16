#!/Users/siida/anaconda3/bin/python
import numpy as np
import sys
from MDAnalysis import Universe
import scipy.spatial.distance as distance

def isResInCryptic(pdbfile, ligname, cutoff=4.0):

    u      = Universe(pdbfile)
    ligand = u.select_atoms(f'resname {ligname}')
    protein= u.select_atoms(f'not resname {ligname}')

    resids = []
    s_cryptic, s_not_cryptic = [], []
    for iatom in ligand:
        for jatom in protein:

            distance = np.linalg.norm(iatom.position - jatom.position)
            if distance <= cutoff and jatom.resname in ['PHE','TRP','HIS','TYR']:
                #print(f'{iatom.name}-{iatom.resname}, {jatom.name}-{jatom.resname}{jatom.resid}, {distance}')
                resids.append(jatom.resid)
                s_cryptic.append(f'{jatom.resname}{jatom.resid}')

            elif distance > cutoff and jatom.resname in ['PHE','TRP','HIS','TYR']:
                s_not_cryptic.append(f'{jatom.resname}{jatom.resid}')

#    print(f'|s_cryptic| = {len(set(s_cryptic))}, |s_not_cryptic| = {len(set(s_not_cryptic))}')
#    print(f's_cryptic = {set(s_cryptic)}, s_not_cryptic = {set(s_not_cryptic)}')

    lis = 'sele resi ' + '+'.join([f'{i}' for i in set(resids)])
#    print(lis)

    return set(s_cryptic), set(s_not_cryptic)

def main():
    S_cryptic, S_not = isResInCryptic(sys.argv[1], 'X0B')

if __name__ == "__main__":
    main()
