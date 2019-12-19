#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from MDAnalysis import Universe

def ratioOfPeToPb(std_SASA, buried_upper_limit=0.1):
    """
    Args:
        std_SASA: [dict] standardized SASA like {PHE105:ASA, ..., TYR200:ASA}
        buried_upper_limit: [float] upper limit of buriedness (standardized SASA) [no unit]
    Return:
        Reb: Ratio of exposed probability to buried one. Small value indicates buried states are more dominant than exposed.
    """
    bin_start, bin_end, interval = 0, 1, 0.01
    hist = np.histogram(std_SASA, bins=[i for i in np.arange(bin_start,bin_end,interval)])
    probs, bin_edges = hist[0], hist[1]
    half_width       = 0.5 * abs(bin_edges[0] - bin_edges[1])
    hist_xpoints     = np.array([edge+half_width for edge in bin_edges[:-1]])

    index_Pe = np.where(hist_xpoints > buried_upper_limit)[0]
    index_Pb = np.where(hist_xpoints <= buried_upper_limit)[0]
    Pe       = [probs[i] for i in index_Pe]
    xe       = hist_xpoints[hist_xpoints > buried_upper_limit]
    Pb       = [probs[j] for j in index_Pb]
    xb       = hist_xpoints[hist_xpoints < buried_upper_limit]

    Reb = np.sum(Pe) / np.sum(Pb) # R_eb, the more it is, the more exposed and vice versa.

    plt.plot(xe, Pe)
    plt.plot(xb, Pb)
    plt.axvline(x=buried_upper_limit, c='black')
    plt.xlim(0,1.0)
    return Reb

def isASAVaried(Reb, Reb_lower=0.1, Reb_upper=50.0):
    """
    Args:
        Reb_lower:
        Reb_upper:
    Returns:
        True or False
    """
    # -- Reb thresholds
    if Reb > Reb_lower and Reb < Reb_upper:
        return True

    elif Reb >= Reb_upper or Reb <= Reb_lower:
        return False

def predictCrypticResidues(standardized_SASA):
    S_aromatic_resname = set(['PHE','TRP','TYR','HIS'])
    S_varied     = []
    S_not_varied = []
    for key in standardized_SASA:
        # -- I just concentrate on aromatic residues.
        # -- the iput pkl file must be modified if you want to check the other residues out, because they've not been normalised.
        if key[0:3] in S_aromatic_resname:
            # -- Calculate the ratio (Reb) of exposed states to buried ones.
            Reb = ratioOfPeToPb(standardized_SASA[key])

            if isASAVaried(Reb):
                S_varied.append(key)

            else:
                S_not_varied.append(key)

    return set(S_varied), set(S_not_varied)

def classifyResiduesIntoTwo(apo_pdb, holo_pdb, ligname, cutoff=4.0):
    S_aromatic_resname = set(['PHE','TRP','TYR','HIS'])
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
        for jatom in apo:
            distance = np.linalg.norm(iatom.position - jatom.position)

            if distance <= cutoff and jatom.resname in S_aromatic_resname:
                #print(f'{iatom.name}-{iatom.resname}, {jatom.name}-{jatom.resname}{jatom.resid}, {distance}')
                resids.append(jatom.resid)
                S_cryptic.append(f'{jatom.resname}{jatom.resid}')
    S_cryptic = set(S_cryptic)
# -- a set of aromatic residue's names are generated here. note that this is specialised for aromatic residues
    S_all_aroma   = set([f'{residue.resname}{residue.resid}' for residue in holo.residues if residue.resname in S_aromatic_resname])
    S_not_cryptic = S_all_aroma - S_cryptic

    return set(S_cryptic), set(S_not_cryptic)

def main():
    # -- input standardized SASA
    f1 = open('st_dict_sasa.pkl','rb')
    standardized_SASA = pickle.load(f1)

    # -- sets of residues are generated
    S_candidate, S_other = predictCrypticResidues(standardized_SASA)
    print(f'Candidates: {S_candidate}')
    print(f'Others    : {S_other}')

    # -- two holo for bcl-xl are considered to get the union of two residue sets.
    apo_pdb = 'template.pdb'
    holo_pdb1 , ligand1 ='3zlr_A.pdb', 'X0B'
    holo_pdb2 , ligand2 ='2yxj_A.pdb', 'N3C'
    S1_cryptic, S1_not_cryptic = classifyResiduesIntoTwo(apo_pdb, holo_pdb1, ligand1, cutoff=3.)
    S2_cryptic, S2_not_cryptic = classifyResiduesIntoTwo(apo_pdb, holo_pdb2, ligand2, cutoff=3.)
    S_cryptic     =  S1_cryptic     | S2_cryptic # -- new set from both
    S_not_cryptic =  S1_not_cryptic | S2_not_cryptic
    print(f'Cryptic     : {S_cryptic}')
    print(f'Not cryptic : {S_not_cryptic}')

    #-- common set of residues is yielded
    N11 = len(S_candidate & S_cryptic)
    N12 = len(S_candidate & S_not_cryptic)
    N21 = len(S_other     & S_cryptic)
    N22 = len(S_other     & S_not_cryptic)
    print(N11, N12, N21, N22)
    contingency_table = [[ N11,N21],[N12,N22]]

    print('--- Contingency table -- ')
    print(f'              Varied   Not Varied')
    print(f'Cryptic     | {contingency_table[0][0]}        {contingency_table[0][1]}')
    print(f'not Cryptic | {contingency_table[1][0]}        {contingency_table[1][1]}')

    # -- fisher's exact test here
    oddsratio, pvalue = stats.fisher_exact(contingency_table)
    print(f'p-value = {pvalue}')

    plt.show()

if __name__=='__main__':
    main()
