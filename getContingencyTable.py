#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt
import isResiInCryptic as cr

def ratio_pe_pb(st_sasa_resn):
    bin_start, bin_end, interval = 0, 1, 0.02
    hist = np.histogram(st_sasa_resn, bins=[i for i in np.arange(bin_start,bin_end,interval)])
    probs, bin_edges = hist[0], hist[1]
    half_width       = 0.5 * abs(bin_edges[0] - bin_edges[1])
    hist_xpoints = np.array([edge+half_width for edge in bin_edges[:-1]])

    boundary_buried_exposed = 0.1 # boundary for standardised ASA [no unit]

    index_Pe = np.where(hist_xpoints > boundary_buried_exposed)[0]
    index_Pb = np.where(hist_xpoints <= boundary_buried_exposed)[0]
    Pe       = [probs[i] for i in index_Pe]
    xe       = hist_xpoints[hist_xpoints > boundary_buried_exposed]
    Pb       = [probs[j] for j in index_Pb]
    xb       = hist_xpoints[hist_xpoints < boundary_buried_exposed]

    Peb = np.sum(Pe) / np.sum(Pb)

    plt.plot(xe, Pe)
    plt.plot(xb, Pb)
    return Peb

def isASAVaried(Reb):
    # -- Reb thresholds
    too_buried_boundary  = 0.02
    too_exposed_boudnary = 0.7
    if Reb > too_buried_boundary and Reb < too_exposed_boudnary:
        return True

    elif Reb >= too_exposed_boudnary or Reb <= too_buried_boundary:
        return False

def count_two_states_buri_and_expos(st_sasa):
    S_varied    = []
    S_not_varied = []
    for key in st_sasa:
        # -- i just concentrate on aromatic residues.
        # -- the iput pkl file must be modified if you want to check the other residues out, because they've not been normalised.
        if key[0:3] in ['PHE','TRP','TYR','HIS']:
            # -- Calculate the ratio (Reb) of exposed states to buried ones.
            Reb = ratio_pe_pb(st_sasa[key])

            if isASAVaried(Reb):
                S_varied.append(key)
                print(key, 'vaired')

            else:
                S_not_varied.append(key)
                print('not vaired')

    return set(S_varied), set(S_not_varied)

def fishersExactTest():
    pass

def main():
    f1, f2  = open('dict_sasa.pkl','rb'), open('st_dict_sasa.pkl','rb')
    sasa    = pickle.load(f1)
    st_sasa = pickle.load(f2)

    S_candidate, S_other = count_two_states_buri_and_expos(st_sasa)
    print(f'Candidates: {S_candidate}')
    print(f'filtered  : {S_other}')

    pdb, ligand ='3zlr_A.pdb', 'X0B'
    S_cryptic, S_not_cryptic = cr.isResInCryptic(pdb, ligand, cutoff=4.0)
    print(f'Cryptic     : {S_cryptic}')
    print(f'Not cryptic : {S_not_cryptic}')

    print(S_candidate.intersection(S_cryptic))
    print(S_candidate.intersection(S_not_cryptic))

#    draw_some(st_sasa)
    #plt.show()

if __name__=='__main__':
    main()
