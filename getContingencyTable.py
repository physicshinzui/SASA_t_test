#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt

def ratio_pe_pb(st_sasa_resn):
    bin_start, bin_end, interval = 0, 1, 0.02
    hist = np.histogram(st_sasa_resn, bins=[i for i in np.arange(bin_start,bin_end,interval)])
    probs, bin_edges = hist[0], hist[1]
    half_width       = 0.5 * abs(bin_edges[0] - bin_edges[1])
    hist_xpoints = np.array([edge+half_width for edge in bin_edges[:-1]])

    boundary_buried_exposed = 0.1 # boundary for standardised ASA [no unit]
    #exposed_upper_cutoff = 0.2
    #buried_upper_cutoff  = 0.1

    index_Pe = np.where(hist_xpoints > boundary_buried_exposed)[0]
    index_Pb = np.where(hist_xpoints <= boundary_buried_exposed)[0]
    Pe       = [probs[i] for i in index_Pe]
    xe       = hist_xpoints[hist_xpoints > boundary_buried_exposed]
    Pb       = [probs[j] for j in index_Pb]
    xb       = hist_xpoints[hist_xpoints < boundary_buried_exposed]

    Peb = np.sum(Pe) / np.sum(Pb)
    #return Peb

    # making a middle cumu prob
    #index_Pmid1 = [index for index in hist_xpoints if (index not in index_Pe) and (index not in index_Pb) ]
    #m1 = np.where(hist_xpoints <= exposed_upper_cutoff,True, False)
    #m2 = np.where(hist_xpoints >= buried_upper_cutoff ,True, False)
    #masks = m1 * m2
    #Pmid = [prob for mask, prob in zip(masks, probs) if mask]
    #xmid = [x for mask, x in zip(masks, hist_xpoints) if mask]
    #print(masks, Pmid)

    plt.plot(xe, Pe)
    #plt.plot(xmid, Pmid)
    plt.plot(xb, Pb)
#    plt.plot(hist_xpoints, probs)
#    plt.hist(st_sasa_resn)

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
    no_of_varied_sasa    = 0
    no_of_nonvaried_sasa = 0
    for key in st_sasa:
        # -- i just concentrate on aromatic residues.
        # -- the iput pkl file must be modified if you want to check the other residues out, because they've not been normalised.
        if key[0:3] == 'PHE' or key[0:3] == 'TRP' or key[0:3] == 'TYR' or key[0:3] == 'HIS':
            # -- Calculate the ratio (Reb) of exposed states to buried ones.
            Reb = ratio_pe_pb(st_sasa[key])
            #print(key, Reb, isASAVaried(Reb))

            if isASAVaried(Reb):
                no_of_varied_sasa += 1
                print(key, 'vaired')

            else:
                no_of_nonvaried_sasa += 1
                print('not vaired')

    return (no_of_varied_sasa, no_of_nonvaried_sasa)

def draw_some(st_sasa):
    ratio_pe_pb(st_sasa['PHE66'])
    ratio_pe_pb(st_sasa['TRP25'])
    #ratio_pe_pb(st_sasa['PHE105'])
    ratio_pe_pb(st_sasa['TRP98'])
    ratio_pe_pb(st_sasa['HIS138'])
    plt.show()

def main():
    f1, f2  = open('dict_sasa.pkl','rb'), open('st_dict_sasa.pkl','rb')
    sasa    = pickle.load(f1)
    st_sasa = pickle.load(f2)

    ctable = count_two_states_buri_and_expos(st_sasa)
    print(f'The number of candidates: {ctable[0]}, filtered: {ctable[1]}')

#    draw_some(st_sasa)
    plt.show()

if __name__=='__main__':
    main()
