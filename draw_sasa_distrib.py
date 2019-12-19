#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt

def draw(std_SASA, buried_upper_limit=0.1):
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

    plt.plot(hist_xpoints, probs)
    plt.axvline(x=buried_upper_limit, c='black')
    plt.xlim(0,1.0)

f1, f2 = open('dict_sasa.pkl','rb'), open('st_dict_sasa.pkl','rb')
sasa    = pickle.load(f1)
st_sasa = pickle.load(f2)

#'PHE97', 'TYR15', 'TYR195', 'TRP137', 'PHE105', 'PHE27'

draw(st_sasa['PHE97'])
#draw(st_sasa['TYR15'])
#draw(st_sasa['TYR195'])
#draw(st_sasa['TRP137'])
draw(st_sasa['PHE105'])
draw(st_sasa['PHE27'])

draw(st_sasa['PHE144'])

plt.show()
