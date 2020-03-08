#!/Users/siida/anaconda3/bin/python
import pickle
import numpy as np
import matplotlib.pyplot as plt

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--sasa_table'   , required=True, help='.pkl')
    p.add_argument('-l', '--residue_list' , nargs='*', required=True, help='e.g. PHE66')
    p.add_argument('-ps', '--plot_style'  , default='plot', help='hist/plot')
    p.add_argument('-b', '--begin'  , type=int, default=0, help='')
    #p.add_argument('-e', '--end'    , type=int, help='')
    args = p.parse_args()
    return args

def draw(std_SASA, label, buried_upper_limit=0.05):
    """
    Args:
        std_SASA: [dict] standardized SASA like {PHE105:ASA, ..., TYR200:ASA}
        buried_upper_limit: [float] upper limit of buriedness (standardized SASA) [no unit]
    Return:
        Reb: Ratio of exposed probability to buried one. Small value indicates buried states are more dominant than exposed.
    """
    bin_start, bin_end, interval = 0, 1, 0.01
    hist = np.histogram(std_SASA, bins=[i for i in np.arange(bin_start,bin_end,interval)])
    probs, bin_edges = hist[0]/np.sum(hist[0]), hist[1]
    half_width       = 0.5 * abs(bin_edges[0] - bin_edges[1])
    hist_xpoints     = np.array([edge+half_width for edge in bin_edges[:-1]])

    plt.plot(hist_xpoints, probs, label=label)
    plt.axvline(x=buried_upper_limit, c='black', linestyle='--')
    plt.xlim(0,1.0)
    plt.ylim(0,0.6)
    plt.xlabel('relative SASA [nm^2]')
    plt.ylabel('Probability')

def main():
    """
    Usage:
        python draw_sasa_distrib.py -f st_dict_sasa.pkl -l PHE66 PHE105

    """
    args = parser()
    sasa_table, residues, style, begin = args.sasa_table, args.residue_list, args.plot_style, args.begin
    f1 = open(sasa_table,'rb')
    sasa    = pickle.load(f1)
    for residue in residues:
        
        n_points = len(sasa[residue][begin:])

        if style == 'hist':
            draw(sasa[residue][begin:], label=residue)

        elif style == 'plot':
            plt.plot(sasa[residue][begin:], label=residue)
            #plt.xticks(np.arange(0,n_points), 
            #           np.arange(begin, begin+n_points) 
            #          )
            plt.ylim(0,1)

    plt.legend()
    plt.tight_layout()

    if style == 'hist':
        plt.savefig('fig_sasa_hist.png', dpi=300)

    elif style == 'plot':
        plt.savefig('fig_sasa_series.png' , dpi=300)

    plt.show()

if __name__ == '__main__':
    main()
