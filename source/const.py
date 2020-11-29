"""
const.py

Helpful constants for working with IGS data, and methods to retrieve them.

"""

import numpy as np
import matplotlib.colors as pltcolors

"""
Define constants (e.g. chromosome sizes, colormaps).

"""

SIZES_HG38 = np.array([0, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 
                           159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
                           115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 
                           63025520, 48129895, 51304566, 155270560, 59373566])

SIZES_MM10 = np.array([0, 195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 
                       145441459, 129401213, 124595110, 130694993, 122082543, 120129022,
                       120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
                       61431566, 171031299, 91744698])

CENTROMERES_HG38 = {1:[121535434,124535434], 2:[92326171,95326171], 3:[90504854,93504854],
                    4:[49660117,52660117], 5:[46405641,49405641], 6:[58830166,61830166],
                    7:[58054331,61054331], 8:[43838887,46838887],9:[47367679,50367679],
                    10:[39254935,42254935], 11:[51644205,54644205], 23:[58632012,61632012]}

DISTINCT_COLORS = [(0,0,1), (0,0,1),(1,0,0),(0,1,0),(1,0.10345,0.72414),(1,0.82759,0),(0,0.34483,0),(0.51724,0.51724,1),
          (0.62069,0.31034,0.27586),(0,1,0.75862),(0,0.51724,0.58621),(0,0,0.48276),(0.58621,0.82759,0.31034),
          (0.96552,0.62069,0.86207),(0.82759,0.068966,1),(0.48276,0.10345,0.41379),(0.96552,0.068966,0.37931),
          (1,0.75862,0.51724),(0.13793,0.13793,0.034483),(0.55172,0.65517,0.48276),(0.96552,0.51724,0.034483),
          (0.51724,0.44828,0),(0.44828,0.96552,1),(0.62069,0.75862,1),(0.44828,0.37931,0.48276)]


DISTINCT_CMAP = pltcolors.ListedColormap(DISTINCT_COLORS[1:], name='distinct_cmap')


"""
Helper functions for getting constants.

"""

def get_genome_sizes(genome):
    """
    Helper function to get chromosome sizes for a genome."
    
    Params:
    --------
        genome: genome requested, string referring to assembly
    
    Returns:
    --------
        sizes: array of numbered chromosome sizes for that genome """
    
    if genome == "hg38":
        return SIZES_HG38
    elif genome == "mm10":
        return SIZES_MM10
    else:
        raise ValueError("Genome not found.")

def get_genome_keys(genome):

    """ 
    Helper function to get data table keys for a genome

    Params:
    -------
        genome: genome requested, string referring to assembly
    
    Returns:   
    --------
        keys: dict of keys """

    if genome == "mm10":
        
        keys = {'x':'x_um_abs',
                'y':'y_um_abs',
                'z':'z_um_abs',
                'pos': 'pos',
                'chr':'chr',
                'cluster':'cluster'
                }
    elif genome == "hg38":
        
        keys = {'x':'x_um',
                'y':'y_um',
                'z':'z_um',
                'pos': 'hg38_pos',
                'chr': 'hg38_chr',
                'cluster':'mle_cluster',
                'dim':['x_um', 'y_um','z_um']
                }
        
    else:
        raise ValueError("Genome not found.")
    
    return keys
