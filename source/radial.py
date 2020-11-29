"""
radial.py

A collection of functions mainly related to manipulating and analyzing the radial
positions of IGS reads in single cells.

"""

import warnings
import numpy as np
import source.const as const
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection


def get_hull(cell, dim=3, genome='hg38'):
    """
    Construct a convex hull of a cell in n-d space.
    
    Params:
    -------
        cell: target cell, dataframe
        dims: number of dimensions for hull
        genome: target genome, for constants
    
    Returns:
    --------
        hull: scipy convex hull object
    """
    
    
    #Init genome-specific dataframe keys
    KEYS = const.get_genome_keys(genome)
    
    #Get spatial position vectors for the clusters
    R = []
    for i in range(dim):
        R.append(cell[KEYS['dim'][i]].values)
    R = np.array(R).T
    
    hull = ConvexHull(R)
    
    return hull

def get_hull_center(hull, dim=3):
    """
    Find the center of a hull's bounding box.

    Params:
    -------
        hull: scipy convhull object
        dim: number of dimensions for hull
    
    Returns:
    --------
        center: bounding box center
    """
    center = np.zeros(dim)
    
    for i in range(dim):
        R_i = hull.points[hull.vertices,i]
        center[i] = (np.max(R_i) + np.min(R_i))/2
        
    return center

def center_cell(cell, origin, dim=3, genome='hg38'):
    """
    Translate cell to a new origin.
    
    Params:
    -------
        cell: target cell, dataframe
        origin: spatial coordinates
        dim: dimensions for translation
        genome: target genome to retrieve constants
    
    Returns:
    --------
        cell: translated cell
    """

    KEYS = const.get_genome_keys(genome)
    
    for index, row in cell.iterrows():
        for i in range(dim):
            cell.at[index, KEYS['dim'][i]] = row[KEYS['dim'][i]] - origin[i]
    
    return cell


def get_r_rel(hull, R):
    """
    Find the relative radial distance of a point relative to the nuclear center
    and periphery. Do this by solving the equation of the plane for the vector 
    defined by the point and all facets of the convex defined by the cell. The
    solution which yields the minimum scaling constant 

    
    Params:
    -------
        hull: scipy convex hull of cell
        R: position of the target point
    
    Returns:
    --------
        r: the relative radial position of the point
        vv: the point on the facet
    """
    
    minC = 1000000 #minimize this
    for (i, e) in enumerate(hull.equations):
        e = e.T
        V,d=e[:-1],e[-1] #plane normal, offset
        if np.dot(V,R) != 0:

            #For scaled point R = C(x1,y1,z1) to be on V = ax + by + cz = - d, 
            #we need s(a*x1 + b*y1 + c*zy1) = - d or s = -d/(V dot R). Minimum
            #C over all candiate planes is smallest scaling factor to get R
            #into a plane.
            
            C = (- d)/np.dot(V,R) #
            if 0 < C and C < minC:
                minC, best_e = C, e

    C = minC
    vv = R*C
    r = np.linalg.norm(R)/np.linalg.norm(vv) #point / point projected into plane
    
    if C == 1000000:
        return -1
    else:
        return(r)


def bootstrap(data, n=5000, func=np.nanmean):
    """
    Generate `n` bootstrap samples, evaluating `func`
    at each resampling. `bootstrap` returns a function,
    which can be called to obtain confidence intervals
    of interest.
    #citation: http://www.jtrive.com/the-empirical-bootstrap-for-confidence-intervals-in-python.html
    """
    simulations = list()
    sample_size = len(data)
    with warnings.catch_warnings(): #catch nanmean runtimewanring
        warnings.simplefilter("ignore", category=RuntimeWarning)
        xbar_init = np.nanmean(data)
    for c in range(n):
        itersample = np.random.choice(data, size=sample_size, replace=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            simulations.append(func(itersample))
    simulations.sort()

    def ci(p):
        """
        Return 2-sided symmetric confidence interval specified
        by p.
        """
        u_pval = (1+p)/2.
        l_pval = (1-u_pval)
        l_indx = int(np.floor(n*l_pval))
        u_indx = int(np.floor(n*u_pval))
        
        return(simulations[l_indx],simulations[u_indx])
    return(ci)

def get_radial_dists(cells,
                     genome='hg38'):
    """
    Get the radial positions of all reads, indexed by chromosome.
    
    Params:
    -------
        cells: a list of the single cells, dataframe
        genome: target genome, to retrieve constants
    
    Returns:
    --------
        R: lists of radial read positions, indexed by chromosome
        
    """
    SIZES = const.get_genome_sizes(genome) #chromosome sizes

    R = [] # to record normed radial distances
    for i in range(len(SIZES)):
        R.append([])

    for cell in cells:
        chr_nums = cell["hg38_chr"].values
        radii = cell["norm_r_2D"].values
    
        for i in range(len(chr_nums)):
            R[chr_nums[i]].append(radii[i])

    return R

def get_radial_statistic(R,
                         func=np.nanmean,
                         genome='hg38'):
    """
    Compute a statistic on the radial data.
   
    Params:
    -------
        R: lists of radial positions, indexed by chromosome
        func: the function used to compute the statistic
        genome: target genome, to retrieve constants
    
    Returns:
    --------
        R_stat: array with a single number per chromosome, i.e. the statistic
    """


    SIZES = const.get_genome_sizes(genome) #chromosome sizes
    
    R_stat = np.zeros(len(SIZES))
    for i in range(len(R)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            R_stat[i] = func(R[i])

    return R_stat
    

def get_radial_CIs(R,R_stat,bounds=0.95, func=np.nanmean):
    """
    Get condience intervals for the radial statistic.
    
    Params:
    -------
        R: the lists of radial measurements indexed by chromosome
        R_stat: the statistic of interest, indexed by chromosome
        bounds: the CI bounds
        func: the function corresponding to the statistic of interest
    Returns:
    --------
        CIs: double sided CIs, indexed by chromosome
    """ 
    CIs = []
    for i in range(len(R)):
        boot = bootstrap(R[i], func=np.nanmean)
        interval = boot(bounds)
        err = np.array([R_stat[i] - interval[0], interval[1] - R_stat[i]])
        CIs.append(err)
    CIs = np.asarray(CIs)

    return CIs

def draw_radial_plot(R_mean, CIs,
                     axbounds=(0.4,0.8),
                     xlabel="Chromosome Size [Mb]",
                     ylabel="Mean rel. radial position"):
    """
    Draw the mean chromosome radial positions as a fn of genomic size.
    
    Params:
    -------
        R_mean: the mean radial chromosome positions
        axbounds, xlabel, ylabel: plot params
    Returns:
    --------
        fig, ax: the plot figure and axes
    """
    sizes = const.SIZES_HG38
    fig = plt.figure()
    ax = fig.add_subplot(111)

    eb1 = ax.errorbar(sizes[1:], R_mean[1:], ls='',yerr=CIs[1:].T, marker='o',
                      capsize=4,markersize=4, markerfacecolor='k', markeredgecolor='k')
    for i in range(1,len(sizes)):

        if i == len(sizes)-1:
            #handle y chromosome edge case
            ax.text(sizes[i], R_mean[i]-1*CIs[i][0], str(i))
        else:
            ax.text(sizes[i], R_mean[i]-2*CIs[i][0], str(i))
            
    ax.set_ylim(axbounds)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    for loc in ['top','right']:
        ax.spines[loc].set_visible(False)

    return fig, ax
