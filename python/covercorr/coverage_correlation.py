import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

'''
Area for rectangles implemented using Segment Tree Class
'''
class SegmentTree:
    def __init__(self, weights):
        # create a segment tree with a given list of weights
        n = len(weights) # number of segments
        m = 2 ** np.ceil(np.log2(n)).astype(int) # number of leaves in a full segement tree
        self.m = m
        self.covered = np.zeros(2 * m, dtype=int) # covering status of each segment
        self.weights = np.zeros(2 * m) # total weights of each segment
        self.scores = np.zeros(2 * m) # weight of covered area in each segment
        self.start = np.zeros(2 * m, dtype=int) # starting point of each segment [start, end)
        self.end = np.zeros(2 * m, dtype=int) # # end point of each segment [start, end)
        
        for node in range(2 * m - 1, 0, -1):
            if node >= m: # leaf node, i.e. corresponding to an elementary segment
                if node - m < n:
                    self.weights[node] = weights[node - m]
                self.start[node] = node - m
                self.end[node] = node - m + 1
            else: # node i has children 2i and 2i+1
                self.weights[node] = self.weights[2 * node] + self.weights[2 * node + 1]
                self.start[node] = self.start[2 * node]
                self.end[node] = self.end[2 * node + 1]
    
    def total_score(self):
        # total length covered by used intervals
        return self.scores[1]
    
    def update(self, l, r, multiplier):
        # update the segment tree by adding in all elementary segments in [l, r), with
        # a `multiplier` weight; if multiplier = -1, it means deletion.
        self._update(1, l, r, multiplier)
        
    def _update(self, node, l, r, multiplier):
        # update the segment tree by adding in all elementary segments in [l, r), with
        # a `multiplier` weight, in descendents of node ; if multiplier = -1, it means deletion.
        
        # update covered counts
        if self.start[node] >= r or self.end[node] <= l: 
            # node outside [l, r), no update needed
            return
        if self.start[node] >= l and self.end[node] <= r:
            # node fully contained in [l, r), update cover counts
            self.covered[node] += multiplier
        else:
            # node partially contained, pass to two children
            self._update(2 * node, l, r, multiplier)
            self._update(2 * node + 1, l, r, multiplier)
            
        # update scores
        if self.covered[node] != 0:
            self.scores[node] = self.weights[node] # entire segment is covered
        elif node >= self.m:
            self.scores[node] = 0 # uncovered leaf node has 0 score
        else:
            # score of uncovered nonleaf node depends on its children
            self.scores[node] = self.scores[2 * node] + self.scores[2 * node + 1] 
    
    def display(self):
        print(pd.DataFrame({'start': self.start, 'end': self.end, 'weights': self.weights,
                            'covered': self.covered, 'scores': self.scores}))

'''
Compute the covered area of rectangles:
It computes the area of union of rectanglges using a sweeping line algorithm.
Rectangles are defined by min and max x and y values.
'''
def covered_area(xmin, xmax, ymin, ymax):
    n = len(xmin)
    
    # Sorting x and y coordinates and ranking them
    x_combined = np.concatenate((xmin, xmax))
    x_order = np.argsort(x_combined)
    x_sorted = x_combined[x_order]
    
    y_combined = np.concatenate((ymin, ymax))
    y_sorted = np.sort(y_combined)
    rank_y = np.argsort(np.argsort(y_combined))
    rank_ymin = rank_y[:n]
    rank_ymax = rank_y[n:]
    
    # create a segment tree with elementary segments defined by y_combined
    ele_seg_len = np.append(np.diff(y_sorted), 0)
    seg_tree = SegmentTree(ele_seg_len)
    
    # create an event data frame indexed by sorted x
    event_type = np.repeat([1, -1], n)
    ymin_repeated = np.tile(ymin, 2)
    ymax_repeated = np.tile(ymax, 2)
    rank_ymin_repeated = np.tile(rank_ymin, 2)
    rank_ymax_repeated = np.tile(rank_ymax, 2)
    
    sweep_events = {
        'x': x_sorted,
        'type': event_type[x_order],
        'ymin': ymin_repeated[x_order],
        'ymax': ymax_repeated[x_order],
        'rank_ymin': rank_ymin_repeated[x_order],
        'rank_ymax': rank_ymax_repeated[x_order],
        'x_width': np.append(np.diff(x_sorted), 0),
        'y_intercept_len': np.zeros(len(x_sorted))
    }
    
#    print(pd.DataFrame(sweep_events))
    
    # iterate through events and compute y intercept length
    for i in range(len(x_sorted) - 1):
        multiplier = sweep_events['type'][i]
        l = sweep_events['rank_ymin'][i]
        r = sweep_events['rank_ymax'][i]
        #print(l, r, multiplier)
        seg_tree.update(l, r, multiplier)
        sweep_events['y_intercept_len'][i] = seg_tree.total_score()
        #seg_tree.display()
    
    #print(pd.DataFrame(sweep_events))
    # compute total area
    total_area = np.sum(sweep_events['y_intercept_len'] * sweep_events['x_width'])
    return total_area
    
'''
Compute the covered volume of hyperrectangles. 
zmin and zmax are each an n x d np.array, where each row is a vector of min/max coordinates for a rectangle
'''
def covered_volume(zmin, zmax):
    if zmin.ndim == 1:
        zmin = zmin[:, np.newaxis]
        zmax = zmax[:, np.newaxis]
        
    n, d = zmin.shape
    if d == 1:
        return covered_area(zmin[:,0], zmax[:,0], np.repeat(0, n), np.repeat(1, n))
    if d == 2:
        return covered_area(zmin[:,0], zmax[:,0], zmin[:, 1], zmax[:, 1])
    
    # for d >= 3, recursively take out one coordinate to perform a sweeping algorithm
    xmin = zmin[:, 0]
    ymin = zmin[:, 1:]
    xmax = zmax[:, 0]
    ymax = zmax[:, 1:]
    
    active = np.repeat(False, n) # keeps track of which hyperrectangles intersect the sweep hyperplane
    
    x_combined = np.concatenate((xmin, xmax))
    x_sorted = np.sort(x_combined)
    x_widths = np.append(np.diff(x_sorted), 0)
    order_x = np.argsort(x_combined)
    multiplier_x = 1 - order_x // n * 2 
    order_x = order_x % n
    
    total_volume = 0
    for i, (order, multiplier) in enumerate(zip(order_x, multiplier_x)):
        flag = multiplier == 1
        active[order] = flag
        slice_area = covered_volume(ymin[active, :], ymax[active, :])
        total_volume += slice_area * x_widths[i]
        
    return total_volume
    
'''
Compute covered volume by partitioning [0,1]^d into b^d blocks, where
b = max integer such that b**d <= n (n = number of rectangles).
For each block, clip rectangles to the block and call covered_volume on the clipped set.
This can accelerate computation when rectangles are small/mostly disjoint.
zmin, zmax: (n x d) arrays with 0 <= zmin < zmax <= 1 assumed.
'''
def covered_volume_partitioned(zmin, zmax):
    if zmin.ndim == 1:
        zmin = zmin[:, np.newaxis]
        zmax = zmax[:, np.newaxis]
    n, d = zmin.shape 
    b = int(np.floor(n ** (1.0 / d)))
    cell_size = 1.0 / b
    total = 0.0
    for idx in itertools.product(range(b), repeat=d):
        block_lo = np.array(idx, dtype=float) * cell_size
        block_hi = block_lo + cell_size
        # Intersect test: zmin < block_hi and zmax > block_lo in all coords
        mask = np.all(zmin < block_hi, axis=1) & np.all(zmax > block_lo, axis=1)
        if not np.any(mask):
            continue
        zmin_sel = zmin[mask].copy()
        zmax_sel = zmax[mask].copy()
        # Clip to block
        np.maximum(zmin_sel, block_lo, out=zmin_sel)
        np.minimum(zmax_sel, block_hi, out=zmax_sel)
        # Remove degenerates
        valid = np.all(zmin_sel < zmax_sel, axis=1)
        if not np.any(valid):
            continue
        total += covered_volume(zmin_sel[valid], zmax_sel[valid])
    return total
    
    
'''
transforming an nxd data matrix x into multivariate rank statistics. This is the 2-Wasserstein distance
optimal transport image of x to a uniform sample in [0,1]^d. Output is an nxd matrix of multivariate ranks
'''
def rank_transform(x):
    n, d = x.shape
    if d == 1: # handle 1d rank separately; 1d optimal transport solution is sorting via rearrangement inequality
        u = np.sort(np.random.rand(n))
        x_rank = np.argsort(np.argsort(x[:, 0]))
        return(u[x_rank][:, np.newaxis])
    
    # for d >= 2, optimal transport solution is computed using the Hungarian algorithm via 
    # scipy.optimize.linear_sum_assignment() with squared Euclidean distance costs
    u = np.random.rand(n, d)
    dist = cdist(x, u, 'sqeuclidean')
    _, opt_assign = linear_sum_assignment(dist)
    return u[opt_assign, :]
    
'''
Function to visualise rectangles and draw the [0,1]^2 square
'''
def plot_rectangles(xmin, xmax, ymin, ymax):
    fig, ax = plt.subplots()
    for (a, c, b, d) in zip(xmin, ymin, xmax, ymax):
        rect = plt.Rectangle((a, c), b - a, d - c, linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    box = plt.Rectangle((0,0), 1, 1, linewidth=2, edgecolor='k', facecolor='none')
    ax.add_patch(box)

    ax.autoscale_view()
    ax.set_aspect('equal', 'box')
    plt.grid(True)
    plt.show()

'''
Split rectangles in [-1,2]^d into smaller rectangles by considering the quotient space R^d / Z^d, 
i.e. rectangles are regarded to reside in [0,1]^d with periodic edges
'''
def split_rectangles(zmin, zmax):
    _, d = zmin.shape
    all_shifts = list(itertools.product([-1,0,1], repeat=d)) # compute all 3^d shifts
    
    # create a dummy row to allow concatenation
    zmin_splitted = np.zeros((1, d))
    zmax_splitted = np.zeros((1, d))
    
    for shift in all_shifts:
        # shift the rectangles
        zmin_shifted = zmin + shift
        zmax_shifted = zmax + shift
        
        # compute intersection with [0,1]^d
        zmin_shifted = np.maximum(zmin_shifted, 0)
        zmax_shifted = np.minimum(zmax_shifted, 1)
        
        idx = np.all(zmin_shifted < zmax_shifted, axis=1)
        zmin_splitted = np.vstack((zmin_splitted, zmin_shifted[idx, :]))
        zmax_splitted = np.vstack((zmax_splitted, zmax_shifted[idx, :]))
        
    # remove the dummy row
    zmin_splitted = zmin_splitted[1:, :]
    zmax_splitted = zmax_splitted[1:, :]
    
    return zmin_splitted, zmax_splitted


'''
exact formula for variance of test stat
'''

def variance_formula(n, d):
    V = ((1 - 2/n)**n - (1 - 1/n)**(2*n)) * n
    multiplier = 1.
    inv_factorial = 1.
    for s in range(1, n+1):
        multiplier *= (1 - (s-1)/n)
        inv_factorial /= s
        V += multiplier * inv_factorial * ((1 - 2/n)**(n-s)) * ((2/(s+1))**d)
    return V

'''
The main function to compute the coverage correlation statistics
'''

def coverage_correlation(x, y, visualise=False):
    if x.ndim == 1: 
        x = x[:, np.newaxis]
    if y.ndim == 1:
        y = y[:, np.newaxis]
    n, dx = x.shape
    ny, dy = y.shape
    assert n == ny, 'x and y should have the same number of observations'
    d = dx + dy
    
    # rank transform x and y
    x_rank = rank_transform(x)
    y_rank = rank_transform(y)
    
    # generate hyperrectangles centred at data points
    eps_x = np.power(n, -1/d) / 2
    eps_y = np.power(n, -1/d) / 2
    zmin = np.hstack((x_rank - eps_x, y_rank - eps_y))
    zmax = np.hstack((x_rank + eps_x, y_rank + eps_y))
        
    # hyperrectangles should 'wrap around' the edges
    zmin_splitted, zmax_splitted = split_rectangles(zmin, zmax)
    
    # to visualise, plot first coordinate of x against last coordinate of y
    # if visualise:
    #    plot_rectangles(zmin_splitted[:, 0], zmax_splitted[:, 0], zmin_splitted[:, -1], zmax_splitted[:, -1])
    
    total_volume = covered_volume_partitioned(zmin_splitted, zmax_splitted)
    excess_vacancy = 1 - np.exp(-1) - total_volume  # 1 - (1 - 1/n)^n is the exact mean
    kappa = excess_vacancy / (1 - np.exp(-1)) # coverage coefficient
    sd = math.sqrt(variance_formula(n, d))  # use variance formula to compute exact variance
    Z = excess_vacancy * math.sqrt(n) / sd  # standardised statistic
    pval = norm.sf(Z)  # Using survival function for upper tail probability to get p-value

    return kappa, pval



