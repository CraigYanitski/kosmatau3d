import logging
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from pprint import pprint


def plot_box(x=0, y=0, z=0, ds=1, ax=None, set_proj=False, **kwargs):
    '''
    This is a short function to plot a box with a shaded surface.
    The input coordinates are the *zero* position of the region
    spanning (x, y, z) to (x+ds, y+dy, z+dz).
    It is useful for representing voxels.
    '''

    if ax=None:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d', elev=20, azim=-140)
    elif not isinstance(ax, Axes3D):
        logging.error('ERROR: argument `ax` must be of type `Axes3D`\n'
                      + f'user supplied type {type(ax)}')
        return

    if not ('color' in kwargs.keys() or 'c' in kwargs.keys()):
        c = 'xkcd:electric blue'

    if not 'alpha' in kwargs.keys():
        alpha = 0.3

    if 'line_color' in kwargs.keys():
        lc = line_color
        kwargs.pop('line_color')
    elif 'lc' in kwargs.keys():
        lc = lc
        kwargs.pop('lc')
    else:
        lc = 'xkcd:black'

    if not ('linewidth' in kwargs.keys() or 'lw' in kwargs.keys()):
        lw = 3
    
    if set_props:
        ax.set_proj_type('persp')
        ax.set_axes_off()

    x1, x2 = np.meshgrid(2*((0, ds),))

    ax.plot_surface(x+x1, y+x2,    z, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(x+x1, y+x2, z+ds, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(x+x1,    y, z+x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(x+x1, y*ds, z+x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface   (x, y+x1, z+x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(x+ds, y+x1, z+x2, color=c, alpha=alpha, **kwargs)

    r = np.array([0, 1])
    for s, e in combinations(np.array(list(product(r+ds, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), c=lc, lw=lw)

    return ax


def plot_clumpy_ism_box(x=0, y=0, z=0, ds=1, N=100, ax=None, set_proj=True, **kwargs):
    '''
    This fuction will populate a rectangular region with ISM clumps.
    The input coordinates are the *zero* position of the region
    spanning (x, y, z) to (x+ds, y+dy, z+dz).
    '''

    if ax=None:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d', elev=20, azim=-140)
    elif not isinstance(ax, Axes3D):
        logging.error('ERROR: argument `ax` must be of type `Axes3D`\n'
                      + f'user supplied type {type(ax)}')
        return

    if isinstance(ds, tuple):
        if len(ds) == 3:
            dx = ds[0]
            dy = ds[1]
            dz = ds[2]
        else:
            logging.warn('WARNING: argument `ds` must be a float or int for a cube, '
                         + 'or a 3-tuple for an arbitrary rectangle\n'
                         + f'user supplied {ds}.')
    elif isinstance(ds, (int, float)):
        dx = ds
        dy = ds
        dz = ds

    if not ('color' in kwargs.keys() or 'c' in kwargs.keys()):
        c = 'xkcd:maroon'

    if not ('size' in kwargs.keys() or 's' in kwargs.keys()):
        c = 3

    if 'surface_alpha' in kwargs.keys():
        sa = surface_alpha
        kwargs.pop('surface_alpha')
    elif 'sa' in kwargs.keys()):
        sa = sa
        kwargs.pop('sa')
    else:
        sa = 0.3

    xc, yc, zc = np.random.rand(3, N)

    ax.scatter(x+dx*xc, y+dy*yc, z+dz*zc, color=c, s=s)
    plot_box(x=x, y=y, z=z, ds=ds, ax=ax, alpha=sa, c=c, lc=c, set_proj=False)

    return ax
