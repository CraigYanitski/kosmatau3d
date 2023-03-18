import logging
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from pprint import pprint


def plot_box(x=0, y=0, z=0, ds=1, ax=None, set_props=True, **kwargs):
    '''
    This is a short function to plot a box with a shaded surface.
    The input coordinates are the *zero* position of the region
    spanning (x, y, z) to (x+ds, y+dy, z+dz).
    It is useful for representing voxels.
    '''

    if ax == None:
        if 'figsize' in kwargs.keys():
            fig = plt.figure(figsize=kwaargs['figsize'])
        else:
            fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d', elev=20, azim=-140, box_aspect=(1, 1, 1))
    elif not isinstance(ax, Axes3D):
        logging.error(' argument `ax` must be of type `Axes3D`\n'
                      + f'user supplied type {type(ax)}')
        return

    if isinstance(ds, tuple):
        if len(ds) == 3:
            dx = ds[0]
            dy = ds[1]
            dz = ds[2]
        else:
            logging.warn(' argument `ds` must be a float or int for a cube, '
                         + 'or a 3-tuple for an arbitrary rectangle\n'
                         + f'user supplied {ds}.')
            dx = ds[0]
            dy = ds[0]
            dz = ds[0]
    elif isinstance(ds, (int, float)):
        dx = ds
        dy = ds
        dz = ds

    if 'color' in kwargs.keys():
        c = kwargs['color']
        kwargs.pop('color')
    elif 'c' in kwargs.keys():
        c = kwargs['c']
        kwargs.pop('c')
    else:
        c = 'xkcd:electric blue'

    if 'alpha' in kwargs.keys():
        alpha = kwargs['alpha']
        kwargs.pop('alpha')
    else:
        alpha = 0.3

    if 'line_color' in kwargs.keys():
        lc = kwargs['line_color']
        kwargs.pop('line_color')
    elif 'lc' in kwargs.keys():
        lc = kwargs['lc']
        kwargs.pop('lc')
    else:
        lc = 'xkcd:black'

    if 'linewidth' in kwargs.keys():
        lw = kwargs['linewidth']
        kwargs.pop('linewidth')
    elif 'lw' in kwargs.keys():
        lw = kwargs['lw']
        kwargs.pop('lw')
    else:
        lw = 3
    
    if set_props:
        ax.set_proj_type('persp')
        ax.set_axis_off()

    x1, x2 = np.meshgrid(*(2*((0, 1),)))
    surf = np.ones_like(x1)

    ax.plot_surface(  x+dx*x1,   y+dy*x2,    z*surf, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(  x+dx*x1,   y+dy*x2, z+dz*surf, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(  x+dx*x1,    y*surf,   z+dz*x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(  x+dx*x1, y+dy*surf,   z+dz*x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(   x*surf,   y+dy*x1,   z+dz*x2, color=c, alpha=alpha, **kwargs)
    ax.plot_surface(x+dx*surf,   y+dy*x1,   z+dz*x2, color=c, alpha=alpha, **kwargs)

    ax.plot3D(   (x, x+dx),       (y, y),       (z, z), c=lc, lw=lw)
    ax.plot3D(   (x, x+dx), (y+dy, y+dy),       (z, z), c=lc, lw=lw)
    ax.plot3D(   (x, x+dx),       (y, y), (z+dz, z+dz), c=lc, lw=lw)
    ax.plot3D(   (x, x+dx), (y+dy, y+dy), (z+dz, z+dz), c=lc, lw=lw)
    ax.plot3D(      (x, x),    (y, y+dy),       (z, z), c=lc, lw=lw)
    ax.plot3D(      (x, x),    (y, y+dy), (z+dz, z+dz), c=lc, lw=lw)
    ax.plot3D((x+dx, x+dx),    (y, y+dy),       (z, z), c=lc, lw=lw)
    ax.plot3D((x+dx, x+dx),    (y, y+dy), (z+dz, z+dz), c=lc, lw=lw)
    ax.plot3D(      (x, x),       (y, y),    (z, z+dz), c=lc, lw=lw)
    ax.plot3D(      (x, x), (y+dy, y+dy),    (z, z+dz), c=lc, lw=lw)
    ax.plot3D((x+dx, x+dx),       (y, y),    (z, z+dz), c=lc, lw=lw)
    ax.plot3D((x+dx, x+dx), (y+dy, y+dy),    (z, z+dz), c=lc, lw=lw)

    return ax


def plot_clumpy_ism_box(x=0, y=0, z=0, ds=1, N=100, ax=None, set_props=True, **kwargs):
    '''
    This fuction will populate a rectangular region with ISM clumps.
    The input coordinates are the *zero* position of the region
    spanning (x, y, z) to (x+ds, y+dy, z+dz).
    '''

    if ax == None:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d', elev=20, azim=-140, box_aspect=(1, 1, 1))
    elif not isinstance(ax, Axes3D):
        logging.error(' argument `ax` must be of type `Axes3D`\n'
                      + f'user supplied type {type(ax)}')
        return

    if isinstance(ds, tuple):
        if len(ds) == 3:
            dx = ds[0]
            dy = ds[1]
            dz = ds[2]
        else:
            logging.warn(' argument `ds` must be a float or int for a cube, '
                         + 'or a 3-tuple for an arbitrary rectangle\n'
                         + f'user supplied {ds}.')
            dx = ds[0]
            dy = ds[0]
            dz = ds[0]
    elif isinstance(ds, (int, float)):
        dx = ds
        dy = ds
        dz = ds

    if 'color' in kwargs.keys():
        c = kwargs['color']
        kwargs.pop('color')
    elif 'c' in kwargs.keys():
        c = kwargs['c']
        kwargs.pop('c')
    else:
        c = 'xkcd:maroon'

    if 'size' in kwargs.keys():
        s = kwargs['s']
        kwargs.pop('s')
    elif 's' in kwargs.keys():
        s = kwargs['s']
        kwargs.pop('s')
    else:
        s = 3

    if 'surface_alpha' in kwargs.keys():
        sa = kwargs['surface_alpha']
        kwargs.pop('surface_alpha')
    elif 'sa' in kwargs.keys():
        sa = kwargs['sa']
        kwargs.pop('sa')
    else:
        sa = 0.3
    
    if set_props:
        ax.set_proj_type('persp')
        ax.set_axis_off()

    xc, yc, zc = np.random.rand(3, N)

    ax.scatter(x+dx*xc, y+dy*yc, z+dz*zc, color=c, s=s)
    plot_box(x=x, y=y, z=z, ds=ds, ax=ax, alpha=sa, c=c, lc=c, set_props=False)

    return ax
