import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as mticker
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

def spectro_plot(eflux, time, xlabel=None, ylabel=None,
                 zlabel=None, yscale=None,
                 channels = None, ax=None, figsize=(10,2),
                 vmin=None, vmax=None, lognorm=True, 
                 cmap=None):

    """Function for plotting an energy spectrogram

    Authors:
        Developped by Alexandre Schulz, Vincent Génot, Illya Plotnikov
        Adapted for SOAR SWA-PAS by Rungployphan Kieokaew
        Affiliation: IRAP

    Description: a Python function for plotting an omni-directional energy spectogram
                 using SWA-PAS L2 "1D differential energy flux" as an input

    Arguments:
        eflux: 1D omni-directional energy differential flux given as 2-D array [time, energy]
        time: data timestamp given as 1-D array [time]
        channels: table of energy channels

    Optional arguments
        xlabel: title of X-axis, e.g., Time (UTC)
        ylabel: title of Y-axis, e.g., Energy (eV)
        zlabel: title of the color bar, e.g., ion eflux, eV/cm2-s-eV
        yscale: type of the Y-axis scale, e.g., "log" or "linear"
        ax: axis argument (matplotlib.axes) for plotting
        figsize: figure size
        vmin: minimum flux count value (for data normalization)
        vmax: maximum flux count value (for data normalization)
        lognorm: option for logarithmic color scale - True or False
        cmap: colormap (matplotlib)

    Returns:
        ax object (matplotlib.axes)
    """

    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=figsize)

    X = eflux

    # channels (constant channels case)
    if channels is None:
        m = X.shape[1]#get the number of energy channel
        y = np.arange(0,m,1)
    else:
        y = channels

    # grid
    x1, y1 = np.meshgrid(time, y, indexing="ij")

    # data bounds
    if vmin is None:
        vmin = np.nanmin(X)
    if vmax is None:
        vmax = np.nanmax(X)

    # colormap
    if not cmap:
        cmap = matplotlib.cm.rainbow.copy()
        cmap.set_bad('White',0.)

    # normalize colormapping
    if lognorm and vmin>0.:
        norm=colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm=None

    c = ax.pcolormesh(x1, y1, X, cmap=cmap, norm=norm, edgecolors="face")
    cbar = matplotlib.pyplot.colorbar(c,ax=ax, norm=norm)
    if zlabel:
        cbar.set_label(zlabel)

    if xlabel:
        ax.set_xlabel(xlabel)

    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_ylim(y.min(), y.max())

    if yscale:
        ax.set_yscale(yscale)

    return ax
