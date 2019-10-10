import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#=================================================

def matplotlib_plot(rows=1, columns=1, figsize=(4,4),
                    font_size=14):
    """Returns a fig and ax object to be opearated on."""
    figsize = (figsize[0]*columns, figsize[1]*rows)
    fig, ax = plt.subplots(rows, columns, figsize=figsize)
    plt.rc('font', size=font_size)
    fig.tight_layout()
    return fig, ax

#=================================================

def piecewise_constant_curve(domain, curve):
    try:
        err_line1 = "len of domain and euler_char_curve must be equal\n---"
        assert len(domain) == len(curve), err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None
        
    new_domain = []
    new_curve = []
    for i, (a, b) in enumerate(zip(domain[:-1], domain[1:])):
        new_domain.extend([a, (b-a)*0.9999+a])
        new_curve.extend([curve[i], curve[i]])

    return np.array(new_domain), np.array(new_curve)

#=================================================
def euler_curve_plot(bins, euler_char_curve, line_color="k",
                     line_width=1, marker="", xlim=None, ylim=None,
                     xticks=[0], yticks=[0], xlabel="Parameter",
                     ylabel="$\chi(K)$", title="", figsize=(3,3),
                     font_size=14):
    """Display an Euler characteristic curve as a piecewise constant line.
    Parameters
    ----------
    bins
        np.ndarray, x-axis values of the plot
    euler_char_curve
        np.ndarray of integers
    line_color
        color string, 'k' by default.
    line_width
       int
    marker
        string, can be 'o', '^'
    xlim, ylim
        lists of two int, limits of plot
    xticks, yticks
        lists, ticks to plot along axes
    xlabel, ylabel
        strings, labels of axes
    title
        string
    figsize
       two-tuple, default is (3,3)
    font_size
       int, default is 14

    """
    
    new_domain, new_curve = piecewise_constant_curve(bins, euler_char_curve)

    fig, ax = matplotlib_plot(1, 1, figsize=figsize,
                              font_size=font_size)

    ax.plot(new_domain, new_curve, color=line_color,
            linewidth=line_width, marker=marker)
    
    if xlim is not None and ylim is not None:
        ax.set(xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
               xlabel=xlabel, ylabel=ylabel, title=title)
    else:
        ax.set(xticks=xticks, yticks=yticks,
               xlabel=xlabel, ylabel=ylabel, title=title)

    return fig, ax

#=================================================

def euler_surface_plot(bins1, bins2, euler_char_surf, n_levels=30,
                       dx=0.05, dy = 0.05, title="-", color_map='RdYlBu',
                       figsize=(3,3), font_size=14):
    """Plot Euler characteristic surface."""

    # Contours levels
    levels = MaxNLocator(nbins=n_levels).tick_values(euler_char_surf.min(), euler_char_surf.max())
    # Color map
    cmap = plt.get_cmap(color_map)

    fig, ax = matplotlib_plot(1, 1, figsize=figsize, font_size=font_size)
    xx, yy = np.meshgrid(bins1, bins2)

    # Plot the surface.
    cf = ax.contourf(xx + dx/2.,
                     yy + dy/2.,
                     euler_char_surf,
                     levels=levels,
                     cmap=cmap)
    fig.colorbar(cf, ax=ax)

    #ax.set(xlabel=xlabel, ylabel=ylabel, xticks=xticks, yticks=yticks)
    #ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    #ax.xaxis.set_major_formatter(ticker.FixedFormatter((name_list)))
    ax.set_title(title)

    return fig, ax

#=================================================

def save(path="unnamed.png", dpi=200):
    plt.savefig(fname=path, dpi=dpi)
