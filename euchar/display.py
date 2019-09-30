import numpy as np
import matplotlib.pyplot as plt

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
def euler_curve_plot(ax_obj, euler_char_curve,
    line_color="k", line_width=1, marker="",
    xlim=[], ylim=[], xticks=[], yticks=[],
    xlabel="Pixel intensity", ylabel="$\chi(K)$", title="",
    label="", no_label=True):
    """Display an Euler characteristic curve as a piecewise constant line."""

    domain = []
    new_curve = []
    for i in range(len(euler_char_curve)):
        domain.extend([i, i+0.999])
        new_curve.extend([euler_char_curve[i], euler_char_curve[i]])
    
    if no_label:
        ax_obj.plot(domain, new_curve, color=line_color,
                    marker=marker, linewidth=line_width)
    else:
        ax_obj.plot(domain, new_curve, color=line_color,
                    marker=marker, linewidth=line_width, label=label)
        
    if len(xlim) != 0 and len(ylim) != 0:
        ax_obj.set(xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
                   xlabel=xlabel, ylabel=ylabel, title=title)
    else:
        ax_obj.set(xticks=xticks, yticks=yticks,
                   xlabel=xlabel, ylabel=ylabel, title=title)
        
#=================================================

def save(path="unnamed.png", dpi=200):
    plt.savefig(fname=path, dpi=dpi)
