import numpy as np
import matplotlib.pyplot as plt
import euchar.utils as u
from scipy.spatial.distance import euclidean

#=================================================
def matplotlib_plot(rows=1, columns=1, figsize=(4,4), facecolor="w",
                    font_size=10, font_family="CMU Serif"):
    """
    Returns a fig and ax object.
    
    Parameters
    ----------
    rows
        int
    cols
        int
    figsize
        two tuple
    font_size
        int
    facecolor
        string, default 'white'.
    font_family
        string
    
    Returns
    -------
    fig, ax
        matplotlib plots objects. ax will be an array of axes
        if rows>1 and/or cols>1
    
    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> from euchar.display import matplotlib_plot
    >>> from numpy.random import rand
    >>> points1, points2 = rand(10, 2), rand(1000, 2)
    >>> fig, ax = matplotlib_plot(1, 2)
    >>> ax[0].scatter(points1[:,0], points1[:,1], c='red')
    >>> ax[1].scatter(points2[:,0], points2[:,1], c='blue')
    >>> plt.show()
    
    """

    figsize = (figsize[0]*columns, figsize[1]*rows)
    fig, ax = plt.subplots(rows, columns, figsize=figsize)
    
    # Global Setting
    plt.rcParams['font.family'] = 'CMU Serif'
    plt.rcParams['font.size'] = str(font_size)              
    plt.rcParams['mathtext.fontset'] = 'cm'       
    plt.rcParams['mathtext.rm'] = 'CMU Serif'      
    plt.rcParams['grid.color'] = '#F0F0F0'
    plt.rcParams['grid.linestyle'] = 'solid'
    plt.rcParams['text.usetex'] = True   
    
    fig.patch.set_facecolor(facecolor)
    
    return fig, ax

#=================================================

def visualize_triangle_alpha_parametrization(triangle, points, figsize=(3,3)):
    """
    Visualize the miniball of a triangle in the plane.

    Parameters
    ----------

    triangle
       np.ndarray of three integers, the indices of the vertices of the 
        triangle in `points`
    points
        np.ndarray of shape (N, 2), points in 2 dimensions
    figsize
        two-tuple
    
    """
    
    triangle = np.array(triangle)
    if np.sum(triangle != -1) != 3:
        print(f"`triangle` must be a triangle. Your input is {triangle},")
        print(f"which has {np.sum(triangle != -1)} vertices.")
        return None
    
    def plot_circle(ax, C, r, n=100):
        circle = []
        for i in range(n+1):
            x = r * np.cos(2*np.pi / n * i) + C[0]
            y = r * np.sin(2*np.pi / n * i) + C[1]
            circle.append([x,y])
        circle = np.array(circle)
        ax.plot(circle[:,0], circle[:,1], "-r", alpha=0.6, zorder=1)
    
    tri = points[triangle]
    fig, ax = matplotlib_plot(1, 1, figsize=figsize)
    C = u.center_triangle(*tri)
    r = u.parameter_triangle(*tri) / 2
    plot_circle(ax, C, r, n=100)
    ax.scatter(tri[:,0], tri[:,1])
    ax.scatter(C[0], C[1], c="k")
    fig.tight_layout()
    
#=================================================
    
def visualize_tetrahedron_alpha_parametrization(tetrahedron, points, figsize=(3,3)):
    """Visualize the miniball of a tetrahedron in 3D."""
    
    def plot_sphere(C, r, tetra=None, n=100, figsize=(5,5)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
    
        print("radius is:", r)
        radius_ok = []
        inside = []
        print("Distances from points to center miniball are:")
        for point in tetra:
            radius_ok.append(np.isclose(r, euclidean(C, point)))
            inside.append(r+0.00000001 >= euclidean(C, point))
            print(euclidean(C, point), end=" | ")
    
        # Make data
        w = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = r * np.outer(np.cos(w), np.sin(v)) + C[0]
        y = r * np.outer(np.sin(w), np.sin(v)) + C[1]
        z = r * np.outer(np.ones(np.size(w)), np.cos(v)) + C[2]

        # Plot the surface
        ax.plot_surface(x, y, z, color='cyan', alpha=0.4)
    
        # Plot tetra vertices
        if tetra is not None:
            ax.scatter(tetra[:,0], tetra[:,1], tetra[:,2], s = 10, c="r")
        #Plot center
        ax.scatter([C[0]], [C[1]], [C[2]], s = 30, c="k")
    
        ax.set(xlim=[-.3, 1.3], ylim=[-.3, 1.3], zlim=[-.3, 1.3])
    
        fig.tight_layout()
    
    tetra = points[tetrahedron]
    r = u.parameter_tetrahedron(*tetra) / 2
    C = u.center_tetrahedron(*tetra)
    plot_sphere(C, r, tetra)

#=================================================

def piecewise_constant_curve(domain, curve):
    """
    Obtain piecewise constant version of a discrete Euler characteristic
    curve.
    """
    new_domain = []
    new_curve = []
    for i, (a, b) in enumerate(zip(domain[:-1], domain[1:])):
        new_domain.extend([a, (b-a)*0.999999 + a])
        new_curve.extend([curve[i], curve[i]])

    new_domain.append(domain[-1])
    new_curve.append(curve[-1])
    
    return np.array(new_domain), np.array(new_curve)

#=================================================

def piecewise_constant_surface(bins1, bins2, surface):
    """
    Obtain the piecewise constant version of a discrete Euler 
    characteristic surface.
    """
    new_bins1, new_bins2 = [], []
    new_surface = np.zeros(shape=(2*len(bins1)-1, 2*len(bins2)-1))

    for i, (a, b) in enumerate(zip(bins1[:-1], bins1[1:])):
        new_bins1.extend([a, (b-a)*0.999999 + a])
    new_bins1.append(bins1[-1])
    
    for j, (c, d) in enumerate(zip(bins2[:-1], bins2[1:])):
        new_bins2.extend([c, (d-c)*0.999999 + c])
    new_bins2.append(bins2[-1])
        
    for i, (a, b) in enumerate(zip(bins1[:-1], bins1[1:])):
        for j, (c, d) in enumerate(zip(bins2[:-1], bins2[1:])):
            val = surface[i][j]
            new_surface[2*i][2*j]     = val
            new_surface[2*i+1][2*j]   = val
            new_surface[2*i][2*j+1]   = val
            new_surface[2*i+1][2*j+1] = val

    final_row = []
    for item in surface[-1,:-1]:
        final_row.extend([item, item])
    final_row.append(surface[-1,-1])

    final_col = []
    for item in surface[:-1,-1]:
        final_col.extend([item, item])
    final_col.append(surface[-1,-1])

    new_surface[-1] = np.array(final_row)
    new_surface[:,-1] = np.array(final_col)

    return np.array(new_bins1), np.array(new_bins2), new_surface
    
#=================================================

def euler_curve_plot(fig, ax, bins, euler_char_curve, 
                     line_color="k", line_width=2,
                     xlabel="Parameter", ylabel="$\chi(K)$", 
                     title="", 
                     xlim=None, ylim=None,
                     size_arrows=None,
                     xticks=None, yticks=None, 
                     xticks_length=None, xticks_width=None,
                     yticks_length=None, yticks_width=None,
                     xticks_locations=None,
                     yticks_locations=None,
                     font_size_ticks=12):
    """
    Display an Euler characteristic curve as a piecewise constant curve. 
    plot.
    
    Parameters
    ----------
    bins,
        np.ndarray, x-axis values of the plot
    euler_char_curve
        np.ndarray of integers
    line_color
       string, default is "k"
    line_width
       int
    xticks, yticks
        lists, ticks to plot along axes
    xlabel, ylabel
        strings, labels of axes
    title
        string
    xlim, ylim
        two-tuples, limits of the plot
    xticks_locations, yticks_locations
        list of two-tuples
    xticks_spacing, yticks_spacing
        two-tuples, delta of space between tick mark and axis
    x_arrow_head_width, x_arrow_head_length
        floats
    y_arrow_head_width, y_arrow_head_length
        floats

    """
    new_bins, new_curve = piecewise_constant_curve(bins, euler_char_curve)
    shape_pixels = fig.get_size_inches()*fig.dpi
    
    ax.plot(new_bins, new_curve, color=line_color,
            linewidth=line_width, marker="", alpha=1, zorder=3)
    
    # limits and title 
    if xlim is None:
        xlim = [min(bins), max(bins)*1.2]
    if ylim is None:
        ylim = [min(euler_char_curve), max(euler_char_curve)*1.2]
    ax.set(xlim=xlim, ylim=ylim, title=title)
    
    # define proportion arrows
    if size_arrows is None:
        size_arrows = int((shape_pixels[0] + shape_pixels[1])/50)
    
    x_arrow_head_width = (ylim[1]-ylim[0]) / shape_pixels[1] * size_arrows / 1.27
    x_arrow_head_length = (xlim[1] - xlim[0]) / shape_pixels[0] * size_arrows * 1.27
    y_arrow_head_width = (xlim[1]-xlim[0]) / shape_pixels[0] * size_arrows / 1.27
    y_arrow_head_length = (ylim[1] - ylim[0]) / shape_pixels[1] * size_arrows* 1.27
    
    # axes off
    ax.axis("off")
    
    
    # x-axis
    ax.arrow(x=xlim[0], y=0, 
             dx=xlim[1] - xlim[0] - x_arrow_head_length, dy=0,
             head_width=x_arrow_head_width,
             head_length=x_arrow_head_length,
             fc="k", zorder=1)
    
    if xticks_length is None:
        xticks_length = (ylim[1] - ylim[0]) / 30
    if xticks_width is None:
        xticks_width = 1
    if xticks_locations is None:
        xticks_locations = [(tick - y_arrow_head_width * 0.8, 
                             -(1.8 * y_arrow_head_length)) for tick in xticks]
        
    for tick, location in zip(xticks, xticks_locations):
        ax.plot([tick, tick], [0, -xticks_length], color="k", linewidth=xticks_width)
        ax.annotate(str(tick), xy=(0,0), xytext=location,
                    fontsize=font_size_ticks)
        
    # y-axis
    ax.arrow(x=bins[0], y=ylim[0], dx=0, 
             dy=ylim[1] - ylim[0] - y_arrow_head_length,
             head_width=y_arrow_head_width,
             head_length=y_arrow_head_length,
             fc="k", zorder=1)
    if yticks_length is None:
        yticks_length = (xlim[1] - xlim[0]) / 30
    if yticks_width is None:
        yticks_width = 1
    if yticks_locations is None:
        yticks_locations = [(-(1.5 * x_arrow_head_length) * len(str(tick)), 
                             tick - x_arrow_head_width * 0.8) for tick in yticks]
    
    for tick, location in zip(yticks, yticks_locations):
        ax.plot([-yticks_length, 0], [tick, tick],
                 color="k", linewidth=xticks_width)
        ax.annotate(str(tick), xy=(0,0), xytext=location,
                   fontsize=font_size_ticks)

#=================================================

def euler_surface_plot(fig, ax, bins1, bins2, euler_char_surf, 
                       n_levels=30, levels=None,
                       min_level=None, max_level=None,
                       xticks=None, yticks=None, colorbar_ticks=None,
                       xlim=None, ylim=None,
                       xlabel="Parametrization 1",
                       ylabel="Parametrization 2", 
                       dx=0.05, dy=0.05,
                       title="", color_map='coolwarm'):
    """
    Display an Euler characteristic surface as a piecewise constant contour 
    plot.
    
    Parameters
    ----------
    bins1, bin2
        np.ndarray, x-axis and y-axis values of the plot
    euler_char_surface
        np.ndarray of integers
    n_levels
        int, number of levels of the contour plot
    min_levels, max_levels
        floats, minimum and maximum level of the contour plot
    xticks, yticks
        lists, ticks to plot along axes
    colorbar_ticks
        list, ticsk of the colorbar
    xlim, ylim
        two-tuples, limits of the plot
    xlabel, ylabel
        strings, labels of axes
    dx, dy
        floats, resolution parameters of the contour plot. 
        Make these smaller to increase resolution.
    title
        string
    color_map
        string, default is 'coolwarm'. Other options are: 'Greys', 
        'PiYG', 'viridis', 'magma', etc.
  
    """
    
    bins1, bins2, euler_char_surf = piecewise_constant_surface(bins1, bins2, euler_char_surf)
    
    # Contours levels
    if min_level is None:
        min_level = euler_char_surf.min()
    if max_level is None:
        max_level = euler_char_surf.max()
    if colorbar_ticks is None:
        colorbar_ticks = [min_level, max_level]
    if levels is None:
        levels = [int(el) for el in np.linspace(min_level, max_level, num=n_levels)]
    
    # Color map
    cmap = plt.get_cmap(color_map)
    
    # Plot the surface.
    xx, yy = np.meshgrid(bins1, bins2)
    cf = ax.contourf(xx + dx/2., yy + dy/2.,
                     euler_char_surf,
                     levels=levels, cmap=cmap)
    fig.colorbar(cf, ax=ax, ticks=colorbar_ticks)
    
    if xticks is None:
        xticks = [min(bins1), max(bins1)]
    if yticks is None:
        yticks = [min(bins2), max(bins2)]
    if xlim is None:
        xim = [min(bins1), max(bins1)]
    if ylim is None:
        ylim = [min(bins2), max(bins2)]
    
    ax.set(title=title,
           xticks=xticks, yticks=yticks,
           xlim=xlim, ylim=ylim,
           xlabel=xlabel, ylabel=ylabel)
