import numpy as np
import matplotlib.pyplot as plt
import euchar.utils as u
from scipy.spatial.distance import euclidean

#=================================================

def matplotlib_plot(rows=1, columns=1, figsize=(4,4),
                    font_size=14, facecolor="azure"):
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
        string, default 'azure'. Also RGB colors '#cccccc'.
    
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
    plt.rc('font', size=font_size)
    fig.patch.set_facecolor(facecolor)
    fig.tight_layout()
    return fig, ax

#=================================================

def visualize_triangle_alpha_parametrization(simplex, points, figsize=(3,3)):
    simplex = np.array(simplex)
    if np.sum(simplex != -1) != 3:
        print(f"`simplex` must be a triangle. Your input is {simplex},")
        print(f"which has {np.sum(simplex != -1)} vertices.")
        return None
    
    def plot_circle(ax, C, r, n=100):
        circle = []
        for i in range(n+1):
            x = r * np.cos(2*np.pi / n * i) + C[0]
            y = r * np.sin(2*np.pi / n * i) + C[1]
            circle.append([x,y])
        circle = np.array(circle)
        ax.plot(circle[:,0], circle[:,1], "-r", alpha=0.6, zorder=1)
    
    tri = points[simplex]
    fig, ax = matplotlib_plot(1, 1, figsize=figsize)
    C = u.center_triangle(*tri)
    r = u.parameter_triangle(*tri) / 2
    plot_circle(ax, C, r, n=100)
    ax.scatter(tri[:,0], tri[:,1])
    ax.scatter(C[0], C[1], c="k")
    fig.tight_layout()
    
#=================================================
    
def visualize_tetrahedron_alpha_parametrization(simplex, points, figsize=(3,3)):
    
    def plot_sphere(C, r, tetra=None, n=100, figsize=(5,5)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
    
        print("radius is:", r)
        radius_ok = []
        inside = []
        print("Distances from points to center miniball are: ", end="")
        for point in tetra:
            radius_ok.append(np.isclose(r, euclidean(C, point)))
            inside.append(r+0.00000001 >= euclidean(C, point))
            print(euclidean(C, point), end=" | ")
    
        print("\n\nIs radius equal to distance points to center?", all(radius_ok))
        print("Are points inside?", all(inside))
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
    
    tetra = points[simplex]
    r = u.parameter_tetrahedron(*tetra) / 2
    C = u.center_tetrahedron(*tetra)
    plot_sphere(C, r, tetra)

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
        new_domain.extend([a, (b-a)*0.999999 + a])
        new_curve.extend([curve[i], curve[i]])

    return np.array(new_domain), np.array(new_curve)

#=================================================

def piecewise_constant_surface(bins1, bins2, surface):
    new_bins1, new_bins2 = [], []
    new_surface = np.zeros(shape=(2*len(bins1)-2, 2*len(bins2)-2))

    for i, (a, b) in enumerate(zip(bins1[:-1], bins1[1:])):
        new_bins1.extend([a, (b-a)*0.999999 + a])
    
    for j, (c, d) in enumerate(zip(bins2[:-1], bins2[1:])):
        new_bins2.extend([c, (d-c)*0.999999 + c])
    
    for i, (a, b) in enumerate(zip(bins1[:-1], bins1[1:])):
        for j, (c, d) in enumerate(zip(bins2[:-1], bins2[1:])):
            val = surface[i][j]
            new_surface[2*i][2*j]     = val
            new_surface[2*i+1][2*j]   = val
            new_surface[2*i][2*j+1]   = val
            new_surface[2*i+1][2*j+1] = val
    return np.array(new_bins1), np.array(new_bins2), new_surface
    
#=================================================

def euler_curve_plot(ax, bins, euler_char_curve, line_color="k",
                     line_width=1, marker="", xlim=None, ylim=None,
                     xticks=[0], yticks=[0], xlabel="Parameter",
                     ylabel="$\chi(K)$", title="", figsize=(3,3),
                     font_size=14):
    """
    Display an Euler characteristic curve as a piecewise constant line.
    
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

    ax.plot(new_domain, new_curve, color=line_color,
            linewidth=line_width, marker=marker)
    
    if xlim is not None and ylim is not None:
        ax.set(xticks=xticks, yticks=yticks, xlim=xlim, ylim=ylim,
               xlabel=xlabel, ylabel=ylabel, title=title)
    else:
        ax.set(xticks=xticks, yticks=yticks,
               xlabel=xlabel, ylabel=ylabel, title=title)

#=================================================

def euler_surface_plot(fig, ax, bins1, bins2, euler_char_surf, 
                       n_levels=30, min_level=None, max_level=None,
                       xticks=None, yticks=None, colorbar_ticks=None,
                       xlabel="Parametrization 1",
                       ylabel="Parametrization 2", 
                       dx=0.05, dy=0.05,
                       title="", color_map='RdYlBu',
                       figsize=(3,3), font_size=14):
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
    xlabel, ylabel
        strings, labels of axes
    dx, dy
        floats, resolution parameters of the contour plot. 
        Make these smaller to increase resolution.
    title
        string
    color_map
        string, default is 'RdYlBu'. Other options are: 'Greys', 'PiYG',
        'viridis', 'magma'.
    figsize
       two-tuple, default is (3,3)
    font_size
       int, default is 14

    """
    
    bins1, bins2, euler_char_surf = piecewise_constant_surface(bins1, bins2, euler_char_surf)
    
    # Contours levels
    if min_level is None:
        min_level = euler_char_surf.min()
    if max_level is None:
        max_level = euler_char_surf.max()
    if colorbar_ticks is None:
        colorbar_ticks = [min_level, max_level]
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
        xticks=[min(bins1), max(bins1)]
    if yticks is None:
        yticks=[min(bins2), max(bins2)]
    
    ax.set(title=title, xticks=xticks, yticks=yticks, 
           xlabel=xlabel, ylabel=ylabel)
