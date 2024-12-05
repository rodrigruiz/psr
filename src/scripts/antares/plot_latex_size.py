'''                                                                                                                                                                                                    to find out the width in latex do: \showthe\textwidth
'''
  
def set_size(width, fraction=1, subplots=(1,1), ratio='golden'):
    """Set figure dimensions to avoid scaling in LaTeX.
   
    Parameters
    ----------
       width: float
             Document textwidth or columnwidth in pts
      fraction: float, optional
             Fraction of the width which you wish the figure to occupy  
    Returns
    -------
        fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction
 23  
 24     # Convert from pt to inches
 25     inches_per_pt = 1 / 72.27
 26  
 27     # Golden ratio to set aesthetic figure height
 28     # https://disq.us/p/2940ij3
 29     if ratio == 'golden':
 30         ratio = (5**.5 - 1) / 2
 31     if ratio == 'double':
 32         ratio = 2
 33     if ratio == 'half':
 34         ratio = 1/2
 35     if ratio == 'equal':
 36         ratio = 1
 37  
 38     # Figure width in inches
 39     fig_width_in = fig_width_pt * inches_per_pt
 40     # Figure height in inches
 41     fig_height_in = fig_width_in * ratio * ( subplots[0] / subplots[1] )
 42  
 43     fig_dim = (fig_width_in, fig_height_in)
 44  
 45     return fig_dim
