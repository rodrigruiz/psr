'''
adapted from:
 https://jwalton.info/Embed-Publication-Matplotlib-Latex/
 
to find out the width in latex do: 
    \showthe\textwidth
after \begin{document} in latex, 
read latex log message (where you usu\ally read the error messages)
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
 
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27
 
    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    if ratio == 'golden':
        ratio = (5**.5 - 1) / 2
    #if ratio == 'double':
    #    ratio = 2
    #if ratio == 'half':
   #     ratio = 1/2
    #if ratio == 'equal':
    #    ratio = 1
 
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio * ( subplots[0] / subplots[1] )
 
    fig_dim = (fig_width_in, fig_height_in)
 
    return fig_dim
