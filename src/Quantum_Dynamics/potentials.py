# defines all relevant potentials
import numpy as np

def harmonic(x_vals, curvature=1, deriv=0):
    """
    Returns the potential energy for the harmonic oscillator.

    The potential energy V(x) is defined as .5*curvature*x^2 everywhere.
    The derivatives can also be calculated by adjusting the parameter deriv. 

    Parameters
    ----------
    x_vals:array_like
        The x-value(s) for which the harmonic potential is calculated.
    curvature:float, optional
        The curvature of the harmonic oscillator, set to 1 if unspecified.
    deriv:int, optional
        The function can calculate the n-th order derivative, where n is given
        by deriv,
        
    
    Returns
    -------
    y_vals:array_like
        The harmonic oscillator potential energy for each value of x_vals.
        Can also be the derivatives of the potential energy if deriv is adjusted.
        This is a scalar if x_vals is a scalar.

    """
    if deriv == 0:
         return .5 * curvature * x_vals ** 2
    if deriv == 1:
        return curvature * x_vals
    if deriv == 2:
        return curvature * np.ones(x_vals.size)
    else:
        return np.zeros(x_vals.size)

def step(x_vals, deriv=0):
    """Returns the potential energy for the step function
    
    The potential energy V(x) is defined as 0 for x < 0, and 1 for
    x >= 0. Can also return the derivative if deriv is adjusted; np.Inf
    can be a possible return value.
        
    Parameters
    ----------
    x_vals:array_like
        The x-value(s) for which the step-potential is calculated
    deriv:int, optional
        The nth order derivative of the step potential is calculated, where n 
        is given by deriv. If deriv = 0,(the default behavior)
        the step potential is returned.
    
    Returns
    --------
    y_vals:array_like
        The step function applied to each value of x_vals. Can also be the derivatives
        of the step function if deriv is adjusted. 
        This is a scalar if x_vals is scalar.    
    

    """
    if deriv == 0:
        return np.piecewise(x_vals, [x_vals < 0, x_vals >= 0], [0, 1])
    if deriv == 1:
        return np.where(x_vals == 0, [np.Inf, 0])
    else:
        return np.zeros(x_vals.size)
def abs(x_vals, deriv=0):
    """Returns the potential energy for the absolute value function
    
    The potential energy V(x) is defined as -x for x < 0, and x for
    x >= 0. Can also return derivatives if deriv is adjusted; np.Inf 
    can be a possibe return value. For the first derivative, it will return
    np.nan at x = 0
        
    Parameters
    ----------
    x_vals:array_like
        The x-value(s) for which the absolute value potential is calculated
    deriv;int, optional
        The nth order derivative of the step potential is calculated, where n 
        is given by deriv. If deriv = 0,(the default behavior)
        the absolute value potential is returned.
     
    
    Returns
    --------
    y_vals:array_like
        The absolute value function applied to each value of x_vals. Can also be the
        derivatives if deriv is adjusted.
        This is a scalar if x_vals is scalar.    


    """
    if(deriv == 0):
        return np.piecewise(
                  x_vals,
                 [x_vals < 0, x_vals >= 0],
                 [lambda x: -x, lambda x: x]
                 )

    if(deriv == 1):
        return np.piecewise(
                x_vals,
                [x_vals < 0, x_vals == 0, x_vals > 0],
                [-1, np.nan, 1]
                )
    if(deriv == 2):
        return np.where(x_vals == 0, [np.Inf, 0])
        #it is understood that np.Inf represents 2 * dirac delta
    else:
        return np.zeros(x_vals.size)
def bathtub(x_vals, l_cut=-3, r_cut=3,
            deriv=0):
    

    """Returns the potential energy for a "bathtub" function
    
    The potential energy V(x) is defined as .5(x-l_cut)^2 for x < l_cut, 0 for 
    l_cut <= x < r_cut, and .5(x-r_cut)^2 for x >= r_cut. In other words,
    it is the left and right halves of the normal harmonic oscillator potential
    separated by a region of 0 potential (as in the particle in a box). Example
    can be found in Mazzitelli2017 in American Journal of Chemical Physics.
        
    Parameters
    ----------
    x_vals:array_like
        The x-value(s) for which the absolute value potential is calculated.
    l_cut:float, optional
        The value which defines the cutoff between the left and middle regions.
    r_cut:float, optional
        The value which defines the cutoff between the middle and right regions.
    deriv:int, optional
        The function calculates the nth derivative, where n is given by
        deriv. For deriv = 0, the bathtub function is calculated.
    Returns
    --------
    y_vals:array_like
        The "bathtub" function applied to each value of x_vals, or its derivatives
        if deriv has been changed.
        This is a scalar if x_vals is scalar.    


    """
    if(deriv == 0):
        return np.piecewise(
                x_vals,
                [x_vals < l_cut, (x_vals >= l_cut) & (x_vals < r_cut) , x_vals >= r_cut],
                [lambda x: .5*(x-l_cut)**2, 0, lambda x: .5*(x-r_cut)**2]
                )
    if(deriv == 1):
        return np.piecewise(
                x_vals,
                [x_vals < l_cut, (x_vals >= l_cut) & (x_vals < r_cut) , x_vals >= r_cut],
                [lambda x: x - l_cut, 0, lambda x: x - r_cut]
                )
    if(deriv == 2):
        return np.piecewise(
                x_vals,
                [x_vals < l_cut, (x_vals >=l_cut) & (x_vals < r_cut) , x_vals >= r_cut],
                [1, 0, 1]
                )
def quadratic_piecewise (x_vals, l_cut=-2, r_cut=2, l_curvature=1, c_curvature=1, r_curvature=1):
    """Returns the potential energy for a piecewise-defined quadratic function.
    
    The potential energy V(x) is defined as .5*l_curvature(x-l_cut)^2 for x < l_cut, 0 for 
    l_cut <= x < r_cut, and .5(x-r_cut)^2 for x >= r_cut. In other words,
    it is the left and right halves of the normal harmonic oscillator potential
    separated by a region of 0 potential (as in the particle in a box). Example
    can be found in Mazzitelli2017 in American Journal of Chemical Physics.
        
    Parameters
    ----------
    x_vals:array_like
        The x-value(s) for which the absolute value potential is calculated.
    l_cut:float
        The value which defines the cutoff between the left and middle regions.
    r_cut:float
        The value which defines the cutoff between the middle and right regions.
    
    Returns
    --------
    y_vals:array_like
        The "bathtub" function applied to each value of x_vals.
        This is a scalar if x_vals is scalar.    


    """
    return# np.piecewise(
           # x_vals,
           # [x_vals < l_cut, (x_vals >= l_cut) & (x_vals < r_cut), x_vals >= r_cut]
           # [lambda x: .5*l_curvature
    
