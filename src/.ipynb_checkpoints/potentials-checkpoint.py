#defines all relevant potentials
import numpy as np


def step (x_vals):
    return (x_vals >=0).astype(int)

def abs (x_vals):
    return np.abs(x_vals)


