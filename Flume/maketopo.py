#!/usr/bin/env python
xlower = 0
x_A    = 100.
xupper = 400.

slope = 0.05
mx     = 1000 + 1
z_A    = 0-slope*x_A


if __name__ == "__main__":
    import numpy as np

    x = np.linspace(xlower, xupper, num=mx, endpoint=True)
    z = - slope*x
    z[x > x_A] = -slope*x_A

    if False:
        from matplotlib import pyplot as plt
        plt.fill_between(x, z, z.min()-10)
        plt.show()

    np.savetxt(fname="topography.data",
               X=np.vstack((x, z)).T,
               header=f"{mx} # number of points",
               comments="")
