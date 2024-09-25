import numpy as np


xRes = yRes = 5

lake_level = 1767

xmin = 2669700.
xmax = 2671000.
ymin = 1170100.
ymax = 1172100.

dam_z = lake_level + 10
dam_thk = 30
limit = 2670561

out_format = "binary"

nx = 50
ny = int(1.62*nx)

amr_ratios = dict(
    x=(2,),
    y=(2,),
    t=(2,),
)


def dam_upstream(x, y, l=limit):
    yd = ymax - 0.3*(x-xmin) - 50000/(x-l) - 300
    yd[x > l] = float("inf") 
    return yd


def dam_downstream(x, y, thk=dam_thk):
    d = dam_upstream(x, y)
    u = dam_upstream(x+thk, y, l=limit-thk) + thk
    return np.maximum(u, d)


if __name__=="__main__":
    from matplotlib import pyplot as plt
    x = np.linspace(xmin, xmax)
    plt.plot(x, dam_downstream(x, 0))
    plt.plot(x, dam_upstream(x, 0))
    plt.show()

