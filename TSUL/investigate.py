import numpy as np
from matplotlib import pyplot as plt

data = """
   2669891.8750000000        1171551.2500000000     
   2669897.6405274114        1171560.3081827844        0.0000000000000000     
   2669880.1526717558        1171527.8427205102        0.0000000000000000     
   2669870.1596113811        1171567.8002125400        0.0000000000000000     
   2669895.1422623177        1171592.7736450585        0.0000000000000000
"""

P, P1, P2, P3, P4 = [np.array([float(e) for e in line.split(" ") if e]) for line in data.split("\n")[1:-1]]

plt.plot(*np.loadtxt("../TOPM/contour1.xy").T)
plt.plot(*np.loadtxt("../TOPM/contour2.xy").T)
xp, yp, zp = np.vstack((P1, P2, P4, P3)).T
plt.fill(xp, yp)
plt.scatter(*P)
plt.show()
