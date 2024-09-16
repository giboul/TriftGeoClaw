#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import clawpack.geoclaw.topotools as topotools

path = os.path.join("bathymetry.xyz")
topo_file = topotools.Topography(path, topo_type=1)
topo_file.read()

fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)

axes.pcolor(topo_file.X, topo_file.Y, topo_file.Z)

plt.show()
