import json
import numpy as np
from matplotlib import pyplot as plt


with open("avalanches.geojson", "r") as file:
    data = json.load(file)

coords = []
for e in data["features"]:
    ix = e['properties']['id']
    points = e['geometry']['coordinates'][0][0]
    for x, y in points:
        coords.append([ix-1, x, y])

coords = np.array(coords)
ix, x, y = coords.T

plt.scatter(x, y, c=ix, cmap="rainbow")
plt.show()
np.savetxt("avalanches.csv", coords)

