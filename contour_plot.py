import numpy as np
from matplotlib import pyplot as plt
from skimage import io
from yaml import safe_load

with open("config.yaml") as file:
    config = safe_load(file)

image = io.imread("swissALTI3D_merged.tif", as_gray=True)[::-1, ::][5*2500:5*3000, 5*1500:5*2000]


print(image)
plt.imshow(image, origin="lower")
plt.contour(image, levels=(1772, 1776, 1777, 1780, 1785, 1790, 1795))
plt.show()
