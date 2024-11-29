n = 50
xmin = -1.5
xmax = 1.5
ymin = -1.5
ymax = 1.5
sea_level = 0.

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    import numpy as np
  
    # (progressing from upper left corner across rows, then down)
    x, y = np.meshgrid(np.linspace(xmin, xmax, n), np.linspace(ymax, ymin, n))
  
    z = np.zeros_like(x)
    z[(0.3*x**2 <= y) & (y <= 0.3*x**2 + 0.2)] = 1 
    z[(y < -1) & (y > -1.2)] = 1
  
    np.savetxt("bathymetry.xyz", np.column_stack((x.flatten(), y.flatten(), z.flatten())))
  
    q0 = np.full_like(z, 0.5)
    q0[y > 0.3*x**2] = 0
    q0[y < -1] = 0
    np.savetxt("qinit.xyz", np.column_stack((x.flatten(), y.flatten(), q0.flatten())))
  
    plt.imshow(z, cmap="inferno")
    q0[q0 <= 0] = float("nan")
    plt.imshow(q0)
    plt.gca().set_xticks((0, n-1), (xmin, xmax))
    plt.gca().set_yticks((0, n-1), (ymax, ymin))
    plt.show()

