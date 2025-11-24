from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import scienceplots
plt.style.use("science")


Ea, Ew, ix = np.loadtxt(Path(__file__).parent/"TriftMaxima.csv", delimiter=",").T

plt.figure(layout="tight")
l, = plt.plot(Ea/1e9, Ew/1e9, 'o', mfc="w", label="Vagues d'impulsion (Trift)")
plt.plot(Ea[ix==11]/1e9, Ew[ix==11]/1e9, 'o', mec=l.get_markeredgecolor(), label="Vague de la Figure 4 (Trift)")
plt.axline((0, 0), slope=10/100, ls="-.", alpha=0.5, label="Zitti et al. (2016)")
plt.axline((0, 0), slope=19/100, ls="-.", c="r", alpha=0.5, label="Meng et al. (2019)")
plt.xlabel(r"$E_\mathrm{avalanche}$ [GJ]")
plt.ylabel(r"$E_\mathrm{vague}$ [GJ]")
plt.xlim(0, None)
plt.ylim(0, None)
plt.legend()
plt.show()
