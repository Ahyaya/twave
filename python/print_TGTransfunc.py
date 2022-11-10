import tqplatform
import matplotlib.pyplot as plt
import numpy as np

Is = 1200
freq = np.logspace(-4, -2.92, 16)
rf = np.array([1.0, 1.5, 2.0, 3.0])
cf = np.array([1.0, 1.5, 2.0, 3.0])
dTb = (tqplatform.twave.c_double*len(freq))()
dTbc1 = (tqplatform.twave.c_double*len(freq))()
dTbc2 = (tqplatform.twave.c_double*len(freq))()
dTbc3 = (tqplatform.twave.c_double*len(freq))()
dTbr1 = (tqplatform.twave.c_double*len(freq))()
dTbr2 = (tqplatform.twave.c_double*len(freq))()
dTbr3 = (tqplatform.twave.c_double*len(freq))()

freq.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))

for ptr in range(len(freq)):
    dTb[ptr] = tqplatform.thermo_dyn(freq[ptr], Is)[2]
    dTbr1[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [rf[1], 1.0])[2]
    dTbr2[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [rf[2], 1.0])[2]
    dTbr3[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [rf[3], 1.0])[2]
    dTbc1[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [1.0, cf[1]])[2]
    dTbc2[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [1.0, cf[2]])[2]
    dTbc3[ptr] = tqplatform.thermo_dyn(freq[ptr], Is, [1.0, cf[3]])[2]

fig=plt.figure()

ax=plt.gca()
ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in", width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)
ax.plot(1000*freq, dTb, "k-", label="current state")
ax.plot(1000*freq, dTbc1, "^-", markersize=4.5, markerfacecolor="none", label="capacity x1.5", color=plt.cm.Paired(4))
ax.plot(1000*freq, dTbc2, "s-", markersize=4.5, markerfacecolor="none", label="capacity x2.0", color=plt.cm.Paired(2))
ax.plot(1000*freq, dTbc3, "v-", markersize=4.5, markerfacecolor="none", label="capacity x3.0", color=plt.cm.Paired(0))
ax.plot(1000*freq, dTbr1, ":.", label="resistance x1.5", color=plt.cm.Paired(5))
ax.plot(1000*freq, dTbr2, "-.", label="resistance x2.0", color=plt.cm.Paired(3))
ax.plot(1000*freq, dTbr3, "--", label="resistance x3.0", color=plt.cm.Paired(1))

ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_yscale("log")
ax.set_xlim(0,1.35)
ax.set_ylim(3e-5,1)
ax.set_xticks([0.1,0.2,0.3,0.5,0.8,1.0,1.2])

ax.set_xlabel("Frequency [$10^{-3}\,\mathrm{Hz}$]", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Temperature Fluctuation $\delta T_b$ [$\mathrm{K/Hz^{1/2}}$]", fontdict={"weight":"bold","size":12})

ax.legend(loc="upper right", frameon=0, fontsize=8, borderpad=1.2)

#plt.savefig("./figs/TGTransfunc.eps")
plt.show()