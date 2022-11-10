import tqplatform
import matplotlib.pyplot as plt
import numpy as np

Is = 1200
freq = np.logspace(-4, -2.92, 16)

dTs = (tqplatform.twave.c_double*len(freq))()
dTb = (tqplatform.twave.c_double*len(freq))()
dTd = (tqplatform.twave.c_double*len(freq))()

freq.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))

for ptr in range(len(freq)):
    hsb0, dTs0, dTb0, dTd0 = tqplatform.thermo_dyn(freq[ptr],Is)
    dTs[ptr] = dTs0
    dTb[ptr] = dTb0
    dTd[ptr] = dTd0


fig=plt.figure()

ax=plt.gca()
ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in", width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)

ax.plot(1000*freq, dTs, "-s", markersize=4, markerfacecolor="none", label="$\delta T_s$ (sunshield top)", color=plt.cm.Set2(0))
ax.plot(1000*freq, dTb, "-o", markersize=4, markerfacecolor="none", label="$\delta T_b$ (sunshield bottom)", color=plt.cm.Set2(1))
ax.plot(1000*freq, dTd, "-^", markersize=4, markerfacecolor="none", label="$\delta T_d$ (top plate)", color=plt.cm.Set2(2))

ax.set_yscale("log")
ax.set_xticks([0.1,0.2,0.3,0.5,0.8,1.0,1.2])

ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_ylim(1e-5,1)
ax.set_xlabel("Frequency [$10^{-3}\,\mathrm{Hz}$]", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Temperature Response [$\mathrm{K/Hz^{1/2}}$]", fontdict={"weight":"bold","size":12})

ax.legend(loc="lower left", ncol=1, frameon=0, fontsize=8, borderpad=1.2)

#plt.savefig("./figs/TmpRespond.eps")
plt.show()