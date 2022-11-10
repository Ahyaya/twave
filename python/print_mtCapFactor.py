import tqplatform
import matplotlib.pyplot as plt
import numpy as np

Is = 1200
freq = 1e-4
cf = np.linspace(1,7,25)
dTb = (tqplatform.twave.c_double*len(cf))()
dTd = (tqplatform.twave.c_double*len(cf))()

for ptr in range(len(cf)):
    hsb0, dTs0, dTb0, dTd0 = tqplatform.thermo_dyn(freq,Is,[1.0, cf[ptr]])
    dTb[ptr] = dTb0
    dTd[ptr] = dTd0

mt = np.interp(0.01,dTd[::-1],cf[::-1])
print(mt)

plt.figure()

ax=plt.gca()
ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in",width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)
ax.plot(cf, dTb, "-.", label="$\mathrm{\delta} T_b$ (sun shield back)", color=plt.cm.Set1(4))
ax.plot(cf, dTd, "-", label="$\mathrm{\delta} T_d$ (deck, top plate)", color=plt.cm.Set1(1))
ax.plot([2,5],[0.01,0.01], "k--", label="request level")
ax.plot([1,1],[dTb[0],dTd[0]], "o", label="current state", color=plt.cm.Set1(2))
ax.plot([mt,mt],[np.interp(mt,cf,dTb),np.interp(mt,cf,dTd)], "^", label="improved state", color=plt.cm.Set1(3))
ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_xlim(0,8)
ax.set_ylim(0,0.12)

ax.set_xlabel("Scale Factor", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Temperature Fluctuation [$\mathrm{K/Hz^{1/2}}$]", fontdict={"weight":"bold","size":12})
ax.legend(loc="upper right", frameon=0, fontsize=8, borderpad=1.6)
#plt.savefig("./figs/mtCapFactor.eps")
plt.show()