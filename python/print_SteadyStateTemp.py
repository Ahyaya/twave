import tqplatform
import matplotlib.pyplot as plt
import numpy as np

Is = np.linspace(700,1500,16)
freq = np.array([0.1,0.2,0.4,0.8,1.6])*1e-3
Ts = (tqplatform.twave.c_double*len(Is))()
Tb = (tqplatform.twave.c_double*len(Is))()
Td = (tqplatform.twave.c_double*len(Is))()
hs01 = (tqplatform.twave.c_double*len(Is))()
hs02 = (tqplatform.twave.c_double*len(Is))()
hs03 = (tqplatform.twave.c_double*len(Is))()
hs04 = (tqplatform.twave.c_double*len(Is))()
hs05 = (tqplatform.twave.c_double*len(Is))()

Is.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))
freq.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))

for ptr in range(len(Is)):
    ts0, tb0, td0 = tqplatform.thermo_stat(Is[ptr])
    Ts[ptr] = ts0
    Tb[ptr] = tb0
    Td[ptr] = td0
    hs01[ptr] = tqplatform.thermo_dyn(freq[0],Is[ptr])[0]
    hs02[ptr] = tqplatform.thermo_dyn(freq[1],Is[ptr])[0]
    hs03[ptr] = tqplatform.thermo_dyn(freq[2],Is[ptr])[0]
    hs04[ptr] = tqplatform.thermo_dyn(freq[3],Is[ptr])[0]
    hs05[ptr] = tqplatform.thermo_dyn(freq[4],Is[ptr])[0]

fig=plt.figure()

ax1=plt.gca()
ax1.minorticks_on()
ax1.tick_params(which="both",top=1)
ax1.tick_params(which="both",axis="both",direction="in",width=2.0)
ax1.tick_params(which="major",length=6.0)
ax1.tick_params(which="minor",length=4.0)
ax1.plot(Is, hs01, "-", label="$h_s$ @0.1mHz", color=plt.cm.Set2(0))
ax1.plot(Is, hs02, "-.", label="$h_s$ @0.2mHz", color=plt.cm.Set2(1))
ax1.plot(Is, hs03, "o--", markersize=4, markerfacecolor="none", label="$h_s$ @0.4mHz", color=plt.cm.Set2(2))
ax1.plot(Is, hs04, "--", label="$h_s$ @0.8mHz", color=plt.cm.Set2(3))
ax1.plot(Is, hs05, ":.", label="$h_s$ @1.6mHz", color=plt.cm.Set2(4))
ax1.spines["top"].set_linewidth(2)
ax1.spines["right"].set_linewidth(2)
ax1.spines["bottom"].set_linewidth(2)
ax1.spines["left"].set_linewidth(2)
labels = ax1.get_xticklabels()+ax1.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")
ax1.set_yscale("log")

ax2=ax1.twinx()
ax2.minorticks_on()
ax2.tick_params(which="both",axis="both",direction="in",width=2.0)
ax2.tick_params(which="major",length=6.0)
ax2.tick_params(which="minor",length=4.0)
ax2.plot(Is, Ts, "^-", linewidth=2, label="${T_s}$ (sun shield front)", color=plt.cm.Set3(3))
ax2.plot(Is, Tb, "d-", linewidth=2, label="${T_b}$ (sun shield back)", color=plt.cm.Set3(5))
ax2.plot(Is, Td, "s-", linewidth=2, label="${T_d}$ (deck, top plate)", color=plt.cm.Set3(4))
ax2.spines["top"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
labels = ax2.get_xticklabels()+ax2.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax1.set_xlim(600,1600)
ax1.set_ylim(5e-5,1)
ax2.set_xlim(600,1600)
ax2.set_ylim(270,420)

ax1.set_xlabel("Incoming Solar Flux [$\mathrm{W/m^2}$]", fontdict={"weight":"bold","size":12})
ax1.set_ylabel("Transfer Function [$\mathrm{m^2 K/W}$]", fontdict={"weight":"bold","size":12})
ax2.set_ylabel("Equilibrium Temperature [$\mathrm{K}$]", fontdict={"weight":"bold","size":12}, rotation=270, labelpad=16)

ax1.legend(loc="upper left", frameon=0, fontsize=8, borderpad=1.2)
ax2.legend(loc="upper right", frameon=0, fontsize=8, borderpad=1.6)
#plt.savefig("./figs/SteadyStateTemp.eps")
plt.show()