import tqplatform
import matplotlib.pyplot as plt
import numpy as np

Is = np.array([900, 1200, 1400])
freq = np.logspace(-4, -2.92, 32)
dTs900 = (tqplatform.twave.c_double*len(freq))()
dTs1200 = (tqplatform.twave.c_double*len(freq))()
dTs1400 = (tqplatform.twave.c_double*len(freq))()
dTb900 = (tqplatform.twave.c_double*len(freq))()
dTb1200 = (tqplatform.twave.c_double*len(freq))()
dTb1400 = (tqplatform.twave.c_double*len(freq))()

Is.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))
freq.ctypes.data_as(tqplatform.twave.POINTER(tqplatform.twave.c_double))

freqsim = np.array([1e-4, 5e-4, 1e-3])
As900 = np.array([0.30925, 0.15726, 0.087472])
Ab900 = np.array([0.091610, 0.024946, 0.0014095])
As1200 = np.array([0.25863, 0.13636, 0.085632])
Ab1200 = np.array([0.075926, 0.021568, 0.0013805])
As1400 = np.array([0.23406, 0.12542, 0.084214])
Ab1400 = np.array([0.068337, 0.019804, 0.001354])

for ptr in range(len(freq)):
    hsb0, dTs0, dTb0 = tqplatform.sunshield_dyn(freq[ptr],Is[0])
    dTs900[ptr] = dTs0
    dTb900[ptr] = dTb0
    hsb0, dTs0, dTb0 = tqplatform.sunshield_dyn(freq[ptr],Is[1])
    dTs1200[ptr] = dTs0
    dTb1200[ptr] = dTb0
    hsb0, dTs0, dTb0 = tqplatform.sunshield_dyn(freq[ptr],Is[2])
    dTs1400[ptr] = dTs0
    dTb1400[ptr] = dTb0

fig=plt.figure()

ax=plt.gca()
ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in", width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)
ax.plot(1000*freq, dTs900, "--", label="$\delta T_s$ ($I_s = 900$)", color=plt.cm.Accent(5))
ax.plot(1000*freq, dTs1200, "-", label="$\delta T_s$ ($I_s = 1200$)", color=plt.cm.Accent(4))
ax.plot(1000*freq, dTs1400, "-.", label="$\delta T_s$ ($I_s = 1400$)", color=plt.cm.Accent(1))
ax.plot(1000*freq, dTb900, "--", label="$\delta T_b$ ($I_s = 900$)", color=plt.cm.Accent(2))
ax.plot(1000*freq, dTb1200, "-", label="$\delta T_b$ ($I_s = 1200$)", color=plt.cm.Accent(0))
ax.plot(1000*freq, dTb1400, "-.", label="$\delta T_b$ ($I_s = 1400$)", color=plt.cm.Accent(6))

ax.plot(1000*freqsim, As900, "v", markersize=4, markerfacecolor="none", label="$\delta {T_s}^G$ ($I_s = 900$)", color=plt.cm.Accent(5))
ax.plot(1000*freqsim, As1200, "o", markersize=4, markerfacecolor="none", label="$\delta {T_s}^G$ ($I_s = 1200$)", color=plt.cm.Accent(4))
ax.plot(1000*freqsim, As1400, "^", markersize=4, markerfacecolor="none", label="$\delta {T_s}^G$ ($I_s = 1400$)", color=plt.cm.Accent(1))
ax.plot(1000*freqsim, Ab900, "v", markersize=4, markerfacecolor="none", label="$\delta {T_b}^G$ ($I_s = 900$)", color=plt.cm.Accent(2))
ax.plot(1000*freqsim, Ab1200, "o", markersize=4, markerfacecolor="none", label="$\delta {T_b}^G$ ($I_s = 1200$)", color=plt.cm.Accent(2))
ax.plot(1000*freqsim, Ab1400, "^", markersize=4, markerfacecolor="none", label="$\delta {T_b}^G$ ($I_s = 1400$)", color=plt.cm.Accent(6))

ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_xlim(0,1.25)
ax.set_ylim(-0.02,0.4)
ax.set_xticks([0.1,0.2,0.3,0.5,0.8,1.0,1.2])

ax.set_xlabel("Frequency [$10^{-3}\,\mathrm{Hz}$]", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Temperature Fluctuation [$\mathrm{K/Hz^{1/2}}$]", fontdict={"weight":"bold","size":12})

ax.legend(loc="upper right", ncol=2, frameon=0, fontsize=8, borderpad=1.2)

#plt.savefig("./figs/SimValid.eps")
plt.show()