import twave

btype = (twave.c_int*3)()
bprop = (twave.twaveBoundaryProp_t*3)()
lprop = (twave.twaveLayerProp_t*3)()

lprop[0].conductivity=11.5
lprop[0].density=50
lprop[0].specHeat=945
lprop[0].thickness=0.010

btype[0]=0

lprop[1].conductivity=5.4e-3
lprop[1].density=64
lprop[1].specHeat=1090
lprop[1].thickness=0.020

btype[1]=2
bprop[1].radiTransCoef=0.026
bprop[1].leakageCoefs[0]=9.57e-3
bprop[1].leakageCoefs[1]=9.57e-3
bprop[1].temperatures[0]=321.90
bprop[1].temperatures[1]=294.83

lprop[2].conductivity=11.5
lprop[2].density=50
lprop[2].specHeat=945
lprop[2].thickness=0.020

btype[2]=2
bprop[2].radiTransCoef=0.05
bprop[2].temperatures[0]=293.15

freq = 1e-4

twave.solve_twave_glvar(freq, lprop, btype, bprop,3)
h0 = twave.get_transfunc_glvar_lf(0, 2, 0)
h1 = twave.get_transfunc_glvar_lf(0, 2, 1)

#print results to console
print("@1e-4 Hz:\nh0(f)={:.6f}\nh1(f)={:.6f}\n".format(
    h0,
    h1)
)
