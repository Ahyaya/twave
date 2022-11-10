import twave
import numpy as np

def thermo_stat(Is, rcScale=[1.0,1.0]):
    rf=rcScale[0]
    Rs=0.010/11.5+0.02/(0.0054/rf)
    alpha=0.64
    epss=0.79
    epsgap=0.025633
    epsv=0.0095692
    epsin=0.05
    Th=293.15
    qs_0=0.0
    qs_1=alpha*Is

    for ptr in range(4):
        qs=np.linspace(qs_0,qs_1,64)
        ts=((alpha*Is-qs)/epss/(5.67e-8))**0.25
        tRs=ts-Rs*qs
        tRs[np.where(tRs<0)]=0
        td4=((epsv+epsgap)*tRs**4-qs/5.67e-8)/epsgap
        td4[np.where(td4<0)]=0
        td=td4**0.25
        qd=5.67e-8*(epsgap*(tRs**4-td4)-epsv*td4)
        func=qd-5.67e-8*epsin*(td4-Th**4)
        pivot=np.min(np.where(func>0))
        if(pivot<1):
            pivot=1
        qs_1=qs[pivot]
        qs_0=qs[pivot-1]

    Qs=np.mean(qs[pivot-1:pivot])
    Ts=((alpha*Is-Qs)/epss/(5.67e-8))**0.25
    TRs=Ts-Rs*Qs
    Td4=((epsv+epsgap)*TRs**4-Qs/5.67e-8)/epsgap
    Td=Td4**0.25
    Qd=5.67e-8*(epsgap*(TRs**4-Td4)-epsv*Td4)
    Tb=TRs
    #Func=Qd-5.67e-8*epsin*(Td4-Th**4)
    return Ts, Tb, Td

def thermo_dyn(freq, Is, rcScale=[1.0,1.0]):
    rf=rcScale[0]
    cf=rcScale[1]
    btype = (twave.c_int*3)()
    bprop = (twave.twaveBoundaryProp_t*3)()
    lprop = (twave.twaveLayerProp_t*3)()

    Ts, Tb, Td = thermo_stat(Is,rcScale)
    
    lprop[0].conductivity=11.5
    lprop[0].density=50
    lprop[0].specHeat=945
    lprop[0].thickness=0.010
    
    btype[0]=0
    
    lprop[1].conductivity=5.4e-3/rf
    lprop[1].density=64
    lprop[1].specHeat=1090*cf
    lprop[1].thickness=0.020
    
    btype[1]=2
    bprop[1].radiTransCoef=0.025633
    bprop[1].leakageCoefs[0]=0.0095692
    bprop[1].leakageCoefs[1]=0.0095692
    bprop[1].temperatures[0]=Tb
    bprop[1].temperatures[1]=Td
    
    lprop[2].conductivity=11.5
    lprop[2].density=50
    lprop[2].specHeat=945
    lprop[2].thickness=0.020
    
    btype[2]=2
    bprop[2].radiTransCoef=0.05
    bprop[2].temperatures[0]=293.15

    twave.solve_twave_glvar(freq, lprop, btype, bprop,3)

    h0c2lf = twave.complex2lf_t()
    h1c2lf = twave.complex2lf_t()
    h0bc2lf = twave.complex2lf_t()
    h1bc2lf = twave.complex2lf_t()

    h0c2lf = twave.get_transfunc_glvar_c2lf(0, 2, 0)
    h1c2lf = twave.get_transfunc_glvar_c2lf(0, 2, 1)
    h0bc2lf = twave.get_transfunc_glvar_c2lf(0, 1, 0)
    h1bc2lf = twave.get_transfunc_glvar_c2lf(0, 1, 1)

    h0 = complex(h0c2lf.real, h0c2lf.imag)
    h1 = complex(h1c2lf.real, h1c2lf.imag)
    h0b = complex(h0bc2lf.real, h0bc2lf.imag)
    h1b = complex(h1bc2lf.real, h1bc2lf.imag)

    alpha=0.64
    epss=0.79
    
    hsd=np.abs(alpha/(1.0/h1+4*5.67e-8*epss*(Ts**3)/h0))
    dTd=0.175*(freq**(-1.0/3))*hsd
    dTs=dTd/np.abs(h0)
    hsb=np.abs(alpha/(1/h1b+4*5.67e-8*epss*Ts**3/h0b))
    dTb=0.175*freq**(-1.0/3)*hsb
    return hsb, dTs, dTb, dTd


def sunshield_stat(Is, rcScale=[1.0,1.0]):
    rf=rcScale[0]
    Rs=0.010/11.5+0.02/(0.0054/rf)
    alpha=0.64
    epss=0.79
    epsb=0.05

    qs_0=0.0
    qs_1=alpha*Is

    for ptr in range(4):
        qs=np.linspace(qs_0,qs_1,64)
        ts=((alpha*Is-qs)/epss/(5.67e-8))**0.25
        tb=ts-Rs*qs
        tb[np.where(tb<0)]=0
        func=qs-(5.67e-8)*epsb*(tb**4)
        pivot=np.min(np.where(func>0))
        if(pivot<1):
            pivot=1
        qs_1=qs[pivot]
        qs_0=qs[pivot-1]

    Qs=np.mean(qs[pivot-1:pivot])
    Ts=((alpha*Is-Qs)/epss/(5.67e-8))**0.25
    Tb=Ts-Rs*Qs
    return Ts, Tb

def sunshield_dyn(freq, Is, rcScale=[1.0,1.0]):
    rf=rcScale[0]
    cf=rcScale[1]
    btype = (twave.c_int*2)()
    bprop = (twave.twaveBoundaryProp_t*2)()
    lprop = (twave.twaveLayerProp_t*2)()

    Ts, Tb = sunshield_stat(Is, rcScale)
    
    lprop[0].conductivity=11.5
    lprop[0].density=50
    lprop[0].specHeat=945
    lprop[0].thickness=0.010
    
    btype[0]=0
    
    lprop[1].conductivity=5.4e-3/rf
    lprop[1].density=64
    lprop[1].specHeat=1090*cf
    lprop[1].thickness=0.020
    
    btype[1]=2
    bprop[1].radiTransCoef=0.05
    bprop[1].temperatures[0]=Tb

    twave.solve_twave_glvar(freq, lprop, btype, bprop,2)

    h0c2lf = twave.complex2lf_t()
    h1c2lf = twave.complex2lf_t()

    h0c2lf = twave.get_transfunc_glvar_c2lf(0, 1, 0)
    h1c2lf = twave.get_transfunc_glvar_c2lf(0, 1, 1)

    h0 = complex(h0c2lf.real, h0c2lf.imag)
    h1 = complex(h1c2lf.real, h1c2lf.imag)

    alpha=0.64
    epss=0.79
    
    hsb=np.abs(alpha/(1.0/h1+4*5.67e-8*epss*(Ts**3)/h0))
    dTb=0.175*(freq**(-1.0/3))*hsb
    dTs=dTb/np.abs(h0)
    
    return hsb, dTs, dTb