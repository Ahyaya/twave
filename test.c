#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "include/twave.h"

int main(){
    twaveSolution_t sol[3];
    twaveBoundaryProp_t bprop[3];
    twaveLayerProp_t lprop[3];
    int btype[3];

    lprop[0].conductivity=11.5;
    lprop[0].density=50;
    lprop[0].specHeat=945;
    lprop[0].thickness=0.010;

    btype[0]=0;

    lprop[1].conductivity=5.4e-3;
    lprop[1].density=64;
    lprop[1].specHeat=1090;
    lprop[1].thickness=0.020;

    btype[1]=2;
    bprop[1].radiTransCoef=0.026;
    bprop[1].leakageCoefs[0]=9.57e-3;
    bprop[1].leakageCoefs[1]=9.57e-3;
    bprop[1].temperatures[0]=321.90;
    bprop[1].temperatures[1]=294.83;

    lprop[2].conductivity=11.5;
    lprop[2].density=50;
    lprop[2].specHeat=945;
    lprop[2].thickness=0.020;

    btype[2]=2;
    bprop[2].radiTransCoef=0.05;
    bprop[2].temperatures[0]=293.15;

    double freq = 1e-4;
    double Ts = 360.62;
    double mu0=sqrt((lprop[0].conductivity)*(lprop[0].density)*(lprop[0].specHeat));

    int rc=solve_twave (sol,freq,lprop,btype,bprop,3);
    if(rc){
        fprintf(stderr,"error code: %d\n",rc);
        return rc;
    }
    
    /*
    int ptr;
    for(ptr=0;ptr<3;++ptr){
        fprintf(stderr,"%lf + %lfi   %lf + %lfi\n",
        creal(sol[ptr].An),
        cimag(sol[ptr].An),
        creal(sol[ptr].Bn),
        cimag(sol[ptr].Bn)
        );
    }*/
    
    double complex h0=get_transfunc(sol,0,2,0);
    double complex h1=get_transfunc(sol,0,2,1);
    double complex hs=0.64/(1/h1+4*(5.67e-8)*0.79*pow(Ts,3)/h0);

    fprintf(stderr,"@1e-4 Hz:\n");
    fprintf(stderr,"h0(f) = %lf\n",cabs(h0));
    fprintf(stderr,"h1(f) = %lf\n",cabs(h1));
    fprintf(stderr,"hs(f) = %lf\n",cabs(hs));
    double dI=0.175*pow(freq,-1.0/3.0);
    fprintf(stderr,"\nSCF induced fluctuation:\n%lf %lf K/sqrt{Hz}\n",dI*cabs(hs/h0),dI*cabs(hs));
    
    return 0;
}