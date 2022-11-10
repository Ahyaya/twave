#include <stdio.h>
#include <math.h>
#include <complex.h> 

#include "../include/twave.h"

static const double pi=3.14159265358979;

static const int boundaryType_default[64]={0};

static const twaveBoundaryProp_t boundaryProp_default[64]={0};

static twaveSolution_t glvar_solution[64];

int solve_twave (twaveSolution_t out[],
                double freq,
                const twaveLayerProp_t layerProp[],
                const int boundaryType[],
                const twaveBoundaryProp_t boundaryProp[],
                int length)
{
    int ptr;
    double mu[length];
    double tdepth[length];
    double sqrtPiFreq=sqrt(pi*freq);

    double _Complex An = CMPLX(1,0);
    double _Complex Bn = CMPLX(0,0);
    double _Complex nextAn, nextBn, mu_ambient, qn; 
    double Rn, epsR, L0, L1, T0, T1;

    if (length > 64 || length <= 0){
        return -1;
    }
    if (boundaryType == NULL){
        boundaryType = boundaryType_default;
    }
    if (boundaryProp == NULL){
        boundaryProp = boundaryProp_default;
    }
    if (out == NULL){
        out=glvar_solution;
    }

    switch (boundaryType[length-1])
    {
    case 0:
        mu_ambient = boundaryProp[length-1].ambientEffus;
        break;
    case 1:
        mu_ambient = boundaryProp[length-1].convectTransCoef/sqrt(2*pi*freq)*cexp(-0.25*pi*I);
        break;
    case 2:
        mu_ambient = 4*(5.67e-8)*(boundaryProp[length-1].radiTransCoef)*pow(boundaryProp[length-1].temperatures[0],3)/sqrt(2*pi*freq)*cexp(-0.25*pi*I);
        break;
    default:
        return -3;
    }

    for(ptr=0;ptr<length;++ptr){
        mu[ptr]=sqrt((layerProp[ptr].conductivity)*(layerProp[ptr].density)*(layerProp[ptr].specHeat));
        tdepth[ptr]=(layerProp[ptr].thickness)/sqrt((layerProp[ptr].conductivity)/((layerProp[ptr].density)*(layerProp[ptr].specHeat)));
        out[ptr].mu=mu[ptr];
        out[ptr].TD=tdepth[ptr];
    }
    nextAn = 0.5*cexp((1.0+1.0*I)*sqrtPiFreq*tdepth[length-1])*((1.0 + mu_ambient/mu[length-1])*An+(1.0 - mu_ambient/mu[length-1])*Bn);
    nextBn = 0.5*cexp(-(1.0+1.0*I)*sqrtPiFreq*tdepth[length-1])*((1-mu_ambient/mu[length-1])*An+(1+mu_ambient/mu[length-1])*Bn);

    An=nextAn; Bn=nextBn;
    out[length-1].An=nextAn;
    out[length-1].Bn=nextBn;

    for(ptr=length-2;ptr>=0;--ptr){
        switch (boundaryType[ptr])
        {
        case 0:
            nextAn = 0.5*cexp((1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*((1+mu[ptr+1]/mu[ptr])*An+(1-mu[ptr+1]/mu[ptr])*Bn);
            nextBn = 0.5*cexp(-(1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*((1-mu[ptr+1]/mu[ptr])*An+(1+mu[ptr+1]/mu[ptr])*Bn);
            break;
        case 1:
            Rn = boundaryProp[ptr].contactResist;
            qn = mu[ptr+1]*cexp(0.25*I*pi)*sqrt(2*pi*freq)*(An-Bn);
            nextAn = 0.5*cexp((1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*(Rn*qn+(1+mu[ptr+1]/mu[ptr])*An+(1-mu[ptr+1]/mu[ptr])*Bn);
            nextBn = 0.5*cexp(-(1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*(Rn*qn+(1-mu[ptr+1]/mu[ptr])*An+(1+mu[ptr+1]/mu[ptr])*Bn);
            break;
        case 2:
            L0=boundaryProp[ptr].leakageCoefs[0];
            L1=boundaryProp[ptr].leakageCoefs[1];
            epsR=boundaryProp[ptr].radiTransCoef;
            T0=boundaryProp[ptr].temperatures[0];
            T1=boundaryProp[ptr].temperatures[1];
            nextAn = 0.5*cexp((1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*((cexp(0.25*I*pi)*mu[ptr+1]*sqrt(2*pi*freq)/(4*(5.67e-8)*epsR*pow(T0,3))+(1+L0/epsR)*mu[ptr+1]/mu[ptr])*(An-Bn)+((1+4*(5.67e-8)*L0*pow(T0,3)/(mu[ptr]*sqrt(2*pi*freq))*cexp(-0.25*I*pi))*(1+L1/epsR)*pow(T1/T0,3)+4*(5.67e-8)*L1*pow(T1,3)*cexp(-0.25*I*pi)/(mu[ptr]*sqrt(2*pi*freq)))*(An+Bn));
            nextBn = 0.5*cexp(-(1.0+1.0*I)*sqrtPiFreq*tdepth[ptr])*((cexp(0.25*I*pi)*mu[ptr+1]*sqrt(2*pi*freq)/(4*(5.67e-8)*epsR*pow(T0,3))-(1+L0/epsR)*mu[ptr+1]/mu[ptr])*(An-Bn)+((1-4*(5.67e-8)*L0*pow(T0,3)/(mu[ptr]*sqrt(2*pi*freq))*cexp(-0.25*I*pi))*(1+L1/epsR)*pow(T1/T0,3)-4*(5.67e-8)*L1*pow(T1,3)*cexp(-0.25*I*pi)/(mu[ptr]*sqrt(2*pi*freq)))*(An+Bn));
            break;
        default:
            return -3;
        }
        An=nextAn; Bn=nextBn;
        out[ptr].An=nextAn;
        out[ptr].Bn=nextBn;
    }
    out[0].freq=freq;

    return 0;
}

double _Complex get_transfunc(twaveSolution_t result[], 
                int inlayer,
                int outlayer, 
                int type)
/*
*   type definition: 0 (default)
*   +1 inlayer concerns flux q instead of temperature H
*   +2 outlayer concerns flux q instead of temperature H
*   +4 inlayer is right boundary
*   +8 outlayer is left boundary
*
*   see also the bit definition:
*   bit     3      2      1      0
*   0       R-out  L-in   H-out  H-in
*   1       L-out  R-in   q-out  q-in
*/
{
    double _Complex Ain, Bin, Aout, Bout, wavein, waveout;
    double freq=result[0].freq;
    if(type&0B0100){
        Ain=(result[inlayer].An)*cexp(-(1+I)*sqrt(pi*freq)*(result[inlayer].TD));
        Bin=(result[inlayer].Bn)*cexp((1+I)*sqrt(pi*freq)*(result[inlayer].TD));
    }else{
        Ain=result[inlayer].An;
        Bin=result[inlayer].Bn;
    }
    if(type&0B1000){
        Aout=result[outlayer].An;
        Bout=result[outlayer].Bn;
    }else{
        Aout=(result[outlayer].An)*cexp(-(1+I)*sqrt(pi*freq)*(result[outlayer].TD));
        Bout=(result[outlayer].Bn)*cexp((1+I)*sqrt(pi*freq)*(result[outlayer].TD));
    }
    if(type&0B0001){
        wavein=(result[inlayer].mu)*sqrt(2*pi*freq)*cexp(0.25*pi*I)*(Ain-Bin);
    }else{
        wavein=Ain+Bin;
    }
    if(type&0B0010){
        waveout=(result[outlayer].mu)*sqrt(2*pi*freq)*cexp(0.25*pi*I)*(Aout-Bout);
    }else{
        waveout=Aout+Bout;
    }
    return waveout/wavein;
}

double get_transfunc_lf(twaveSolution_t result[], 
                int inlayer,
                int outlayer, 
                int type)
{
    return cabs(get_transfunc(result, inlayer, outlayer, type));
}

double get_transfunc_glvar_lf(int inlayer,
                int outlayer, 
                int type)
{
    return cabs(get_transfunc(glvar_solution, inlayer, outlayer, type));
}

complex2lf_t get_transfunc_glvar_c2lf(int inlayer, int outlayer, int type)
{
    double _Complex h = get_transfunc(glvar_solution, inlayer, outlayer, type);
    complex2lf_t transfunc={creal(h), cimag(h)};
    return transfunc;
}

int solve_twave_glvar (double freq,
                const twaveLayerProp_t layerProp[],
                const int boundaryType[],
                const twaveBoundaryProp_t boundaryProp[],
                int length)
{
    if(solve_twave(glvar_solution, freq, layerProp, boundaryType, boundaryProp, length)){
        fprintf(stderr,"error call: solve_twave();\n");
        return -1;
    }
    return 0;
}
