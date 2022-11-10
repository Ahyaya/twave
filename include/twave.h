#ifndef TWAVE
#define TWAVE
#endif

struct complex2lf_t
{
    double real;
    double imag;
};

typedef struct complex2lf_t complex2lf_t;

struct twaveSolution_t
{
    double _Complex An;
    double _Complex Bn;
    double mu;
    double TD;
    double freq;
};

typedef struct twaveSolution_t twaveSolution_t;

struct twaveBoundaryProp_t
{
    double radiTransCoef;
    double convectTransCoef;
    double contactResist;
    double ambientEffus;
    double temperatures[2];
    double leakageCoefs[2];
};

typedef struct twaveBoundaryProp_t twaveBoundaryProp_t;

struct twaveLayerProp_t
{
    double conductivity;
    double density;
    double specHeat;
    double thickness;
};

typedef struct twaveLayerProp_t twaveLayerProp_t;

int solve_twave (twaveSolution_t out[],
                double freq,
                const twaveLayerProp_t layerProp[],
                const int boundaryType[],
                const twaveBoundaryProp_t boundaryProp[],
                int length);


double _Complex get_transfunc(twaveSolution_t result[], 
                int inlayer,
                int outlayer, 
                int type);
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


int solve_twave_glvar (double freq,
                const twaveLayerProp_t layerProp[],
                const int boundaryType[],
                const twaveBoundaryProp_t boundaryProp[],
                int length);


double get_transfunc_lf(twaveSolution_t result[], 
                int inlayer,
                int outlayer, 
                int type);


complex2lf_t get_transfunc_glvar_c2lf(int inlayer,
                int outlayer, 
                int type);


double get_transfunc_glvar_lf(int inlayer,
                int outlayer, 
                int type);

