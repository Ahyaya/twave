from ctypes import *

#dynamically load C library libtwave.so as object twave,
libtwave=cdll.LoadLibrary('../lib/libtwave.so')

#Begin of C2py
#=====================================================================
#set a map to translate variable types from C function to Python function
class complex2lf_t(Structure):
    _fields_ = (("real",c_double),("imag",c_double))
class twaveBoundaryProp_t(Structure):
    _fields_ = (("radiDissipateCoef",c_double),("radiTransCoefs",c_double*2),("convectTransCoef",c_double),("contactResist",c_double),("ambientEffus",c_double),("temperatures",c_double*2),("leakageCoefs",c_double*2))
class twaveLayerProp_t(Structure):
    _fields_ = (("conductivity",c_double),("density",c_double),("specHeat",c_double),("thickness",c_double))

libtwave.get_transfunc_glvar_lf.argtypes = (c_int, c_int, c_int)
libtwave.get_transfunc_glvar_lf.restype = c_double
libtwave.get_transfunc_glvar_c2lf.argtypes = (c_int, c_int, c_int)
libtwave.get_transfunc_glvar_c2lf.restype = complex2lf_t
libtwave.solve_twave_glvar.argtypes = (c_double, POINTER(twaveLayerProp_t), POINTER(c_int), POINTER(twaveBoundaryProp_t), c_int)
libtwave.solve_twave_glvar.restype = c_int

solve_twave_glvar = libtwave.solve_twave_glvar
#int solve_twave_glvar (double freq,
#                const twaveLayerProp_t layerProp[],
#                const int boundaryType[],
#                const twaveBoundaryProp_t boundaryProp[],
#                int length);

get_transfunc_glvar_lf = libtwave.get_transfunc_glvar_lf
#double get_transfunc_glvar_lf(int inlayer,
#                int outlayer, 
#                int type);

get_transfunc_glvar_c2lf = libtwave.get_transfunc_glvar_c2lf
#complex2lf_t get_transfunc_glvar_c2lf(int inlayer,
#                int outlayer, 
#                int type);


#End of C2Py
#=====================================================================
