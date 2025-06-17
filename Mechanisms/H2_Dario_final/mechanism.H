#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  H
1  N
2  O
3  Ar
4  He
*/

// Species
#define H_ID 0
#define H2_ID 1
#define O_ID 2
#define OH_ID 3
#define H2O_ID 4
#define O2_ID 5
#define N2_ID 6
#define HO2_ID 7
#define H2O2_ID 8
#define OHV_ID 9
#define AR_ID 10
#define HE_ID 11

#define NUM_GAS_ELEMENTS 5// Elements in the homogeneous phase
#define NUM_GAS_SPECIES 12// Species in the homogeneous phase
#define NUM_GAS_REACTIONS 35// Reactions in the homogeneous phase

#define SITE_DENSITY 0.000000E+00// mol/cm^2

#define NUM_SURFACE_ELEMENTS 0// Additional elements in heterogeneous phase
#define NUM_SURFACE_SPECIES 0// Species in the heterogeneous phase
#define NUM_SURFACE_REACTIONS 0// Reactions in the heterogeneous phase

#define NUM_ELEMENTS (NUM_GAS_ELEMENTS + NUM_SURFACE_ELEMENTS)
#define NUM_SPECIES (NUM_GAS_SPECIES + NUM_SURFACE_SPECIES)
#define NUM_REACTIONS (NUM_GAS_REACTIONS + NUM_SURFACE_REACTIONS)

#define NUM_IONS 0

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp 
void atomicWeight(amrex::Real *  awt);
//  MISC 
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int * _rmap);
void CKINU(const int i, int &nspec, int * ki, int * nu);
void CKKFKR(const amrex::Real P, const amrex::Real T,const amrex::Real * x, amrex::Real *  q_f, amrex::Real *  q_r);
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r,amrex::Real *  sc, amrex::Real T);
//  SPARSE INFORMATION 
void SPARSITY_INFO(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, const int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, const int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base);


// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
mm = 5;
kk = 12;
ii = 35;
nfit = -1; // Why do you need this anyway ? 
}

//  inverse molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[12]={
0.9920634920634921,// H
0.4960317460317460,// H2
0.0625039064941559,// O
0.0587993179279120,// OH
0.0555092978073827,// H2O
0.0312519532470779,// O2
0.0356964374955379,// N2
0.0302975216627280,// HO2
0.0293996589639560,// H2O2
0.0587993179279120,// OHV
0.0250312891113892,// AR
0.2498374807187924,// HE
};
#endif
const amrex::Real h_global_imw[12]={
0.9920634920634921,// H
0.4960317460317460,// H2
0.0625039064941559,// O
0.0587993179279120,// OH
0.0555092978073827,// H2O
0.0312519532470779,// O2
0.0356964374955379,// N2
0.0302975216627280,// HO2
0.0293996589639560,// H2O2
0.0587993179279120,// OHV
0.0250312891113892,// AR
0.2498374807187924,// HE
};

//  molecular weights 
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[12]={
1.008000,// H
2.016000,// H2
15.999000,// O
17.007000,// OH
18.015000,// H2O
31.998000,// O2
28.014000,// N2
33.006000,// HO2
34.014000,// H2O2
17.007000,// OHV
39.950000,// AR
4.002602,// HE
};
#endif
const amrex::Real h_global_mw[12]={
1.008000,// H
2.016000,// H2
15.999000,// O
17.007000,// OH
18.015000,// H2O
31.998000,// O2
28.014000,// N2
33.006000,// HO2
34.014000,// H2O2
17.007000,// OHV
39.950000,// AR
4.002602,// HE
};

//  inverse molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
imw_new[0] = 0.9920634920634921;// H
imw_new[1] = 0.4960317460317460;// H2
imw_new[2] = 0.0625039064941559;// O
imw_new[3] = 0.0587993179279120;// OH
imw_new[4] = 0.0555092978073827;// H2O
imw_new[5] = 0.0312519532470779;// O2
imw_new[6] = 0.0356964374955379;// N2
imw_new[7] = 0.0302975216627280;// HO2
imw_new[8] = 0.0293996589639560;// H2O2
imw_new[9] = 0.0587993179279120;// OHV
imw_new[10] = 0.0250312891113892;// AR
imw_new[11] = 0.2498374807187924;// HE
}

//  inverse molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real imw(const int n){
#if AMREX_DEVICE_COMPILE
return global_imw[n];
#else
return h_global_imw[n];
#endif
}
//  molecular weights 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
mw_new[0] = 1.008000;// H
mw_new[1] = 2.016000;// H2
mw_new[2] = 15.999000;// O
mw_new[3] = 17.007000;// OH
mw_new[4] = 18.015000;// H2O
mw_new[5] = 31.998000;// O2
mw_new[6] = 28.014000;// N2
mw_new[7] = 33.006000;// HO2
mw_new[8] = 34.014000;// H2O2
mw_new[9] = 17.007000;// OHV
mw_new[10] = 39.950000;// AR
mw_new[11] = 4.002602;// HE
}

//  molecular weight 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mw(const int n){
#if AMREX_DEVICE_COMPILE
return global_mw[n];
#else
return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
 ru  = 8.31446261815324e+07; 
 ruc = 1.98721558317399615845; 
 pa  = 1.01325e+06; 
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
// species 1: H2
species[1] =
+1.72113605e+00 +4.79214245e-03 * T-1.01482243e-05  * T2+8.96936949e-09 * T3-2.70663102e-12 * T4
;
// species 2: O
species[2] =
+2.10652772e+00 -2.75690396e-03 * T+5.11435812e-06  * T2-4.29494187e-09 * T3+1.34756569e-12 * T4
;
// species 3: OH
species[3] =
+2.98454038e+00 -2.33807954e-03 * T+4.43232593e-06  * T2-3.65814533e-09 * T3+1.27094854e-12 * T4
;
// species 4: H2O
species[4] =
+3.23947087e+00 -2.38193734e-03 * T+7.53145778e-06  * T2-6.70039089e-09 * T3+2.27801570e-12 * T4
;
// species 5: O2
species[5] =
+2.65493298e+00 -1.91768034e-03 * T+6.68974292e-06  * T2-5.89495614e-09 * T3+1.66341808e-12 * T4
;
// species 6: N2
species[6] =
+2.58760040e+00 -6.02547100e-04 * T+8.98331257e-07  * T2+7.54921682e-10 * T3-7.07467604e-13 * T4
;
// species 7: HO2
species[7] =
+2.87761522e+00 -1.15984888e-03 * T+1.06552565e-05  * T2-1.16818377e-08 * T3+4.03564321e-12 * T4
;
// species 8: H2O2
species[8] =
+2.82819116e+00 +3.27308040e-03 * T+5.58298816e-06  * T2-8.21779307e-09 * T3+3.05493365e-12 * T4
;
// species 9: OHV
species[9] =
+2.63725927e+00 +1.85145041e-04 * T-1.67631272e-06  * T2+2.38736757e-09 * T3-8.43208359e-13 * T4
;
// species 10: AR
species[10] =
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
// species 11: HE
species[11] =
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
}
else {
// species 0: H
species[0] =
+1.50000000e+00 
;
// species 1: H2
species[1] =
+1.93286575e+00 +8.26608026e-04 * T-1.46402364e-07  * T2+1.54100414e-11 * T3-6.88804800e-16 * T4
;
// species 2: O
species[2] =
+1.54363697e+00 -2.73162486e-05 * T-4.19029520e-09  * T2+4.95481845e-12 * T3-4.79553694e-16 * T4
;
// species 3: OH
species[3] =
+1.83853033e+00 +1.10741289e-03 * T-2.94000209e-07  * T2+4.20698729e-11 * T3-2.42289890e-15 * T4
;
// species 4: H2O
species[4] =
+1.67703890e+00 +2.97318160e-03 * T-7.73768890e-07  * T2+9.44335140e-11 * T3-4.26899910e-15 * T4
;
// species 5: O2
species[5] =
+2.66096065e+00 +6.56365811e-04 * T-1.41149627e-07  * T2+2.05797935e-11 * T3-1.29913436e-15 * T4
;
// species 6: N2
species[6] =
+1.95257637e+00 +1.39690040e-03 * T-4.92631603e-07  * T2+7.86010195e-11 * T3-4.60755204e-15 * T4
;
// species 7: HO2
species[7] =
+3.17228741e+00 +1.88117627e-03 * T-3.46277286e-07  * T2+1.94657549e-11 * T3+1.76256905e-16 * T4
;
// species 8: H2O2
species[8] =
+3.57977305e+00 +4.05326003e-03 * T-1.29844730e-06  * T2+1.98211400e-10 * T3-1.13968792e-14 * T4
;
// species 9: OHV
species[9] =
+1.88273000e+00 +1.01397430e-03 * T-2.27687700e-07  * T2+2.17468300e-11 * T3-5.12630500e-16 * T4
;
// species 10: AR
species[10] =
+1.50000000e+00 
;
// species 11: HE
species[11] =
+1.50000000e+00 
;
}
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
// species 1: H2
species[1] =
+2.72113605e+00 +4.79214245e-03 * T-1.01482243e-05 * T2+8.96936949e-09 * T3-2.70663102e-12 * T4
;
// species 2: O
species[2] =
+3.10652772e+00 -2.75690396e-03 * T+5.11435812e-06 * T2-4.29494187e-09 * T3+1.34756569e-12 * T4
;
// species 3: OH
species[3] =
+3.98454038e+00 -2.33807954e-03 * T+4.43232593e-06 * T2-3.65814533e-09 * T3+1.27094854e-12 * T4
;
// species 4: H2O
species[4] =
+4.23947087e+00 -2.38193734e-03 * T+7.53145778e-06 * T2-6.70039089e-09 * T3+2.27801570e-12 * T4
;
// species 5: O2
species[5] =
+3.65493298e+00 -1.91768034e-03 * T+6.68974292e-06 * T2-5.89495614e-09 * T3+1.66341808e-12 * T4
;
// species 6: N2
species[6] =
+3.58760040e+00 -6.02547100e-04 * T+8.98331257e-07 * T2+7.54921682e-10 * T3-7.07467604e-13 * T4
;
// species 7: HO2
species[7] =
+3.87761522e+00 -1.15984888e-03 * T+1.06552565e-05 * T2-1.16818377e-08 * T3+4.03564321e-12 * T4
;
// species 8: H2O2
species[8] =
+3.82819116e+00 +3.27308040e-03 * T+5.58298816e-06 * T2-8.21779307e-09 * T3+3.05493365e-12 * T4
;
// species 9: OHV
species[9] =
+3.63725927e+00 +1.85145041e-04 * T-1.67631272e-06 * T2+2.38736757e-09 * T3-8.43208359e-13 * T4
;
// species 10: AR
species[10] =
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
// species 11: HE
species[11] =
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
;
}
else {
// species 0: H
species[0] =
+2.50000000e+00 
;
// species 1: H2
species[1] =
+2.93286575e+00 +8.26608026e-04 * T-1.46402364e-07 * T2+1.54100414e-11 * T3-6.88804800e-16 * T4
;
// species 2: O
species[2] =
+2.54363697e+00 -2.73162486e-05 * T-4.19029520e-09 * T2+4.95481845e-12 * T3-4.79553694e-16 * T4
;
// species 3: OH
species[3] =
+2.83853033e+00 +1.10741289e-03 * T-2.94000209e-07 * T2+4.20698729e-11 * T3-2.42289890e-15 * T4
;
// species 4: H2O
species[4] =
+2.67703890e+00 +2.97318160e-03 * T-7.73768890e-07 * T2+9.44335140e-11 * T3-4.26899910e-15 * T4
;
// species 5: O2
species[5] =
+3.66096065e+00 +6.56365811e-04 * T-1.41149627e-07 * T2+2.05797935e-11 * T3-1.29913436e-15 * T4
;
// species 6: N2
species[6] =
+2.95257637e+00 +1.39690040e-03 * T-4.92631603e-07 * T2+7.86010195e-11 * T3-4.60755204e-15 * T4
;
// species 7: HO2
species[7] =
+4.17228741e+00 +1.88117627e-03 * T-3.46277286e-07 * T2+1.94657549e-11 * T3+1.76256905e-16 * T4
;
// species 8: H2O2
species[8] =
+4.57977305e+00 +4.05326003e-03 * T-1.29844730e-06 * T2+1.98211400e-10 * T3-1.13968792e-14 * T4
;
// species 9: OHV
species[9] =
+2.88273000e+00 +1.01397430e-03 * T-2.27687700e-07 * T2+2.17468300e-11 * T3-5.12630500e-16 * T4
;
// species 10: AR
species[10] =
+2.50000000e+00 
;
// species 11: HE
species[11] =
+2.50000000e+00 
;
}
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+2.547366000000000e+04 * invT+2.946682850000000e+00 -2.500000000000000e+00 * logT-2.740033625000000e-16 * T+2.665623750000000e-19 * T2-1.552516841666667e-22 * T3+3.672898590000000e-26 * T4
;
// species 1: H2
species[1] =
-9.511467020000000e+02 * invT+3.559656969000000e+00 -2.721136050000000e+00 * logT-2.396071225000000e-03 * T+1.691370716666666e-06 * T2-7.474474574999999e-10 * T3+1.353315510000000e-13 * T4
;
// species 2: O
species[2] =
+2.912770090000000e+04 * invT+8.052918099999999e-01 -3.106527720000000e+00 * logT+1.378451980000000e-03 * T-8.523930200000001e-07 * T2+3.579118225000000e-10 * T3-6.737828450000000e-14 * T4
;
// species 3: OH
species[3] =
+3.369554460000000e+03 * invT+4.058480713800000e+00 -3.984540380000000e+00 * logT+1.169039770000000e-03 * T-7.387209883333333e-07 * T2+3.048454441666667e-10 * T3-6.354742700000001e-14 * T4
;
// species 4: H2O
species[4] =
-3.029732530000000e+04 * invT+5.253373550000000e+00 -4.239470870000000e+00 * logT+1.190968670000000e-03 * T-1.255242963333333e-06 * T2+5.583659075000000e-10 * T3-1.139007850000000e-13 * T4
;
// species 5: O2
species[5] =
-1.052703670000000e+03 * invT-5.176797700000000e-01 -3.654932980000000e+00 * logT+9.588401700000001e-04 * T-1.114957153333333e-06 * T2+4.912463450000000e-10 * T3-8.317090400000000e-14 * T4
;
// species 6: N2
species[6] =
-1.051964560000000e+03 * invT+8.486600300000000e-01 -3.587600400000000e+00 * logT+3.012735500000000e-04 * T-1.497218761666667e-07 * T2-6.291014016666667e-11 * T3+3.537338020000000e-14 * T4
;
// species 7: HO2
species[7] =
+3.014058970000000e+02 * invT-1.551889300000000e+00 -3.877615220000000e+00 * logT+5.799244400000000e-04 * T-1.775876083333333e-06 * T2+9.734864750000000e-10 * T3-2.017821605000000e-13 * T4
;
// species 8: H2O2
species[8] =
-1.766382310000000e+04 * invT-1.411878630000000e+00 -3.828191160000000e+00 * logT-1.636540200000000e-03 * T-9.304980266666667e-07 * T2+6.848160891666667e-10 * T3-1.527466825000000e-13 * T4
;
// species 9: OHV
species[9] =
+5.002130060000000e+04 * invT+2.278371300000000e+00 -3.637259270000000e+00 * logT-9.257252050000000e-05 * T+2.793854533333333e-07 * T2-1.989472975000000e-10 * T3+4.216041795000000e-14 * T4
;
// species 10: AR
species[10] =
-7.453750000000000e+02 * invT-1.879674910000000e+00 -2.500000000000000e+00 * logT-2.740033625000000e-16 * T+2.665623750000000e-19 * T2-1.552516841666667e-22 * T3+3.672898590000000e-26 * T4
;
// species 11: HE
species[11] =
-7.453750000000000e+02 * invT+1.571276026000000e+00 -2.500000000000000e+00 * logT-2.740033625000000e-16 * T+2.665623750000000e-19 * T2-1.552516841666667e-22 * T3+3.672898590000000e-26 * T4
;
}
else {
// species 0: H
species[0] =
+2.547366000000000e+04 * invT+2.946682850000000e+00 -2.500000000000000e+00 * logT
;
// species 1: H2
species[1] =
-8.167484240000000e+02 * invT+3.961172590000000e+00 -2.932865750000000e+00 * logT-4.133040130000000e-04 * T+2.440039400000000e-08 * T2-1.284170116666667e-12 * T3+3.444024000000000e-17 * T4
;
// species 2: O
species[2] =
+2.922661550000000e+04 * invT-2.379309400000000e+00 -2.543636970000000e+00 * logT+1.365812430000000e-05 * T+6.983825333333333e-10 * T2-4.129015375000000e-13 * T3+2.397768470000000e-17 * T4
;
// species 3: OH
species[3] =
+3.697880840000000e+03 * invT-3.006494770000000e+00 -2.838530330000000e+00 * logT-5.537064450000000e-04 * T+4.900003483333333e-08 * T2-3.505822741666666e-12 * T3+1.211449450000000e-16 * T4
;
// species 4: H2O
species[4] =
-2.988629310000000e+04 * invT-4.205080219999999e+00 -2.677038900000000e+00 * logT-1.486590800000000e-03 * T+1.289614816666667e-07 * T2-7.869459500000001e-12 * T3+2.134499550000000e-16 * T4
;
// species 5: O2
species[5] =
-1.214730780000000e+03 * invT+2.442514700000000e-01 -3.660960650000000e+00 * logT-3.281829055000000e-04 * T+2.352493783333333e-08 * T2-1.714982791666666e-12 * T3+6.495671800000001e-17 * T4
;
// species 6: N2
species[6] =
-9.245018400000000e+02 * invT-2.918713730000000e+00 -2.952576370000000e+00 * logT-6.984502000000001e-04 * T+8.210526716666667e-08 * T2-6.550084958333333e-12 * T3+2.303776020000000e-16 * T4
;
// species 7: HO2
species[7] =
+3.516660660000000e+01 * invT+1.210132220000000e+00 -4.172287410000000e+00 * logT-9.405881350000000e-04 * T+5.771288100000000e-08 * T2-1.622146241666667e-12 * T3-8.812845250000000e-18 * T4
;
// species 8: H2O2
species[8] =
-1.800241800000000e+04 * invT+3.909661134000000e+00 -4.579773050000000e+00 * logT-2.026630015000000e-03 * T+2.164078833333333e-07 * T2-1.651761666666667e-11 * T3+5.698439600000000e-16 * T4
;
// species 9: OHV
species[9] =
+5.030140630000000e+04 * invT-2.712986050000000e+00 -2.882730000000000e+00 * logT-5.069871500000000e-04 * T+3.794795000000000e-08 * T2-1.812235833333334e-12 * T3+2.563152500000000e-17 * T4
;
// species 10: AR
species[10] =
-7.453750000000000e+02 * invT-1.879674910000000e+00 -2.500000000000000e+00 * logT
;
// species 11: HE
species[11] =
-7.453750000000000e+02 * invT+1.571276026000000e+00 -2.500000000000000e+00 * logT
;
}
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+2.54736600e+04 * invT+1.94668285e+00 -2.50000000e+00 * logT-2.74003363e-16 * T+2.66562375e-19 * T2-1.55251684e-22 * T3+3.67289859e-26 * T4
;
// species 1: H2
species[1] =
-9.51146702e+02 * invT+2.55965697e+00 -2.72113605e+00 * logT-2.39607123e-03 * T+1.69137072e-06 * T2-7.47447457e-10 * T3+1.35331551e-13 * T4
;
// species 2: O
species[2] =
+2.91277009e+04 * invT-1.94708190e-01 -3.10652772e+00 * logT+1.37845198e-03 * T-8.52393020e-07 * T2+3.57911823e-10 * T3-6.73782845e-14 * T4
;
// species 3: OH
species[3] =
+3.36955446e+03 * invT+3.05848071e+00 -3.98454038e+00 * logT+1.16903977e-03 * T-7.38720988e-07 * T2+3.04845444e-10 * T3-6.35474270e-14 * T4
;
// species 4: H2O
species[4] =
-3.02973253e+04 * invT+4.25337355e+00 -4.23947087e+00 * logT+1.19096867e-03 * T-1.25524296e-06 * T2+5.58365908e-10 * T3-1.13900785e-13 * T4
;
// species 5: O2
species[5] =
-1.05270367e+03 * invT-1.51767977e+00 -3.65493298e+00 * logT+9.58840170e-04 * T-1.11495715e-06 * T2+4.91246345e-10 * T3-8.31709040e-14 * T4
;
// species 6: N2
species[6] =
-1.05196456e+03 * invT-1.51339970e-01 -3.58760040e+00 * logT+3.01273550e-04 * T-1.49721876e-07 * T2-6.29101402e-11 * T3+3.53733802e-14 * T4
;
// species 7: HO2
species[7] =
+3.01405897e+02 * invT-2.55188930e+00 -3.87761522e+00 * logT+5.79924440e-04 * T-1.77587608e-06 * T2+9.73486475e-10 * T3-2.01782160e-13 * T4
;
// species 8: H2O2
species[8] =
-1.76638231e+04 * invT-2.41187863e+00 -3.82819116e+00 * logT-1.63654020e-03 * T-9.30498027e-07 * T2+6.84816089e-10 * T3-1.52746682e-13 * T4
;
// species 9: OHV
species[9] =
+5.00213006e+04 * invT+1.27837130e+00 -3.63725927e+00 * logT-9.25725205e-05 * T+2.79385453e-07 * T2-1.98947298e-10 * T3+4.21604179e-14 * T4
;
// species 10: AR
species[10] =
-7.45375000e+02 * invT-2.87967491e+00 -2.50000000e+00 * logT-2.74003363e-16 * T+2.66562375e-19 * T2-1.55251684e-22 * T3+3.67289859e-26 * T4
;
// species 11: HE
species[11] =
-7.45375000e+02 * invT+5.71276026e-01 -2.50000000e+00 * logT-2.74003363e-16 * T+2.66562375e-19 * T2-1.55251684e-22 * T3+3.67289859e-26 * T4
;
}
else {
// species 0: H
species[0] =
+2.54736600e+04 * invT+1.94668285e+00 -2.50000000e+00 * logT
;
// species 1: H2
species[1] =
-8.16748424e+02 * invT+2.96117259e+00 -2.93286575e+00 * logT-4.13304013e-04 * T+2.44003940e-08 * T2-1.28417012e-12 * T3+3.44402400e-17 * T4
;
// species 2: O
species[2] =
+2.92266155e+04 * invT-3.37930940e+00 -2.54363697e+00 * logT+1.36581243e-05 * T+6.98382533e-10 * T2-4.12901538e-13 * T3+2.39776847e-17 * T4
;
// species 3: OH
species[3] =
+3.69788084e+03 * invT-4.00649477e+00 -2.83853033e+00 * logT-5.53706445e-04 * T+4.90000348e-08 * T2-3.50582274e-12 * T3+1.21144945e-16 * T4
;
// species 4: H2O
species[4] =
-2.98862931e+04 * invT-5.20508022e+00 -2.67703890e+00 * logT-1.48659080e-03 * T+1.28961482e-07 * T2-7.86945950e-12 * T3+2.13449955e-16 * T4
;
// species 5: O2
species[5] =
-1.21473078e+03 * invT-7.55748530e-01 -3.66096065e+00 * logT-3.28182906e-04 * T+2.35249378e-08 * T2-1.71498279e-12 * T3+6.49567180e-17 * T4
;
// species 6: N2
species[6] =
-9.24501840e+02 * invT-3.91871373e+00 -2.95257637e+00 * logT-6.98450200e-04 * T+8.21052672e-08 * T2-6.55008496e-12 * T3+2.30377602e-16 * T4
;
// species 7: HO2
species[7] =
+3.51666066e+01 * invT+2.10132220e-01 -4.17228741e+00 * logT-9.40588135e-04 * T+5.77128810e-08 * T2-1.62214624e-12 * T3-8.81284525e-18 * T4
;
// species 8: H2O2
species[8] =
-1.80024180e+04 * invT+2.90966113e+00 -4.57977305e+00 * logT-2.02663002e-03 * T+2.16407883e-07 * T2-1.65176167e-11 * T3+5.69843960e-16 * T4
;
// species 9: OHV
species[9] =
+5.03014063e+04 * invT-3.71298605e+00 -2.88273000e+00 * logT-5.06987150e-04 * T+3.79479500e-08 * T2-1.81223583e-12 * T3+2.56315250e-17 * T4
;
// species 10: AR
species[10] =
-7.45375000e+02 * invT-2.87967491e+00 -2.50000000e+00 * logT
;
// species 11: HE
species[11] =
-7.45375000e+02 * invT+5.71276026e-01 -2.50000000e+00 * logT
;
}
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4+2.54736600e+04 * invT
;
// species 1: H2
species[1] =
+1.72113605e+00 +2.39607123e-03 * T-3.38274143e-06 * T2+2.24234237e-09 * T3-5.41326204e-13 * T4-9.51146702e+02 * invT
;
// species 2: O
species[2] =
+2.10652772e+00 -1.37845198e-03 * T+1.70478604e-06 * T2-1.07373547e-09 * T3+2.69513138e-13 * T4+2.91277009e+04 * invT
;
// species 3: OH
species[3] =
+2.98454038e+00 -1.16903977e-03 * T+1.47744198e-06 * T2-9.14536332e-10 * T3+2.54189708e-13 * T4+3.36955446e+03 * invT
;
// species 4: H2O
species[4] =
+3.23947087e+00 -1.19096867e-03 * T+2.51048593e-06 * T2-1.67509772e-09 * T3+4.55603140e-13 * T4-3.02973253e+04 * invT
;
// species 5: O2
species[5] =
+2.65493298e+00 -9.58840170e-04 * T+2.22991431e-06 * T2-1.47373904e-09 * T3+3.32683616e-13 * T4-1.05270367e+03 * invT
;
// species 6: N2
species[6] =
+2.58760040e+00 -3.01273550e-04 * T+2.99443752e-07 * T2+1.88730421e-10 * T3-1.41493521e-13 * T4-1.05196456e+03 * invT
;
// species 7: HO2
species[7] =
+2.87761522e+00 -5.79924440e-04 * T+3.55175217e-06 * T2-2.92045942e-09 * T3+8.07128642e-13 * T4+3.01405897e+02 * invT
;
// species 8: H2O2
species[8] =
+2.82819116e+00 +1.63654020e-03 * T+1.86099605e-06 * T2-2.05444827e-09 * T3+6.10986730e-13 * T4-1.76638231e+04 * invT
;
// species 9: OHV
species[9] =
+2.63725927e+00 +9.25725205e-05 * T-5.58770907e-07 * T2+5.96841892e-10 * T3-1.68641672e-13 * T4+5.00213006e+04 * invT
;
// species 10: AR
species[10] =
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
;
// species 11: HE
species[11] =
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
;
}
else {
// species 0: H
species[0] =
+1.50000000e+00 +2.54736600e+04 * invT
;
// species 1: H2
species[1] =
+1.93286575e+00 +4.13304013e-04 * T-4.88007880e-08 * T2+3.85251035e-12 * T3-1.37760960e-16 * T4-8.16748424e+02 * invT
;
// species 2: O
species[2] =
+1.54363697e+00 -1.36581243e-05 * T-1.39676507e-09 * T2+1.23870461e-12 * T3-9.59107388e-17 * T4+2.92266155e+04 * invT
;
// species 3: OH
species[3] =
+1.83853033e+00 +5.53706445e-04 * T-9.80000697e-08 * T2+1.05174682e-11 * T3-4.84579780e-16 * T4+3.69788084e+03 * invT
;
// species 4: H2O
species[4] =
+1.67703890e+00 +1.48659080e-03 * T-2.57922963e-07 * T2+2.36083785e-11 * T3-8.53799820e-16 * T4-2.98862931e+04 * invT
;
// species 5: O2
species[5] =
+2.66096065e+00 +3.28182906e-04 * T-4.70498757e-08 * T2+5.14494837e-12 * T3-2.59826872e-16 * T4-1.21473078e+03 * invT
;
// species 6: N2
species[6] =
+1.95257637e+00 +6.98450200e-04 * T-1.64210534e-07 * T2+1.96502549e-11 * T3-9.21510408e-16 * T4-9.24501840e+02 * invT
;
// species 7: HO2
species[7] =
+3.17228741e+00 +9.40588135e-04 * T-1.15425762e-07 * T2+4.86643873e-12 * T3+3.52513810e-17 * T4+3.51666066e+01 * invT
;
// species 8: H2O2
species[8] =
+3.57977305e+00 +2.02663002e-03 * T-4.32815767e-07 * T2+4.95528500e-11 * T3-2.27937584e-15 * T4-1.80024180e+04 * invT
;
// species 9: OHV
species[9] =
+1.88273000e+00 +5.06987150e-04 * T-7.58959000e-08 * T2+5.43670750e-12 * T3-1.02526100e-16 * T4+5.03014063e+04 * invT
;
// species 10: AR
species[10] =
+1.50000000e+00 -7.45375000e+02 * invT
;
// species 11: HE
species[11] =
+1.50000000e+00 -7.45375000e+02 * invT
;
}
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4+2.54736600e+04 * invT
;
// species 1: H2
species[1] =
+2.72113605e+00 +2.39607123e-03 * T-3.38274143e-06 * T2+2.24234237e-09 * T3-5.41326204e-13 * T4-9.51146702e+02 * invT
;
// species 2: O
species[2] =
+3.10652772e+00 -1.37845198e-03 * T+1.70478604e-06 * T2-1.07373547e-09 * T3+2.69513138e-13 * T4+2.91277009e+04 * invT
;
// species 3: OH
species[3] =
+3.98454038e+00 -1.16903977e-03 * T+1.47744198e-06 * T2-9.14536332e-10 * T3+2.54189708e-13 * T4+3.36955446e+03 * invT
;
// species 4: H2O
species[4] =
+4.23947087e+00 -1.19096867e-03 * T+2.51048593e-06 * T2-1.67509772e-09 * T3+4.55603140e-13 * T4-3.02973253e+04 * invT
;
// species 5: O2
species[5] =
+3.65493298e+00 -9.58840170e-04 * T+2.22991431e-06 * T2-1.47373904e-09 * T3+3.32683616e-13 * T4-1.05270367e+03 * invT
;
// species 6: N2
species[6] =
+3.58760040e+00 -3.01273550e-04 * T+2.99443752e-07 * T2+1.88730421e-10 * T3-1.41493521e-13 * T4-1.05196456e+03 * invT
;
// species 7: HO2
species[7] =
+3.87761522e+00 -5.79924440e-04 * T+3.55175217e-06 * T2-2.92045942e-09 * T3+8.07128642e-13 * T4+3.01405897e+02 * invT
;
// species 8: H2O2
species[8] =
+3.82819116e+00 +1.63654020e-03 * T+1.86099605e-06 * T2-2.05444827e-09 * T3+6.10986730e-13 * T4-1.76638231e+04 * invT
;
// species 9: OHV
species[9] =
+3.63725927e+00 +9.25725205e-05 * T-5.58770907e-07 * T2+5.96841892e-10 * T3-1.68641672e-13 * T4+5.00213006e+04 * invT
;
// species 10: AR
species[10] =
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
;
// species 11: HE
species[11] =
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
;
}
else {
// species 0: H
species[0] =
+2.50000000e+00 +2.54736600e+04 * invT
;
// species 1: H2
species[1] =
+2.93286575e+00 +4.13304013e-04 * T-4.88007880e-08 * T2+3.85251035e-12 * T3-1.37760960e-16 * T4-8.16748424e+02 * invT
;
// species 2: O
species[2] =
+2.54363697e+00 -1.36581243e-05 * T-1.39676507e-09 * T2+1.23870461e-12 * T3-9.59107388e-17 * T4+2.92266155e+04 * invT
;
// species 3: OH
species[3] =
+2.83853033e+00 +5.53706445e-04 * T-9.80000697e-08 * T2+1.05174682e-11 * T3-4.84579780e-16 * T4+3.69788084e+03 * invT
;
// species 4: H2O
species[4] =
+2.67703890e+00 +1.48659080e-03 * T-2.57922963e-07 * T2+2.36083785e-11 * T3-8.53799820e-16 * T4-2.98862931e+04 * invT
;
// species 5: O2
species[5] =
+3.66096065e+00 +3.28182906e-04 * T-4.70498757e-08 * T2+5.14494837e-12 * T3-2.59826872e-16 * T4-1.21473078e+03 * invT
;
// species 6: N2
species[6] =
+2.95257637e+00 +6.98450200e-04 * T-1.64210534e-07 * T2+1.96502549e-11 * T3-9.21510408e-16 * T4-9.24501840e+02 * invT
;
// species 7: HO2
species[7] =
+4.17228741e+00 +9.40588135e-04 * T-1.15425762e-07 * T2+4.86643873e-12 * T3+3.52513810e-17 * T4+3.51666066e+01 * invT
;
// species 8: H2O2
species[8] =
+4.57977305e+00 +2.02663002e-03 * T-4.32815767e-07 * T2+4.95528500e-11 * T3-2.27937584e-15 * T4-1.80024180e+04 * invT
;
// species 9: OHV
species[9] =
+2.88273000e+00 +5.06987150e-04 * T-7.58959000e-08 * T2+5.43670750e-12 * T3-1.02526100e-16 * T4+5.03014063e+04 * invT
;
// species 10: AR
species[10] =
+2.50000000e+00 -7.45375000e+02 * invT
;
// species 11: HE
species[11] =
+2.50000000e+00 -7.45375000e+02 * invT
;
}
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real logT = log(T);


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+2.50000000e+00 * logT+5.48006725e-16 * T-7.99687125e-19 * T2+6.21006737e-22 * T3-1.83644929e-25 * T4-4.46682850e-01 
;
// species 1: H2
species[1] =
+2.72113605e+00 * logT+4.79214245e-03 * T-5.07411215e-06 * T2+2.98978983e-09 * T3-6.76657755e-13 * T4-8.38520919e-01 
;
// species 2: O
species[2] =
+3.10652772e+00 * logT-2.75690396e-03 * T+2.55717906e-06 * T2-1.43164729e-09 * T3+3.36891422e-13 * T4+2.30123591e+00 
;
// species 3: OH
species[3] =
+3.98454038e+00 * logT-2.33807954e-03 * T+2.21616296e-06 * T2-1.21938178e-09 * T3+3.17737135e-13 * T4-7.39403338e-02 
;
// species 4: H2O
species[4] =
+4.23947087e+00 * logT-2.38193734e-03 * T+3.76572889e-06 * T2-2.23346363e-09 * T3+5.69503925e-13 * T4-1.01390268e+00 
;
// species 5: O2
species[5] =
+3.65493298e+00 * logT-1.91768034e-03 * T+3.34487146e-06 * T2-1.96498538e-09 * T3+4.15854520e-13 * T4+4.17261275e+00 
;
// species 6: N2
species[6] =
+3.58760040e+00 * logT-6.02547100e-04 * T+4.49165629e-07 * T2+2.51640561e-10 * T3-1.76866901e-13 * T4+2.73894037e+00 
;
// species 7: HO2
species[7] =
+3.87761522e+00 * logT-1.15984888e-03 * T+5.32762825e-06 * T2-3.89394590e-09 * T3+1.00891080e-12 * T4+5.42950452e+00 
;
// species 8: H2O2
species[8] =
+3.82819116e+00 * logT+3.27308040e-03 * T+2.79149408e-06 * T2-2.73926436e-09 * T3+7.63733412e-13 * T4+5.24006979e+00 
;
// species 9: OHV
species[9] =
+3.63725927e+00 * logT+1.85145041e-04 * T-8.38156360e-07 * T2+7.95789190e-10 * T3-2.10802090e-13 * T4+1.35888797e+00 
;
// species 10: AR
species[10] =
+2.50000000e+00 * logT+5.48006725e-16 * T-7.99687125e-19 * T2+6.21006737e-22 * T3-1.83644929e-25 * T4+4.37967491e+00 
;
// species 11: HE
species[11] =
+2.50000000e+00 * logT+5.48006725e-16 * T-7.99687125e-19 * T2+6.21006737e-22 * T3-1.83644929e-25 * T4+9.28723974e-01 
;
}
else {
// species 0: H
species[0] =
+2.50000000e+00 * logT-4.46682850e-01 
;
// species 1: H2
species[1] =
+2.93286575e+00 * logT+8.26608026e-04 * T-7.32011820e-08 * T2+5.13668047e-12 * T3-1.72201200e-16 * T4-1.02830684e+00 
;
// species 2: O
species[2] =
+2.54363697e+00 * logT-2.73162486e-05 * T-2.09514760e-09 * T2+1.65160615e-12 * T3-1.19888423e-16 * T4+4.92294637e+00 
;
// species 3: OH
species[3] =
+2.83853033e+00 * logT+1.10741289e-03 * T-1.47000104e-07 * T2+1.40232910e-11 * T3-6.05724725e-16 * T4+5.84502510e+00 
;
// species 4: H2O
species[4] =
+2.67703890e+00 * logT+2.97318160e-03 * T-3.86884445e-07 * T2+3.14778380e-11 * T3-1.06724977e-15 * T4+6.88211912e+00 
;
// species 5: O2
species[5] =
+3.66096065e+00 * logT+6.56365811e-04 * T-7.05748135e-08 * T2+6.85993117e-12 * T3-3.24783590e-16 * T4+3.41670918e+00 
;
// species 6: N2
species[6] =
+2.95257637e+00 * logT+1.39690040e-03 * T-2.46315802e-07 * T2+2.62003398e-11 * T3-1.15188801e-15 * T4+5.87129010e+00 
;
// species 7: HO2
species[7] =
+4.17228741e+00 * logT+1.88117627e-03 * T-1.73138643e-07 * T2+6.48858497e-12 * T3+4.40642262e-17 * T4+2.96215519e+00 
;
// species 8: H2O2
species[8] =
+4.57977305e+00 * logT+4.05326003e-03 * T-6.49223650e-07 * T2+6.60704667e-11 * T3-2.84921980e-15 * T4+6.70111916e-01 
;
// species 9: OHV
species[9] =
+2.88273000e+00 * logT+1.01397430e-03 * T-1.13843850e-07 * T2+7.24894333e-12 * T3-1.28157625e-16 * T4+5.59571605e+00 
;
// species 10: AR
species[10] =
+2.50000000e+00 * logT+4.37967491e+00 
;
// species 11: HE
species[11] =
+2.50000000e+00 * logT+9.28723974e-01 
;
}
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, const amrex::Real T)
{
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
species[0] =
+5.48006725e-16 -3.19874850e-18 * T+5.58906063e-21 * T2-2.93831887e-24 * T3
;
// species 1: H2
species[1] =
+4.79214245e-03 -2.02964486e-05 * T+2.69081085e-08 * T2-1.08265241e-11 * T3
;
// species 2: O
species[2] =
-2.75690396e-03 +1.02287162e-05 * T-1.28848256e-08 * T2+5.39026276e-12 * T3
;
// species 3: OH
species[3] =
-2.33807954e-03 +8.86465186e-06 * T-1.09744360e-08 * T2+5.08379416e-12 * T3
;
// species 4: H2O
species[4] =
-2.38193734e-03 +1.50629156e-05 * T-2.01011727e-08 * T2+9.11206280e-12 * T3
;
// species 5: O2
species[5] =
-1.91768034e-03 +1.33794858e-05 * T-1.76848684e-08 * T2+6.65367232e-12 * T3
;
// species 6: N2
species[6] =
-6.02547100e-04 +1.79666251e-06 * T+2.26476505e-09 * T2-2.82987042e-12 * T3
;
// species 7: HO2
species[7] =
-1.15984888e-03 +2.13105130e-05 * T-3.50455131e-08 * T2+1.61425728e-11 * T3
;
// species 8: H2O2
species[8] =
+3.27308040e-03 +1.11659763e-05 * T-2.46533792e-08 * T2+1.22197346e-11 * T3
;
// species 9: OHV
species[9] =
+1.85145041e-04 -3.35262544e-06 * T+7.16210271e-09 * T2-3.37283344e-12 * T3
;
// species 10: AR
species[10] =
+5.48006725e-16 -3.19874850e-18 * T+5.58906063e-21 * T2-2.93831887e-24 * T3
;
// species 11: HE
species[11] =
+5.48006725e-16 -3.19874850e-18 * T+5.58906063e-21 * T2-2.93831887e-24 * T3
;
}
else {
// species 0: H
species[0] =
0.0
;
// species 1: H2
species[1] =
+8.26608026e-04 -2.92804728e-07 * T+4.62301242e-11 * T2-2.75521920e-15 * T3
;
// species 2: O
species[2] =
-2.73162486e-05 -8.38059040e-09 * T+1.48644553e-11 * T2-1.91821478e-15 * T3
;
// species 3: OH
species[3] =
+1.10741289e-03 -5.88000418e-07 * T+1.26209619e-10 * T2-9.69159560e-15 * T3
;
// species 4: H2O
species[4] =
+2.97318160e-03 -1.54753778e-06 * T+2.83300542e-10 * T2-1.70759964e-14 * T3
;
// species 5: O2
species[5] =
+6.56365811e-04 -2.82299254e-07 * T+6.17393805e-11 * T2-5.19653744e-15 * T3
;
// species 6: N2
species[6] =
+1.39690040e-03 -9.85263206e-07 * T+2.35803059e-10 * T2-1.84302082e-14 * T3
;
// species 7: HO2
species[7] =
+1.88117627e-03 -6.92554572e-07 * T+5.83972647e-11 * T2+7.05027620e-16 * T3
;
// species 8: H2O2
species[8] =
+4.05326003e-03 -2.59689460e-06 * T+5.94634200e-10 * T2-4.55875168e-14 * T3
;
// species 9: OHV
species[9] =
+1.01397430e-03 -4.55375400e-07 * T+6.52404900e-11 * T2-2.05052200e-15 * T3
;
// species 10: AR
species[10] =
0.0
;
// species 11: HE
species[11] =
0.0
;
}
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
amrex::Real result = 0; 
amrex::Real cpor[12]; //  temporary storage
cp_R(cpor, T);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*cpor[id];
}

cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
amrex::Real result = 0.0; 

// compute Cp/R at the given temperature
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
result += y[0] * (
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+2.72113605e+00 +4.79214245e-03 * T-1.01482243e-05 * T2+8.96936949e-09 * T3-2.70663102e-12 * T4
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+3.10652772e+00 -2.75690396e-03 * T+5.11435812e-06 * T2-4.29494187e-09 * T3+1.34756569e-12 * T4
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+3.98454038e+00 -2.33807954e-03 * T+4.43232593e-06 * T2-3.65814533e-09 * T3+1.27094854e-12 * T4
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+4.23947087e+00 -2.38193734e-03 * T+7.53145778e-06 * T2-6.70039089e-09 * T3+2.27801570e-12 * T4
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+3.65493298e+00 -1.91768034e-03 * T+6.68974292e-06 * T2-5.89495614e-09 * T3+1.66341808e-12 * T4
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+3.58760040e+00 -6.02547100e-04 * T+8.98331257e-07 * T2+7.54921682e-10 * T3-7.07467604e-13 * T4
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+3.87761522e+00 -1.15984888e-03 * T+1.06552565e-05 * T2-1.16818377e-08 * T3+4.03564321e-12 * T4
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+3.82819116e+00 +3.27308040e-03 * T+5.58298816e-06 * T2-8.21779307e-09 * T3+3.05493365e-12 * T4
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+3.63725927e+00 +1.85145041e-04 * T-1.67631272e-06 * T2+2.38736757e-09 * T3-8.43208359e-13 * T4
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+2.50000000e+00 +5.48006725e-16 * T-1.59937425e-18 * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.2498374807187924
;
}
else {
// species 0: H
result += y[0] * (
+2.50000000e+00 
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+2.93286575e+00 +8.26608026e-04 * T-1.46402364e-07 * T2+1.54100414e-11 * T3-6.88804800e-16 * T4
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+2.54363697e+00 -2.73162486e-05 * T-4.19029520e-09 * T2+4.95481845e-12 * T3-4.79553694e-16 * T4
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+2.83853033e+00 +1.10741289e-03 * T-2.94000209e-07 * T2+4.20698729e-11 * T3-2.42289890e-15 * T4
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+2.67703890e+00 +2.97318160e-03 * T-7.73768890e-07 * T2+9.44335140e-11 * T3-4.26899910e-15 * T4
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+3.66096065e+00 +6.56365811e-04 * T-1.41149627e-07 * T2+2.05797935e-11 * T3-1.29913436e-15 * T4
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+2.95257637e+00 +1.39690040e-03 * T-4.92631603e-07 * T2+7.86010195e-11 * T3-4.60755204e-15 * T4
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+4.17228741e+00 +1.88117627e-03 * T-3.46277286e-07 * T2+1.94657549e-11 * T3+1.76256905e-16 * T4
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+4.57977305e+00 +4.05326003e-03 * T-1.29844730e-06 * T2+1.98211400e-10 * T3-1.13968792e-14 * T4
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+2.88273000e+00 +1.01397430e-03 * T-2.27687700e-07 * T2+2.17468300e-11 * T3-5.12630500e-16 * T4
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+2.50000000e+00 
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+2.50000000e+00 
)* 0.2498374807187924
;
}

cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
amrex::Real result = 0; 
amrex::Real cvor[12]; //  temporary storage
cv_R(cvor, T);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*cvor[id];
}

cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS(const amrex::Real T, const amrex::Real y[],  amrex::Real& cvbs)
{
amrex::Real result = 0.0; 
// compute Cv/R at the given temperature
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
result += y[0] * (
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+1.72113605e+00 +4.79214245e-03 * T-1.01482243e-05  * T2+8.96936949e-09 * T3-2.70663102e-12 * T4
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+2.10652772e+00 -2.75690396e-03 * T+5.11435812e-06  * T2-4.29494187e-09 * T3+1.34756569e-12 * T4
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+2.98454038e+00 -2.33807954e-03 * T+4.43232593e-06  * T2-3.65814533e-09 * T3+1.27094854e-12 * T4
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+3.23947087e+00 -2.38193734e-03 * T+7.53145778e-06  * T2-6.70039089e-09 * T3+2.27801570e-12 * T4
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+2.65493298e+00 -1.91768034e-03 * T+6.68974292e-06  * T2-5.89495614e-09 * T3+1.66341808e-12 * T4
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+2.58760040e+00 -6.02547100e-04 * T+8.98331257e-07  * T2+7.54921682e-10 * T3-7.07467604e-13 * T4
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+2.87761522e+00 -1.15984888e-03 * T+1.06552565e-05  * T2-1.16818377e-08 * T3+4.03564321e-12 * T4
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+2.82819116e+00 +3.27308040e-03 * T+5.58298816e-06  * T2-8.21779307e-09 * T3+3.05493365e-12 * T4
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+2.63725927e+00 +1.85145041e-04 * T-1.67631272e-06  * T2+2.38736757e-09 * T3-8.43208359e-13 * T4
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+1.50000000e+00 +5.48006725e-16 * T-1.59937425e-18  * T2+1.86302021e-21 * T3-7.34579718e-25 * T4
)* 0.2498374807187924
;
}
else {
// species 0: H
result += y[0] * (
+1.50000000e+00 
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+1.93286575e+00 +8.26608026e-04 * T-1.46402364e-07  * T2+1.54100414e-11 * T3-6.88804800e-16 * T4
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+1.54363697e+00 -2.73162486e-05 * T-4.19029520e-09  * T2+4.95481845e-12 * T3-4.79553694e-16 * T4
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+1.83853033e+00 +1.10741289e-03 * T-2.94000209e-07  * T2+4.20698729e-11 * T3-2.42289890e-15 * T4
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+1.67703890e+00 +2.97318160e-03 * T-7.73768890e-07  * T2+9.44335140e-11 * T3-4.26899910e-15 * T4
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+2.66096065e+00 +6.56365811e-04 * T-1.41149627e-07  * T2+2.05797935e-11 * T3-1.29913436e-15 * T4
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+1.95257637e+00 +1.39690040e-03 * T-4.92631603e-07  * T2+7.86010195e-11 * T3-4.60755204e-15 * T4
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+3.17228741e+00 +1.88117627e-03 * T-3.46277286e-07  * T2+1.94657549e-11 * T3+1.76256905e-16 * T4
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+3.57977305e+00 +4.05326003e-03 * T-1.29844730e-06  * T2+1.98211400e-10 * T3-1.13968792e-14 * T4
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+1.88273000e+00 +1.01397430e-03 * T-2.27687700e-07  * T2+2.17468300e-11 * T3-5.12630500e-16 * T4
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+1.50000000e+00 
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+1.50000000e+00 
)* 0.2498374807187924
;
}

cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
amrex::Real result = 0; 
amrex::Real hml[12]; //  temporary storage
amrex::Real RT = 8.31446261815324e+07*T; // R*T
speciesEnthalpy(hml, T);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*hml[id];
}

hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS(const amrex::Real T, const amrex::Real y[],  amrex::Real& hbms)
{
amrex::Real result = 0.0;
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
result += y[0] * (
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4+2.54736600e+04 * invT
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+2.72113605e+00 +2.39607123e-03 * T-3.38274143e-06 * T2+2.24234237e-09 * T3-5.41326204e-13 * T4-9.51146702e+02 * invT
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+3.10652772e+00 -1.37845198e-03 * T+1.70478604e-06 * T2-1.07373547e-09 * T3+2.69513138e-13 * T4+2.91277009e+04 * invT
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+3.98454038e+00 -1.16903977e-03 * T+1.47744198e-06 * T2-9.14536332e-10 * T3+2.54189708e-13 * T4+3.36955446e+03 * invT
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+4.23947087e+00 -1.19096867e-03 * T+2.51048593e-06 * T2-1.67509772e-09 * T3+4.55603140e-13 * T4-3.02973253e+04 * invT
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+3.65493298e+00 -9.58840170e-04 * T+2.22991431e-06 * T2-1.47373904e-09 * T3+3.32683616e-13 * T4-1.05270367e+03 * invT
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+3.58760040e+00 -3.01273550e-04 * T+2.99443752e-07 * T2+1.88730421e-10 * T3-1.41493521e-13 * T4-1.05196456e+03 * invT
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+3.87761522e+00 -5.79924440e-04 * T+3.55175217e-06 * T2-2.92045942e-09 * T3+8.07128642e-13 * T4+3.01405897e+02 * invT
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+3.82819116e+00 +1.63654020e-03 * T+1.86099605e-06 * T2-2.05444827e-09 * T3+6.10986730e-13 * T4-1.76638231e+04 * invT
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+3.63725927e+00 +9.25725205e-05 * T-5.58770907e-07 * T2+5.96841892e-10 * T3-1.68641672e-13 * T4+5.00213006e+04 * invT
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+2.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
)* 0.2498374807187924
;
}
else {
// species 0: H
result += y[0] * (
+2.50000000e+00 +2.54736600e+04 * invT
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+2.93286575e+00 +4.13304013e-04 * T-4.88007880e-08 * T2+3.85251035e-12 * T3-1.37760960e-16 * T4-8.16748424e+02 * invT
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+2.54363697e+00 -1.36581243e-05 * T-1.39676507e-09 * T2+1.23870461e-12 * T3-9.59107388e-17 * T4+2.92266155e+04 * invT
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+2.83853033e+00 +5.53706445e-04 * T-9.80000697e-08 * T2+1.05174682e-11 * T3-4.84579780e-16 * T4+3.69788084e+03 * invT
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+2.67703890e+00 +1.48659080e-03 * T-2.57922963e-07 * T2+2.36083785e-11 * T3-8.53799820e-16 * T4-2.98862931e+04 * invT
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+3.66096065e+00 +3.28182906e-04 * T-4.70498757e-08 * T2+5.14494837e-12 * T3-2.59826872e-16 * T4-1.21473078e+03 * invT
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+2.95257637e+00 +6.98450200e-04 * T-1.64210534e-07 * T2+1.96502549e-11 * T3-9.21510408e-16 * T4-9.24501840e+02 * invT
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+4.17228741e+00 +9.40588135e-04 * T-1.15425762e-07 * T2+4.86643873e-12 * T3+3.52513810e-17 * T4+3.51666066e+01 * invT
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+4.57977305e+00 +2.02663002e-03 * T-4.32815767e-07 * T2+4.95528500e-11 * T3-2.27937584e-15 * T4-1.80024180e+04 * invT
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+2.88273000e+00 +5.06987150e-04 * T-7.58959000e-08 * T2+5.43670750e-12 * T3-1.02526100e-16 * T4+5.03014063e+04 * invT
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+2.50000000e+00 -7.45375000e+02 * invT
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+2.50000000e+00 -7.45375000e+02 * invT
)* 0.2498374807187924
;
}

const amrex::Real RT = 8.31446261815324e+07*T; // R*T

hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
amrex::Real result = 0; 
amrex::Real uml[12]; //  temporary energy array
amrex::Real RT = 8.31446261815324e+07*T; // R*T
speciesInternalEnergy(uml, T);

// perform dot product
for (int id = 0; id < 12; ++id) {
result += x[id]*uml[id];
}

ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
amrex::Real result = 0.0;
const amrex::Real T2 = T*T;
const amrex::Real T3 = T*T*T;
const amrex::Real T4 = T*T*T*T;
const amrex::Real invT = 1.0 / T;


// species with midpoint at T=1000 kelvin
if (T < 1000) {
// species 0: H
result += y[0] * (
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4+2.54736600e+04 * invT
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+1.72113605e+00 +2.39607123e-03 * T-3.38274143e-06 * T2+2.24234237e-09 * T3-5.41326204e-13 * T4-9.51146702e+02 * invT
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+2.10652772e+00 -1.37845198e-03 * T+1.70478604e-06 * T2-1.07373547e-09 * T3+2.69513138e-13 * T4+2.91277009e+04 * invT
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+2.98454038e+00 -1.16903977e-03 * T+1.47744198e-06 * T2-9.14536332e-10 * T3+2.54189708e-13 * T4+3.36955446e+03 * invT
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+3.23947087e+00 -1.19096867e-03 * T+2.51048593e-06 * T2-1.67509772e-09 * T3+4.55603140e-13 * T4-3.02973253e+04 * invT
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+2.65493298e+00 -9.58840170e-04 * T+2.22991431e-06 * T2-1.47373904e-09 * T3+3.32683616e-13 * T4-1.05270367e+03 * invT
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+2.58760040e+00 -3.01273550e-04 * T+2.99443752e-07 * T2+1.88730421e-10 * T3-1.41493521e-13 * T4-1.05196456e+03 * invT
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+2.87761522e+00 -5.79924440e-04 * T+3.55175217e-06 * T2-2.92045942e-09 * T3+8.07128642e-13 * T4+3.01405897e+02 * invT
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+2.82819116e+00 +1.63654020e-03 * T+1.86099605e-06 * T2-2.05444827e-09 * T3+6.10986730e-13 * T4-1.76638231e+04 * invT
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+2.63725927e+00 +9.25725205e-05 * T-5.58770907e-07 * T2+5.96841892e-10 * T3-1.68641672e-13 * T4+5.00213006e+04 * invT
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+1.50000000e+00 +2.74003363e-16 * T-5.33124750e-19 * T2+4.65755052e-22 * T3-1.46915944e-25 * T4-7.45375000e+02 * invT
)* 0.2498374807187924
;
}
else {
// species 0: H
result += y[0] * (
+1.50000000e+00 +2.54736600e+04 * invT
)* 0.9920634920634921
;
// species 1: H2
result += y[1] * (
+1.93286575e+00 +4.13304013e-04 * T-4.88007880e-08 * T2+3.85251035e-12 * T3-1.37760960e-16 * T4-8.16748424e+02 * invT
)* 0.4960317460317460
;
// species 2: O
result += y[2] * (
+1.54363697e+00 -1.36581243e-05 * T-1.39676507e-09 * T2+1.23870461e-12 * T3-9.59107388e-17 * T4+2.92266155e+04 * invT
)* 0.0625039064941559
;
// species 3: OH
result += y[3] * (
+1.83853033e+00 +5.53706445e-04 * T-9.80000697e-08 * T2+1.05174682e-11 * T3-4.84579780e-16 * T4+3.69788084e+03 * invT
)* 0.0587993179279120
;
// species 4: H2O
result += y[4] * (
+1.67703890e+00 +1.48659080e-03 * T-2.57922963e-07 * T2+2.36083785e-11 * T3-8.53799820e-16 * T4-2.98862931e+04 * invT
)* 0.0555092978073827
;
// species 5: O2
result += y[5] * (
+2.66096065e+00 +3.28182906e-04 * T-4.70498757e-08 * T2+5.14494837e-12 * T3-2.59826872e-16 * T4-1.21473078e+03 * invT
)* 0.0312519532470779
;
// species 6: N2
result += y[6] * (
+1.95257637e+00 +6.98450200e-04 * T-1.64210534e-07 * T2+1.96502549e-11 * T3-9.21510408e-16 * T4-9.24501840e+02 * invT
)* 0.0356964374955379
;
// species 7: HO2
result += y[7] * (
+3.17228741e+00 +9.40588135e-04 * T-1.15425762e-07 * T2+4.86643873e-12 * T3+3.52513810e-17 * T4+3.51666066e+01 * invT
)* 0.0302975216627280
;
// species 8: H2O2
result += y[8] * (
+3.57977305e+00 +2.02663002e-03 * T-4.32815767e-07 * T2+4.95528500e-11 * T3-2.27937584e-15 * T4-1.80024180e+04 * invT
)* 0.0293996589639560
;
// species 9: OHV
result += y[9] * (
+1.88273000e+00 +5.06987150e-04 * T-7.58959000e-08 * T2+5.43670750e-12 * T3-1.02526100e-16 * T4+5.03014063e+04 * invT
)* 0.0587993179279120
;
// species 10: AR
result += y[10] * (
+1.50000000e+00 -7.45375000e+02 * invT
)* 0.0250312891113892
;
// species 11: HE
result += y[11] * (
+1.50000000e+00 -7.45375000e+02 * invT
)* 0.2498374807187924
;
}

const amrex::Real RT = 8.31446261815324e+07*T; // R*T

ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real& sbml)
{
amrex::Real result = 0; 
// Log of normalized pressure in cgs units dynes/cm^2 by Patm
amrex::Real logPratio = log ( P / 1013250.0 ); 
amrex::Real sor[12]; //  temporary storage
speciesEntropy(sor, T);

// Compute Eq 42
for (int id = 0; id < 12; ++id) {
result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
}

sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real& sbms)
{
amrex::Real result = 0; 
// Log of normalized pressure in cgs units dynes/cm^2 by Patm
amrex::Real logPratio = log ( P / 1013250.0 ); 
amrex::Real sor[12]; //  temporary storage
amrex::Real x[12]; //  need a ytx conversion
amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
YOW += y[i]* imw(i);
}
// Now compute y to x conversion
x[0] = y[0]/(1.008000*YOW); 
x[1] = y[1]/(2.016000*YOW); 
x[2] = y[2]/(15.999000*YOW); 
x[3] = y[3]/(17.007000*YOW); 
x[4] = y[4]/(18.015000*YOW); 
x[5] = y[5]/(31.998000*YOW); 
x[6] = y[6]/(28.014000*YOW); 
x[7] = y[7]/(33.006000*YOW); 
x[8] = y[8]/(34.014000*YOW); 
x[9] = y[9]/(17.007000*YOW); 
x[10] = y[10]/(39.950000*YOW); 
x[11] = y[11]/(4.002602*YOW); 
speciesEntropy(sor, T);
// Perform computation in Eq 42 and 43
for (int i = 0; i < 12; i++)
{
result += x[i]*(sor[i]-log((x[i]+1e-100))-logPratio);
}
// Scale by R/W
sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(const amrex::Real e, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
const int maxiter = 5000;
const amrex::Real tol = 1.e-12;
#else
const int maxiter = 200;
const amrex::Real tol = 1.e-6;
#endif
amrex::Real tmin = 90;// max lower bound for thermo def
amrex::Real tmax = 4000;// min upper bound for thermo def
amrex::Real e1,emin,emax,cv,t1,dt;
CKUBMS(tmin, y, emin);
CKUBMS(tmax, y, emax);
if (e < emin) {
// Linear Extrapolation below tmin
CKCVBS(tmin, y, cv);
t = tmin - (emin-e)/cv;
ierr = 1;
return;
}
if (e > emax) {
// Linear Extrapolation above tmax
CKCVBS(tmax, y, cv);
t = tmax - (emax-e)/cv;
ierr = 1;
return;
}
t1 = t;
if (t1 < tmin || t1 > tmax) {
t1 = tmin + (tmax-tmin)/(emax-emin)*(e-emin);
}
for (int i = 0; i < maxiter; ++i) {
CKUBMS(t1,y,e1);
CKCVBS(t1,y,cv);
dt = (e - e1) / cv;
if (dt > 100.) { dt = 100.; }
else if (dt < -100.) { dt = -100.; }
else if (fabs(dt) < tol) {break;}
t1 += dt;
}
t = t1;
ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(const amrex::Real h, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
const int maxiter = 5000;
const amrex::Real tol = 1.e-12;
#else
const int maxiter = 200;
const amrex::Real tol = 1.e-6;
#endif
amrex::Real tmin = 90;// max lower bound for thermo def
amrex::Real tmax = 4000;// min upper bound for thermo def
amrex::Real h1,hmin,hmax,cp,t1,dt;
CKHBMS(tmin, y, hmin);
CKHBMS(tmax, y, hmax);
if (h < hmin) {
// Linear Extrapolation below tmin
CKCPBS(tmin, y, cp);
t = tmin - (hmin-h)/cp;
ierr = 1;
return;
}
if (h > hmax) {
// Linear Extrapolation above tmax
CKCPBS(tmax, y, cp);
t = tmax - (hmax-h)/cp;
ierr = 1;
return;
}
t1 = t;
if (t1 < tmin || t1 > tmax) {
t1 = tmin + (tmax-tmin)/(hmax-hmin)*(h-hmin);
}
for (int i = 0; i < maxiter; ++i) {
CKHBMS(t1,y,h1);
CKCPBS(t1,y,cp);
dt = (h - h1) / cp;
if (dt > 100.) { dt = 100.; }
else if (dt < -100.) { dt = -100.; }
else if (fabs(dt) < tol) {break;}
t1 += dt;
}
t = t1;
ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(const amrex::Real rho, const amrex::Real T, const amrex::Real x[], amrex::Real& P)
{
amrex::Real XW = 0;//  To hold mean molecular wt
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W

}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real& P)
{
amrex::Real YOW = 0;//  for computing mean MW

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}
P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W

}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(const amrex::Real rho, const amrex::Real T, const amrex::Real c[],  amrex::Real& P)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*1.008000; // H
W += c[1]*2.016000; // H2
W += c[2]*15.999000; // O
W += c[3]*17.007000; // OH
W += c[4]*18.015000; // H2O
W += c[5]*31.998000; // O2
W += c[6]*28.014000; // N2
W += c[7]*33.006000; // HO2
W += c[8]*34.014000; // H2O2
W += c[9]*17.007000; // OHV
W += c[10]*39.950000; // AR
W += c[11]*4.002602; // HE

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W

}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real& rho)
{
amrex::Real XW = 0;//  To hold mean molecular wt
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)

}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real& rho)
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

rho = P / (8.31446261815324e+07 * T * YOW);// rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(const amrex::Real P, const amrex::Real T,const amrex::Real c[], amrex::Real& rho)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*1.008000; // H
W += c[1]*2.016000; // H2
W += c[2]*15.999000; // O
W += c[3]*17.007000; // OH
W += c[4]*18.015000; // H2O
W += c[5]*31.998000; // O2
W += c[6]*28.014000; // N2
W += c[7]*33.006000; // HO2
W += c[8]*34.014000; // H2O2
W += c[9]*17.007000; // OHV
W += c[10]*39.950000; // AR
W += c[11]*4.002602; // HE

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)

}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real wt[])
{
get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(const amrex::Real y[], amrex::Real& wtm)
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(const amrex::Real x[],  amrex::Real& wtm)
{
amrex::Real XW = 0;//  see Eq 4 in CK Manual
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
wtm = XW;

}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(const amrex::Real c[],  amrex::Real& wtm)
{
// See Eq 5 in CK Manual
amrex::Real W = 0;
amrex::Real sumC = 0;
W += c[0]*1.008000; // H
W += c[1]*2.016000; // H2
W += c[2]*15.999000; // O
W += c[3]*17.007000; // OH
W += c[4]*18.015000; // H2O
W += c[5]*31.998000; // O2
W += c[6]*28.014000; // N2
W += c[7]*33.006000; // HO2
W += c[8]*34.014000; // H2O2
W += c[9]*17.007000; // OHV
W += c[10]*39.950000; // AR
W += c[11]*4.002602; // HE

for (int id = 0; id < 12; ++id) {
sumC += c[id];
}
//  CK provides no guard against division by zero
wtm = W/sumC;

}

// get Cp/R as a function of T 
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(const amrex::Real T, amrex::Real cpor[])
{
cp_R(cpor, T);
}

// get H/RT as a function of T 
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(const amrex::Real T, amrex::Real hort[])
{
speciesEnthalpy(hort, T);
}

// get S/R as a function of T 
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(const amrex::Real T, amrex::Real sor[])
{
speciesEntropy(sor, T);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(const amrex::Real y[], amrex::Real x[])
{
amrex::Real YOW = 0;

for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}

amrex::Real YOWINV = 1.0/YOW;

for (int i = 0; i < 12; i++)
{
x[i] = y[i]*imw(i)*YOWINV;
}
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real c[])
{
amrex::Real YOW = 0;
amrex::Real PWORT;

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
c[i] = y[i]*imw(i);
}
for (int i = 0; i < 12; i++)
{
YOW += c[i];
}

// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// Now compute conversion

for (int i = 0; i < 12; i++)
{
c[i] = PWORT * y[i] * imw(i);
}
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(const amrex::Real rho, amrex::Real /*T*/, const amrex::Real y[],  amrex::Real c[])
{

for (int i = 0; i < 12; i++)
{
c[i] = rho  * y[i] * imw(i);
}
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(const amrex::Real x[], amrex::Real y[])
{
amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
// Compute mean molecular wt first
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
// Now compute conversion
amrex::Real XWinv = 1.0/XW;
y[0] = x[0]*1.008000*XWinv; 
y[1] = x[1]*2.016000*XWinv; 
y[2] = x[2]*15.999000*XWinv; 
y[3] = x[3]*17.007000*XWinv; 
y[4] = x[4]*18.015000*XWinv; 
y[5] = x[5]*31.998000*XWinv; 
y[6] = x[6]*28.014000*XWinv; 
y[7] = x[7]*33.006000*XWinv; 
y[8] = x[8]*34.014000*XWinv; 
y[9] = x[9]*17.007000*XWinv; 
y[10] = x[10]*39.950000*XWinv; 
y[11] = x[11]*4.002602*XWinv; 

}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real c[])
{
amrex::Real PORT = P/(8.31446261815324e+07 * T); // P/RT

// Compute conversion, see Eq 10
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*PORT;
}

}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(const amrex::Real rho, const amrex::Real /*T*/, const amrex::Real x[], amrex::Real c[])
{
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
ROW = rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*ROW;
}

}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(const amrex::Real c[], amrex::Real x[])
{
amrex::Real sumC = 0; 

// compute sum of c 
for (int id = 0; id < 12; ++id) {
sumC += c[id];
}

//  See Eq 13 
amrex::Real sumCinv = 1.0/sumC;
for (int id = 0; id < 12; ++id) {
x[id] = c[id]*sumCinv;
}

}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(const amrex::Real c[], amrex::Real y[])
{
amrex::Real CW = 0; // See Eq 12 in CK Manual
// compute denominator in eq 12 first
CW += c[0]*1.008000; // H
CW += c[1]*2.016000; // H2
CW += c[2]*15.999000; // O
CW += c[3]*17.007000; // OH
CW += c[4]*18.015000; // H2O
CW += c[5]*31.998000; // O2
CW += c[6]*28.014000; // N2
CW += c[7]*33.006000; // HO2
CW += c[8]*34.014000; // H2O2
CW += c[9]*17.007000; // OHV
CW += c[10]*39.950000; // AR
CW += c[11]*4.002602; // HE
// Now compute conversion
amrex::Real CWinv = 1.0/CW;
y[0] = c[0]*1.008000*CWinv; 
y[1] = c[1]*2.016000*CWinv; 
y[2] = c[2]*15.999000*CWinv; 
y[3] = c[3]*17.007000*CWinv; 
y[4] = c[4]*18.015000*CWinv; 
y[5] = c[5]*31.998000*CWinv; 
y[6] = c[6]*28.014000*CWinv; 
y[7] = c[7]*33.006000*CWinv; 
y[8] = c[8]*34.014000*CWinv; 
y[9] = c[9]*17.007000*CWinv; 
y[10] = c[10]*39.950000*CWinv; 
y[11] = c[11]*4.002602*CWinv; 

}

// get specific heat at constant volume as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(const amrex::Real T, amrex::Real cvml[])
{
cv_R(cvml, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
cvml[id] *= 8.31446261815324e+07;
}
}

// get specific heat at constant pressure as a 
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(const amrex::Real T, amrex::Real cpml[])
{
cp_R(cpml, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
cpml[id] *= 8.31446261815324e+07;
}
}

// get internal energy as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(const amrex::Real T, amrex::Real uml[])
{
amrex::Real RT = 8.31446261815324e+07*T; // R*T
speciesInternalEnergy(uml, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
uml[id] *= RT;
}
}

// get enthalpy as a function 
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(const amrex::Real T, amrex::Real hml[])
{
amrex::Real RT = 8.31446261815324e+07*T; // R*T
speciesEnthalpy(hml, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
hml[id] *= RT;
}
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(const amrex::Real T, amrex::Real sml[])
{
speciesEntropy(sml, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
sml[id] *= 8.31446261815324e+07;
}
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
cv_R(cvms, T);
// multiply by R/molecularweight
cvms[0] *= 8.248474819596468e+07; // H
cvms[1] *= 4.124237409798234e+07; // H2
cvms[2] *= 5.196863940342046e+06; // O
cvms[3] *= 4.888847308845322e+06; // OH
cvms[4] *= 4.615299815794193e+06; // H2O
cvms[5] *= 2.598431970171023e+06; // O2
cvms[6] *= 2.967966951578939e+06; // N2
cvms[7] *= 2.519076112874398e+06; // HO2
cvms[8] *= 2.444423654422661e+06; // H2O2
cvms[9] *= 4.888847308845322e+06; // OHV
cvms[10] *= 2.081217176008320e+06; // AR
cvms[11] *= 2.077264394049980e+07; // HE
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
cp_R(cpms, T);
// multiply by R/molecularweight
cpms[0] *= 8.248474819596468e+07; // H
cpms[1] *= 4.124237409798234e+07; // H2
cpms[2] *= 5.196863940342046e+06; // O
cpms[3] *= 4.888847308845322e+06; // OH
cpms[4] *= 4.615299815794193e+06; // H2O
cpms[5] *= 2.598431970171023e+06; // O2
cpms[6] *= 2.967966951578939e+06; // N2
cpms[7] *= 2.519076112874398e+06; // HO2
cpms[8] *= 2.444423654422661e+06; // H2O2
cpms[9] *= 4.888847308845322e+06; // OHV
cpms[10] *= 2.081217176008320e+06; // AR
cpms[11] *= 2.077264394049980e+07; // HE
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(const amrex::Real T, amrex::Real ums[])
{
amrex::Real RT = 8.31446261815324e+07*T; // R*T

speciesInternalEnergy(ums, T);

for (int i = 0; i < 12; i++)
{
ums[i] *= RT*imw(i);
}
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(const amrex::Real T, amrex::Real hms[])
{
amrex::Real RT = 8.31446261815324e+07*T; // R*T

speciesEnthalpy(hms, T);

for (int i = 0; i < 12; i++)
{
hms[i] *= RT*imw(i);
}
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(const amrex::Real T, amrex::Real sms[])
{
speciesEntropy(sms, T);
// multiply by R/molecularweight
sms[0] *= 8.248474819596468e+07; // H
sms[1] *= 4.124237409798234e+07; // H2
sms[2] *= 5.196863940342046e+06; // O
sms[3] *= 4.888847308845322e+06; // OH
sms[4] *= 4.615299815794193e+06; // H2O
sms[5] *= 2.598431970171023e+06; // O2
sms[6] *= 2.967966951578939e+06; // N2
sms[7] *= 2.519076112874398e+06; // HO2
sms[8] *= 2.444423654422661e+06; // H2O2
sms[9] *= 4.888847308845322e+06; // OHV
sms[10] *= 2.081217176008320e+06; // AR
sms[11] *= 2.077264394049980e+07; // HE
}

// GPU version of productionRate: no more use of thermo namespace vectors 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, const amrex::Real * sc, const amrex::Real * /*sc_qss*/,const amrex::Real T, const amrex::Real invT, const amrex::Real logT)
{

// reaction 18: H2O2 (+M) <=> 2 OH (+M)
qf[0] = sc[8];
qr[0] = (sc[3] * sc[3]);

// reaction 30: H + O2 (+M) <=> HO2 (+M)
qf[1] = sc[0]*sc[5];
qr[1] = sc[7];

// reaction 31: H + O2 (+AR) <=> HO2 (+AR)
qf[2] = sc[0]*sc[5]*sc[10];
qr[2] = sc[7]*sc[10];

// reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)
qf[3] = sc[0]*sc[4]*sc[5];
qr[3] = sc[4]*sc[7];

// reaction 33: H + O2 (+H2) <=> HO2 (+H2)
qf[4] = sc[0]*sc[1]*sc[5];
qr[4] = sc[1]*sc[7];

// reaction 34: H + O2 (+HE) <=> HO2 (+HE)
qf[5] = sc[0]*sc[5]*sc[11];
qr[5] = sc[7]*sc[11];

// reaction 0: H2 + M <=> 2 H + M
qf[6] = sc[1];
qr[6] = (sc[0] * sc[0]);

// reaction 1: 2 O + M <=> O2 + M
qf[7] = (sc[2] * sc[2]);
qr[7] = sc[5];

// reaction 5: H + OH + M <=> H2O + M
qf[8] = sc[0]*sc[3];
qr[8] = sc[4];

// reaction 6: H + OH + HE <=> H2O + HE
qf[9] = sc[0]*sc[3]*sc[11];
qr[9] = sc[4]*sc[11];

// reaction 8: H + O + M <=> OH + M
qf[10] = sc[0]*sc[2];
qr[10] = sc[3];

// reaction 9: H + O + M <=> OHV + M
qf[11] = sc[0]*sc[2];
qr[11] = sc[9];

// reaction 2: H2 + O <=> H + OH
qf[12] = sc[1]*sc[2];
qr[12] = sc[0]*sc[3];

// reaction 3: H2 + OH <=> H + H2O
qf[13] = sc[1]*sc[3];
qr[13] = sc[0]*sc[4];

// reaction 4: H + O2 <=> O + OH
qf[14] = sc[0]*sc[5];
qr[14] = sc[2]*sc[3];

// reaction 7: H2O + O <=> 2 OH
qf[15] = sc[2]*sc[4];
qr[15] = (sc[3] * sc[3]);

// reaction 10: H2O + OHV <=> H2O + OH
qf[16] = sc[4]*sc[9];
qr[16] = sc[3]*sc[4];

// reaction 11: H2 + OHV <=> H2 + OH
qf[17] = sc[1]*sc[9];
qr[17] = sc[1]*sc[3];

// reaction 12: N2 + OHV <=> N2 + OH
qf[18] = sc[6]*sc[9];
qr[18] = sc[3]*sc[6];

// reaction 13: OH + OHV <=> 2 OH
qf[19] = sc[3]*sc[9];
qr[19] = (sc[3] * sc[3]);

// reaction 14: H + OHV <=> H + OH
qf[20] = sc[0]*sc[9];
qr[20] = sc[0]*sc[3];

// reaction 15: AR + OHV <=> AR + OH
qf[21] = sc[9]*sc[10];
qr[21] = sc[3]*sc[10];

// reaction 16: OHV <=> OH
qf[22] = sc[9];
qr[22] = sc[3];

// reaction 17: O2 + OHV <=> O2 + OH
qf[23] = sc[5]*sc[9];
qr[23] = sc[3]*sc[5];

// reaction 19: H + H2O2 <=> H2O + OH
qf[24] = sc[0]*sc[8];
qr[24] = sc[3]*sc[4];

// reaction 20: H + H2O2 <=> H2 + HO2
qf[25] = sc[0]*sc[8];
qr[25] = sc[1]*sc[7];

// reaction 21: H2O2 + O <=> HO2 + OH
qf[26] = sc[2]*sc[8];
qr[26] = sc[3]*sc[7];

// reaction 22: H2O2 + OH <=> H2O + HO2
qf[27] = sc[3]*sc[8];
qr[27] = sc[4]*sc[7];

// reaction 23: H2O2 + OH <=> H2O + HO2
qf[28] = sc[3]*sc[8];
qr[28] = sc[4]*sc[7];

// reaction 24: H + HO2 <=> 2 OH
qf[29] = sc[0]*sc[7];
qr[29] = (sc[3] * sc[3]);

// reaction 25: H + HO2 <=> H2 + O2
qf[30] = sc[0]*sc[7];
qr[30] = sc[1]*sc[5];

// reaction 26: HO2 + O <=> O2 + OH
qf[31] = sc[2]*sc[7];
qr[31] = sc[3]*sc[5];

// reaction 27: HO2 + OH <=> H2O + O2
qf[32] = sc[3]*sc[7];
qr[32] = sc[4]*sc[5];

// reaction 28: 2 HO2 <=> H2O2 + O2
qf[33] = (sc[7] * sc[7]);
qr[33] = sc[5]*sc[8];

// reaction 29: 2 HO2 <=> H2O2 + O2
qf[34] = (sc[7] * sc[7]);
qr[34] = sc[5]*sc[8];

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int i = 0; i < 12; ++i) {
mixture += sc[i];
}

// compute the Gibbs free energy
amrex::Real g_RT[12];
gibbs(g_RT, T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
amrex::Real refC = 101325 / 8.31446 * invT;
amrex::Real refCinv = 1 / refC;

// Evaluate the kfs
amrex::Real k_f, Corr;
amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;


// reaction 0:  H2 + M <=> 2 H + M
k_f = 45800000000000
           * exp((-1.4) * logT - (52535.8199100127) * invT);
Corr  = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.17)*sc[11];
qf[6] *= Corr * k_f;
qr[6] *= Corr * k_f * exp(-(-2.000000*g_RT[0] + g_RT[1])) * (refCinv);
// reaction 1:  2 O + M <=> O2 + M
k_f = 6160
           * exp((-0.5) * logT);
Corr  = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.17)*sc[10] + (-0.17)*sc[11];
qf[7] *= Corr * k_f;
qr[7] *= Corr * k_f * exp(-(2.000000*g_RT[2] - g_RT[5])) * (refC);
// reaction 2:  H2 + O <=> H + OH
k_f = 0.0508
           * exp((2.7) * logT - (3150.13632793755) * invT);
qf[12] *= k_f;
qr[12] *= k_f * exp(-(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]));
// reaction 3:  H2 + OH <=> H + H2O
k_f = 43800000
           * exp(-(3517.48449397499) * invT);
qf[13] *= k_f;
qr[13] *= k_f * exp(-(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]));
// reaction 4:  H + O2 <=> O + OH
k_f = 104000000
           * exp(-(7803.88405329959) * invT);
qf[14] *= k_f;
qr[14] *= k_f * exp(-(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]));
// reaction 5:  H + OH + M <=> H2O + M
k_f = 35000000000
           * exp((-2) * logT);
Corr  = mixture + (-0.27)*sc[1] + (2.65)*sc[4] + (-0.62)*sc[10] + (-1)*sc[11];
qf[8] *= Corr * k_f;
qr[8] *= Corr * k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC);
// reaction 6:  H + OH + HE <=> H2O + HE
k_f = 404000000000000
           * exp((-3.234) * logT);
qf[9] *= k_f;
qr[9] *= k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC);
// reaction 7:  H2O + O <=> 2 OH
k_f = 2.97
           * exp((2.02) * logT - (6743.10332178324) * invT);
qf[15] *= k_f;
qr[15] *= k_f * exp(-(g_RT[2] - 2.000000*g_RT[3] + g_RT[4]));
// reaction 8:  H + O + M <=> OH + M
k_f = 4710000
           * exp((-1) * logT);
Corr  = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.25)*sc[10] + (-0.25)*sc[11];
qf[10] *= Corr * k_f;
qr[10] *= Corr * k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[3])) * (refC);
// reaction 9:  H + O + M <=> OHV + M
k_f = 15
           * exp(-(3006.7195781832) * invT);
Corr  = mixture + (5.5)*sc[4] + (-0.6)*sc[5] + (-0.6)*sc[6] + (-0.65)*sc[10];
qf[11] *= Corr * k_f;
qr[11] *= Corr * k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[9])) * (refC);
// reaction 10:  H2O + OHV <=> H2O + OH
k_f = 5930000
           * exp((0.5) * logT - (-432.766332592059) * invT);
qf[16] *= k_f;
qr[16] *= k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[4] + g_RT[9]));
// reaction 11:  H2 + OHV <=> H2 + OH
k_f = 2950000
           * exp((0.5) * logT - (-223.428199617296) * invT);
qf[17] *= k_f;
qr[17] *= k_f * exp(-(g_RT[1] - g_RT[1] - g_RT[3] + g_RT[9]));
// reaction 12:  N2 + OHV <=> N2 + OH
k_f = 108000
           * exp((0.5) * logT - (-624.995098929462) * invT);
qf[18] *= k_f;
qr[18] *= k_f * exp(-(-g_RT[3] + g_RT[6] - g_RT[6] + g_RT[9]));
// reaction 13:  OH + OHV <=> 2 OH
k_f = 6010000
           * exp((0.5) * logT - (-384.457532674806) * invT);
qf[19] *= k_f;
qr[19] *= k_f * exp(-(g_RT[3] - 2.000000*g_RT[3] + g_RT[9]));
// reaction 14:  H + OHV <=> H + OH
k_f = 1310000
           * exp((0.5) * logT - (-84.0371831893882) * invT);
qf[20] *= k_f;
qr[20] *= k_f * exp(-(g_RT[0] - g_RT[0] - g_RT[3] + g_RT[9]));
// reaction 15:  AR + OHV <=> AR + OH
k_f = 1690000
           * exp(-(2080.80091310252) * invT);
qf[21] *= k_f;
qr[21] *= k_f * exp(-(-g_RT[3] + g_RT[9] + g_RT[10] - g_RT[10]));
// reaction 16:  OHV <=> OH
k_f = 1450000
           ;
qf[22] *= k_f;
qr[22] *= k_f * exp(-(-g_RT[3] + g_RT[9]));
// reaction 17:  O2 + OHV <=> O2 + OH
k_f = 2100000
           * exp((0.5) * logT - (-240.537566254656) * invT);
qf[23] *= k_f;
qr[23] *= k_f * exp(-(-g_RT[3] + g_RT[5] - g_RT[5] + g_RT[9]));
// reaction 18:  H2O2 (+M) <=> 2 OH (+M)
k_f = 2000000000000
           * exp((0.9) * logT - (24531.3092413143) * invT);
Corr  = mixture + (2.7)*sc[1] + (4.1)*sc[4] + (-0.21)*sc[5] + (4.2)*sc[8] + (-0.32)*sc[10] + (-0.56)*sc[11];
redP = Corr / k_f * 2.35e+18 
           * exp(-2.293  * logT - (24528.2899413195) *invT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    (0.57)*exp(-T * 1e+30)
    + 0.43 * exp(-T * 1e-30)
    + 0.0);
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[0] *= Corr * k_f;
qr[0] *= Corr * k_f * exp(-(-2.000000*g_RT[3] + g_RT[8])) * (refCinv);
// reaction 19:  H + H2O2 <=> H2O + OH
k_f = 24100000
           * exp(-(1997.77016324474) * invT);
qf[24] *= k_f;
qr[24] *= k_f * exp(-(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[8]));
// reaction 20:  H + H2O2 <=> H2 + HO2
k_f = 21500
           * exp((1) * logT - (3019.29999482832) * invT);
qf[25] *= k_f;
qr[25] *= k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT[8]));
// reaction 21:  H2O2 + O <=> HO2 + OH
k_f = 9.55
           * exp((2) * logT - (1997.77016324474) * invT);
qf[26] *= k_f;
qr[26] *= k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[8]));
// reaction 22:  H2O2 + OH <=> H2O + HO2
k_f = 1740000
           * exp(-(160.022899725901) * invT);
qf[27] *= k_f;
qr[27] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]));
// reaction 23:  H2O2 + OH <=> H2O + HO2
k_f = 75900000
           * exp(-(3657.88194373451) * invT);
qf[28] *= k_f;
qr[28] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]));
// reaction 24:  H + HO2 <=> 2 OH
k_f = 70800000
           * exp(-(148.448916412392) * invT);
qf[29] *= k_f;
qr[29] *= k_f * exp(-(g_RT[0] - 2.000000*g_RT[3] + g_RT[7]));
// reaction 25:  H + HO2 <=> H2 + O2
k_f = 11400
           * exp((1.083) * logT - (278.782032855815) * invT);
qf[30] *= k_f;
qr[30] *= k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[7]));
// reaction 26:  HO2 + O <=> O2 + OH
k_f = 32500000
           ;
qf[31] *= k_f;
qr[31] *= k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[7]));
// reaction 27:  HO2 + OH <=> H2O + O2
k_f = 31648300
           * exp((1.380154522) * logT - (15400.9555526519) * invT);
qf[32] *= k_f;
qr[32] *= k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7]));
// reaction 28:  2 HO2 <=> H2O2 + O2
k_f = 12100
           * exp((0.422) * logT - (-745.26388205679) * invT);
qf[33] *= k_f;
qr[33] *= k_f * exp(-(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]));
// reaction 29:  2 HO2 <=> H2O2 + O2
k_f = 16900000000
           * exp((-0.681) * logT - (6507.59792218663) * invT);
qf[34] *= k_f;
qr[34] *= k_f * exp(-(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]));
// reaction 30:  H + O2 (+M) <=> HO2 (+M)
k_f = 1025000
           * exp((0.604) * logT - (-121.325538125518) * invT);
Corr  = mixture + (-1)*sc[1] + (-1)*sc[4] + (-1)*sc[10] + (-1)*sc[11];
redP = Corr / k_f * 27200000 
           * exp(-1.29  * logT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    (0.542)*exp(-T * 10000000000)
    + 0.458 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[1] *= Corr * k_f;
qr[1] *= Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC);
// reaction 31:  H + O2 (+AR) <=> HO2 (+AR)
k_f = 564000
           * exp((0.691) * logT - (-359.497986050892) * invT);
Corr  = mixture;
redP = Corr / k_f * 1190000 
           * exp(-1.04  * logT - (0.0553538332385192) *invT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    (0.565)*exp(-T * 10000000000)
    + 0.435 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[2] *= Corr * k_f;
qr[2] *= Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC);
// reaction 32:  H + O2 (+H2O) <=> HO2 (+H2O)
k_f = 2150000
           * exp((0.292) * logT - (-297.904266156394) * invT);
Corr  = mixture;
redP = Corr / k_f * 22800000 
           * exp(-0.977  * logT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    + exp(-T * 1e-30)
    + exp(-1e+30 * invT));
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[3] *= Corr * k_f;
qr[3] *= Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC);
// reaction 33:  H + O2 (+H2) <=> HO2 (+H2)
k_f = 4230000
           * exp((0.646) * logT - (702.339500463647) * invT);
Corr  = mixture;
redP = Corr / k_f * 828000000 
           * exp(-1.72  * logT - (-1.36371716433079) *invT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    (0.647)*exp(-T * 10000000000)
    + 0.353 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[4] *= Corr * k_f;
qr[4] *= Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC);
// reaction 34:  H + O2 (+HE) <=> HO2 (+HE)
k_f = 1025000
           * exp((0.604) * logT - (-121.325538125518) * invT);
Corr  = mixture;
redP = Corr / k_f * 27200000 
           * exp(-1.29  * logT);
F = redP / (1.0 + redP);
logPred = log10(redP);
logFcent = log10(
    (0.542)*exp(-T * 10000000000)
    + 0.458 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
troe_c = -0.4 - 0.67 * logFcent;
troe_n = 0.75 - 1.27 * logFcent;
troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));
Corr = F * F_troe;
qf[5] *= Corr * k_f;
qr[5] *= Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC);


}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, const amrex::Real * sc, const amrex::Real T)
{
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
const amrex::Real refC = 101325 / 8.31446 * invT;
const amrex::Real refCinv = 1 / refC;

for (int i = 0; i < 12; ++i) {
wdot[i] = 0.0;
}

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int i = 0; i < 12; ++i) {
mixture += sc[i];
}

// compute the Gibbs free energy
amrex::Real g_RT[12];
gibbs(g_RT, T);

{
// reaction 18:  H2O2 (+M) <=> 2 OH (+M)
const amrex::Real k_f = 2000000000000
           * exp((0.9) * logT - (24531.3092413143) * invT);
amrex::Real Corr = mixture + (2.7)*sc[1] + (4.1)*sc[4] + (-0.21)*sc[5] + (4.2)*sc[8] + (-0.32)*sc[10] + (-0.56)*sc[11];
const amrex::Real redP = Corr / k_f * 2.35e+18 
           * exp(-2.293 * logT - 24528.2899413195 * invT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    0.57 * exp(-T * 1e+30)
    + 0.43 * exp(-T * 1e-30)
    + 0.0);
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[8]);
const amrex::Real qr = Corr * k_f * exp(-(-2.000000*g_RT[3] + g_RT[8])) * (refCinv) * ((sc[3] * sc[3]));
const amrex::Real qdot = qf - qr;
wdot[3] += 2.000000 * qdot;
wdot[8] -= qdot;
}

{
// reaction 30:  H + O2 (+M) <=> HO2 (+M)
const amrex::Real k_f = 1025000
           * exp((0.604) * logT - (-121.325538125518) * invT);
amrex::Real Corr = mixture + (-1)*sc[1] + (-1)*sc[4] + (-1)*sc[10] + (-1)*sc[11];
const amrex::Real redP = Corr / k_f * 27200000 
           * exp(-1.29 * logT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    0.542 * exp(-T * 10000000000)
    + 0.458 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[0]*sc[5]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] -= qdot;
wdot[7] += qdot;
}

{
// reaction 31:  H + O2 (+AR) <=> HO2 (+AR)
const amrex::Real k_f = 564000
           * exp((0.691) * logT - (-359.497986050892) * invT);
amrex::Real Corr = mixture;
const amrex::Real redP = Corr / k_f * 1190000 
           * exp(-1.04 * logT - 0.0553538332385192 * invT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    0.565 * exp(-T * 10000000000)
    + 0.435 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[0]*sc[5]*sc[10]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7]*sc[10]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] -= qdot;
wdot[7] += qdot;
}

{
// reaction 32:  H + O2 (+H2O) <=> HO2 (+H2O)
const amrex::Real k_f = 2150000
           * exp((0.292) * logT - (-297.904266156394) * invT);
amrex::Real Corr = std::max(sc[4],1.e-16); //sc[4];  //mixture;
const amrex::Real redP = Corr / k_f * 22800000 
           * exp(-0.977 * logT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    + 1 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[0]*sc[5]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (sc[7]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] -= qdot;
wdot[7] += qdot;
}

{
// reaction 33:  H + O2 (+H2) <=> HO2 (+H2)
const amrex::Real k_f = 4230000
           * exp((0.646) * logT - (702.339500463647) * invT);
amrex::Real Corr = std::max(sc[1],1.e-16); //mixture;
const amrex::Real redP = Corr / k_f * 828000000 
           * exp(-1.72 * logT - -1.36371716433079 * invT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    0.647 * exp(-T * 10000000000)
    + 0.353 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[0]*sc[5]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (sc[7]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] -= qdot;
wdot[7] += qdot;
}

{
// reaction 34:  H + O2 (+HE) <=> HO2 (+HE)
const amrex::Real k_f = 1025000
           * exp((0.604) * logT - (-121.325538125518) * invT);
amrex::Real Corr = mixture;
const amrex::Real redP = Corr / k_f * 27200000 
           * exp(-1.29 * logT);
const amrex::Real F = redP / (1.0 + redP);
const amrex::Real logPred = log10(redP);
const amrex::Real logFcent = log10(
    0.542 * exp(-T * 10000000000)
    + 0.458 * exp(-T * 1e-30)
    + exp(-1e+30 * invT));
const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
const amrex::Real F_troe = exp(M_LN10 *logFcent / (1.0 + troe * troe));
Corr = F * F_troe;
const amrex::Real qf = Corr * k_f * (sc[0]*sc[5]*sc[11]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7]*sc[11]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[5] -= qdot;
wdot[7] += qdot;
}

{
// reaction 0:  H2 + M <=> 2 H + M
const amrex::Real k_f = 45800000000000
           * exp((-1.4) * logT - (52535.8199100127) * invT);
const amrex::Real Corr = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.17)*sc[11];
const amrex::Real qf = Corr * k_f * (sc[1]);
const amrex::Real qr = Corr * k_f * exp(-(-2.000000*g_RT[0] + g_RT[1])) * (refCinv) * ((sc[0] * sc[0]));
const amrex::Real qdot = qf - qr;
wdot[0] += 2.000000 * qdot;
wdot[1] -= qdot;
}

{
// reaction 1:  2 O + M <=> O2 + M
const amrex::Real k_f = 6160
           * exp((-0.5) * logT);
const amrex::Real Corr = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.17)*sc[10] + (-0.17)*sc[11];
const amrex::Real qf = Corr * k_f * ((sc[2] * sc[2]));
const amrex::Real qr = Corr * k_f * exp(-(2.000000*g_RT[2] - g_RT[5])) * (refC) * (sc[5]);
const amrex::Real qdot = qf - qr;
wdot[2] -= 2.000000 * qdot;
wdot[5] += qdot;
}

{
// reaction 5:  H + OH + M <=> H2O + M
const amrex::Real k_f = 35000000000
           * exp((-2) * logT);
const amrex::Real Corr = mixture + (-0.27)*sc[1] + (2.65)*sc[4] + (-0.62)*sc[10] + (-1)*sc[11];
const amrex::Real qf = Corr * k_f * (sc[0]*sc[3]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC) * (sc[4]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] -= qdot;
wdot[4] += qdot;
}

{
// reaction 6:  H + OH + HE <=> H2O + HE
const amrex::Real k_f = 404000000000000
           * exp((-3.234) * logT);
const amrex::Real qf = k_f * (sc[0]*sc[3]*sc[11]);
const amrex::Real qr = k_f * exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC) * (sc[4]*sc[11]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] -= qdot;
wdot[4] += qdot;
}

{
// reaction 8:  H + O + M <=> OH + M
const amrex::Real k_f = 4710000
           * exp((-1) * logT);
const amrex::Real Corr = mixture + (1.5)*sc[1] + (11)*sc[4] + (-0.25)*sc[10] + (-0.25)*sc[11];
const amrex::Real qf = Corr * k_f * (sc[0]*sc[2]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[3])) * (refC) * (sc[3]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[2] -= qdot;
wdot[3] += qdot;
}

{
// reaction 9:  H + O + M <=> OHV + M
const amrex::Real k_f = 15
           * exp(-(3006.7195781832) * invT);
const amrex::Real Corr = mixture + (5.5)*sc[4] + (-0.6)*sc[5] + (-0.6)*sc[6] + (-0.65)*sc[10];
const amrex::Real qf = Corr * k_f * (sc[0]*sc[2]);
const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] + g_RT[2] - g_RT[9])) * (refC) * (sc[9]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[2] -= qdot;
wdot[9] += qdot;
}

{
// reaction 2:  H2 + O <=> H + OH
const amrex::Real k_f = 0.0508
           * exp((2.7) * logT - (3150.13632793755) * invT);
const amrex::Real qf = k_f * (sc[1]*sc[2]);
const amrex::Real qr = k_f * exp(-(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3])) * (sc[0]*sc[3]);
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[1] -= qdot;
wdot[2] -= qdot;
wdot[3] += qdot;
}

{
// reaction 3:  H2 + OH <=> H + H2O
const amrex::Real k_f = 43800000
           * exp(-(3517.48449397499) * invT);
const amrex::Real qf = k_f * (sc[1]*sc[3]);
const amrex::Real qr = k_f * exp(-(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4])) * (sc[0]*sc[4]);
const amrex::Real qdot = qf - qr;
wdot[0] += qdot;
wdot[1] -= qdot;
wdot[3] -= qdot;
wdot[4] += qdot;
}

{
// reaction 4:  H + O2 <=> O + OH
const amrex::Real k_f = 104000000
           * exp(-(7803.88405329959) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[5]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5])) * (sc[2]*sc[3]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[2] += qdot;
wdot[3] += qdot;
wdot[5] -= qdot;
}

{
// reaction 7:  H2O + O <=> 2 OH
const amrex::Real k_f = 2.97
           * exp((2.02) * logT - (6743.10332178324) * invT);
const amrex::Real qf = k_f * (sc[2]*sc[4]);
const amrex::Real qr = k_f * exp(-(g_RT[2] - 2.000000*g_RT[3] + g_RT[4])) * ((sc[3] * sc[3]));
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[3] += 2.000000 * qdot;
wdot[4] -= qdot;
}

{
// reaction 10:  H2O + OHV <=> H2O + OH
const amrex::Real k_f = 5930000
           * exp((0.5) * logT - (-432.766332592059) * invT);
const amrex::Real qf = k_f * (sc[4]*sc[9]);
const amrex::Real qr = k_f * exp(-(-g_RT[3] + g_RT[4] - g_RT[4] + g_RT[9])) * (sc[3]*sc[4]);
const amrex::Real qdot = qf - qr;
wdot[3] += qdot;
wdot[4] -= qdot;
wdot[4] += qdot;
wdot[9] -= qdot;
}

{
// reaction 11:  H2 + OHV <=> H2 + OH
const amrex::Real k_f = 2950000
           * exp((0.5) * logT - (-223.428199617296) * invT);
const amrex::Real qf = k_f * (sc[1]*sc[9]);
const amrex::Real qr = k_f * exp(-(g_RT[1] - g_RT[1] - g_RT[3] + g_RT[9])) * (sc[1]*sc[3]);
const amrex::Real qdot = qf - qr;
wdot[1] -= qdot;
wdot[1] += qdot;
wdot[3] += qdot;
wdot[9] -= qdot;
}

{
// reaction 12:  N2 + OHV <=> N2 + OH
const amrex::Real k_f = 108000
           * exp((0.5) * logT - (-624.995098929462) * invT);
const amrex::Real qf = k_f * (sc[6]*sc[9]);
const amrex::Real qr = k_f * exp(-(-g_RT[3] + g_RT[6] - g_RT[6] + g_RT[9])) * (sc[3]*sc[6]);
const amrex::Real qdot = qf - qr;
wdot[3] += qdot;
wdot[6] -= qdot;
wdot[6] += qdot;
wdot[9] -= qdot;
}

{
// reaction 13:  OH + OHV <=> 2 OH
const amrex::Real k_f = 6010000
           * exp((0.5) * logT - (-384.457532674806) * invT);
const amrex::Real qf = k_f * (sc[3]*sc[9]);
const amrex::Real qr = k_f * exp(-(g_RT[3] - 2.000000*g_RT[3] + g_RT[9])) * ((sc[3] * sc[3]));
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[3] += 2.000000 * qdot;
wdot[9] -= qdot;
}

{
// reaction 14:  H + OHV <=> H + OH
const amrex::Real k_f = 1310000
           * exp((0.5) * logT - (-84.0371831893882) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[9]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[0] - g_RT[3] + g_RT[9])) * (sc[0]*sc[3]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[0] += qdot;
wdot[3] += qdot;
wdot[9] -= qdot;
}

{
// reaction 15:  AR + OHV <=> AR + OH
const amrex::Real k_f = 1690000
           * exp(-(2080.80091310252) * invT);
const amrex::Real qf = k_f * (sc[9]*sc[10]);
const amrex::Real qr = k_f * exp(-(-g_RT[3] + g_RT[9] + g_RT[10] - g_RT[10])) * (sc[3]*sc[10]);
const amrex::Real qdot = qf - qr;
wdot[3] += qdot;
wdot[9] -= qdot;
wdot[10] -= qdot;
wdot[10] += qdot;
}

{
// reaction 16:  OHV <=> OH
const amrex::Real k_f = 1450000
           ;
const amrex::Real qf = k_f * (sc[9]);
const amrex::Real qr = k_f * exp(-(-g_RT[3] + g_RT[9])) * (sc[3]);
const amrex::Real qdot = qf - qr;
wdot[3] += qdot;
wdot[9] -= qdot;
}

{
// reaction 17:  O2 + OHV <=> O2 + OH
const amrex::Real k_f = 2100000
           * exp((0.5) * logT - (-240.537566254656) * invT);
const amrex::Real qf = k_f * (sc[5]*sc[9]);
const amrex::Real qr = k_f * exp(-(-g_RT[3] + g_RT[5] - g_RT[5] + g_RT[9])) * (sc[3]*sc[5]);
const amrex::Real qdot = qf - qr;
wdot[3] += qdot;
wdot[5] -= qdot;
wdot[5] += qdot;
wdot[9] -= qdot;
}

{
// reaction 19:  H + H2O2 <=> H2O + OH
const amrex::Real k_f = 24100000
           * exp(-(1997.77016324474) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[8]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[8])) * (sc[3]*sc[4]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] += qdot;
wdot[4] += qdot;
wdot[8] -= qdot;
}

{
// reaction 20:  H + H2O2 <=> H2 + HO2
const amrex::Real k_f = 21500
           * exp((1) * logT - (3019.29999482832) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[8]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT[8])) * (sc[1]*sc[7]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[1] += qdot;
wdot[7] += qdot;
wdot[8] -= qdot;
}

{
// reaction 21:  H2O2 + O <=> HO2 + OH
const amrex::Real k_f = 9.55
           * exp((2) * logT - (1997.77016324474) * invT);
const amrex::Real qf = k_f * (sc[2]*sc[8]);
const amrex::Real qr = k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[8])) * (sc[3]*sc[7]);
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[3] += qdot;
wdot[7] += qdot;
wdot[8] -= qdot;
}

{
// reaction 22:  H2O2 + OH <=> H2O + HO2
const amrex::Real k_f = 1740000
           * exp(-(160.022899725901) * invT);
const amrex::Real qf = k_f * (sc[3]*sc[8]);
const amrex::Real qr = k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4]*sc[7]);
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[4] += qdot;
wdot[7] += qdot;
wdot[8] -= qdot;
}

{
// reaction 23:  H2O2 + OH <=> H2O + HO2
const amrex::Real k_f = 75900000
           * exp(-(3657.88194373451) * invT);
const amrex::Real qf = k_f * (sc[3]*sc[8]);
const amrex::Real qr = k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4]*sc[7]);
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[4] += qdot;
wdot[7] += qdot;
wdot[8] -= qdot;
}

{
// reaction 24:  H + HO2 <=> 2 OH
const amrex::Real k_f = 70800000
           * exp(-(148.448916412392) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[7]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - 2.000000*g_RT[3] + g_RT[7])) * ((sc[3] * sc[3]));
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[3] += 2.000000 * qdot;
wdot[7] -= qdot;
}

{
// reaction 25:  H + HO2 <=> H2 + O2
const amrex::Real k_f = 11400
           * exp((1.083) * logT - (278.782032855815) * invT);
const amrex::Real qf = k_f * (sc[0]*sc[7]);
const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[7])) * (sc[1]*sc[5]);
const amrex::Real qdot = qf - qr;
wdot[0] -= qdot;
wdot[1] += qdot;
wdot[5] += qdot;
wdot[7] -= qdot;
}

{
// reaction 26:  HO2 + O <=> O2 + OH
const amrex::Real k_f = 32500000
           ;
const amrex::Real qf = k_f * (sc[2]*sc[7]);
const amrex::Real qr = k_f * exp(-(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[7])) * (sc[3]*sc[5]);
const amrex::Real qdot = qf - qr;
wdot[2] -= qdot;
wdot[3] += qdot;
wdot[5] += qdot;
wdot[7] -= qdot;
}

//Non-Arrhenius reaction: HO2 + OH <=> H2O + O2
//
{
const amrex::Real A = 59068.3359e-6; //bimolecular reactions are scaled by a factor of 1e-6
const amrex::Real beta = -2.354484303;
const amrex::Real T_ev = T/11595;
const amrex::Real E = 0.136645427;
const amrex::Real logTev = log(T_ev);

amrex::Real bi[9]={
	-81.9478771,  // b1
	-117.4423425, // b2
	-55.6034362,  // b3
	-13.26229425, // b4
	-1.862299153, // b5
	-0.162012483, // b6
	-0.008642048, // b7
	-0.000260552, // b8
	-3.42655E-06, // b9
};

amrex::Real sum_coeff = 0.0;
amrex::Real logTe_pow = 1.0;
for (int i = 0; i < 9; ++i) {
    sum_coeff += bi[i] * logTe_pow;
    logTe_pow *= logTev;
}

const amrex::Real k_f = A*pow(T_ev,beta)*exp(E/T_ev+sum_coeff);
const amrex::Real qf = k_f * (sc[3]*sc[7]);
const amrex::Real qr = k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7])) * (sc[4]*sc[5]);
const amrex::Real qdot = qf - qr;

wdot[3] -= qdot; //OH
wdot[4] += qdot; //H2O
wdot[5] += qdot; //O2
wdot[7] -= qdot; //HO2

}

{
// reaction 27:  HO2 + OH <=> H2O + O2
const amrex::Real k_f = 31648300
           * exp((1.380154522) * logT - (15400.9555526519) * invT);
const amrex::Real qf = k_f * (sc[3]*sc[7]);
const amrex::Real qr = k_f * exp(-(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7])) * (sc[4]*sc[5]);
const amrex::Real qdot = qf - qr;
wdot[3] -= qdot;
wdot[4] += qdot;
wdot[5] += qdot;
wdot[7] -= qdot;
}

{
// reaction 28:  2 HO2 <=> H2O2 + O2
const amrex::Real k_f = 12100
           * exp((0.422) * logT - (-745.26388205679) * invT);
const amrex::Real qf = k_f * ((sc[7] * sc[7]));
const amrex::Real qr = k_f * exp(-(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8])) * (sc[5]*sc[8]);
const amrex::Real qdot = qf - qr;
wdot[5] += qdot;
wdot[7] -= 2.000000 * qdot;
wdot[8] += qdot;
}

{
// reaction 29:  2 HO2 <=> H2O2 + O2
const amrex::Real k_f = 16900000000
           * exp((-0.681) * logT - (6507.59792218663) * invT);
const amrex::Real qf = k_f * ((sc[7] * sc[7]));
const amrex::Real qr = k_f * exp(-(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8])) * (sc[5]*sc[8]);
const amrex::Real qdot = qf - qr;
wdot[5] += qdot;
wdot[7] -= 2.000000 * qdot;
wdot[8] += qdot;
}


}


// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

// convert to SI
for (int id = 0; id < 12; ++id) {
C[id] *= 1.0e6;
}

// convert to chemkin units
productionRate(wdot, C, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
C[id] *= 1.0e-6;
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P, const amrex::Real T,const amrex::Real y[], amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real YOW = 0; 
amrex::Real PWORT; 

// Compute inverse of mean molecular wt first
for (int i = 0; i < 12; i++)
{
YOW += y[i]*imw(i);
}
// PW/RT (see Eq. 7)
PWORT = P/(YOW * 8.31446261815324e+07 * T); 
// multiply by 1e6 so c goes to SI
PWORT *= 1e6; 
// Now compute conversion (and go to SI)
for (int i = 0; i < 12; i++)
{
c[i] = PWORT * y[i]*imw(i);
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P, const amrex::Real T,const amrex::Real x[], amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho, const amrex::Real T, const amrex::Real y[],  amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage

// See Eq 8 with an extra 1e6 so c goes to SI
for (int i = 0; i < 12; i++)
{
c[i] = 1e6 * rho * y[i]*imw(i);
}

// call productionRate
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho, const amrex::Real T, const amrex::Real x[],  amrex::Real wdot[])
{
amrex::Real c[12]; // temporary storage
amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
amrex::Real ROW; 
// Compute mean molecular wt first
XW += x[0]*1.008000; // H
XW += x[1]*2.016000; // H2
XW += x[2]*15.999000; // O
XW += x[3]*17.007000; // OH
XW += x[4]*18.015000; // H2O
XW += x[5]*31.998000; // O2
XW += x[6]*28.014000; // N2
XW += x[7]*33.006000; // HO2
XW += x[8]*34.014000; // H2O2
XW += x[9]*17.007000; // OHV
XW += x[10]*39.950000; // AR
XW += x[11]*4.002602; // HE
// Extra 1e6 factor to take c to SI
ROW = 1e6*rho / XW;

// Compute conversion, see Eq 11
for (int id = 0; id < 12; ++id) {
c[id] = x[id]*ROW;
}

// convert to chemkin units
productionRate(wdot, c, T);

// convert to chemkin units
for (int id = 0; id < 12; ++id) {
wdot[id] *= 1.0e-6;
}
}

//  species unit charge number 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
kcharge[0] = 0;// H
kcharge[1] = 0;// H2
kcharge[2] = 0;// O
kcharge[3] = 0;// OH
kcharge[4] = 0;// H2O
kcharge[5] = 0;// O2
kcharge[6] = 0;// N2
kcharge[7] = 0;// HO2
kcharge[8] = 0;// H2O2
kcharge[9] = 0;// OHV
kcharge[10] = 0;// AR
kcharge[11] = 0;// HE
}

//  species charge per unit mass 
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

int kchrg[12];
CKCHRG(kchrg);

for (int id = 0; id < 12; ++id) {
zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
}
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian_precond(amrex::Real *  J, const amrex::Real *  sc, const amrex::Real T, const int HP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
for (int i=0; i<169; i++) {
J[i] = 0.0;
}

amrex::Real wdot[12];
for (auto& val : wdot) {
val = 0.0;
}

const amrex::Real invT = 1.0 / T;
const amrex::Real invT2 = invT * invT;
const amrex::Real logT = log(T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
amrex::Real refC = 101325 / 8.31446 / T;
amrex::Real refCinv = 1.0 / refC;

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int k = 0; k < 12; ++k) {
mixture += sc[k];
}

// compute the Gibbs free energy
amrex::Real g_RT[12];
gibbs(g_RT, T);

// compute the species enthalpy
amrex::Real h_RT[12];
speciesEnthalpy(h_RT, T);

amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
amrex::Real dqdci, dcdc_fac, dqdc[12];
amrex::Real Pr, fPr, F, k_0, logPr;
amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
const amrex::Real ln10 = log(10.0);
const amrex::Real log10e = 1.0/log(10.0);
// reaction 18: H2O2 (+M) <=> 2 OH (+M)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 3.7 - 1)*sc[1] + ( 5.1 - 1)*sc[4] + ( 0.79 - 1)*sc[5] + ( 5.2 - 1)*sc[8] + ( 0.68 - 1)*sc[10] + ( 0.44 - 1)*sc[11];
// forward
phi_f = sc[8];
k_f = 2000000000000
            * exp(0.9 * logT - (24531.3092413143) * invT);
dlnkfdT = 0.9 * invT + (24531.3092413143) * invT2;
// pressure-fall-off
k_0 = 2.35e+24 * exp(-2.293 * logT - (24528.2899413195) * invT);
Pr = 1e-6 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -2.293 * invT + (24528.2899413195) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.43)*exp(-T/1e-30);
Fcent2 = 0.43 * exp(-T/1e+30);
Fcent3 = 0.;
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-30
    -Fcent2/1e+30
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = refC * exp(-2.000000*g_RT[3] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[8]) + (2.000000*h_RT[3]) - 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[3] += 2 * q;// OH
wdot[8] -= q;// H2O2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac;
dqdc[1] = 3.7*dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
dqdc[4] = 5.1*dcdc_fac;
dqdc[5] = 0.79*dcdc_fac;
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac;
dqdc[8] = 5.2*dcdc_fac + k_f;
dqdc[9] = dcdc_fac;
dqdc[10] = 0.68*dcdc_fac;
dqdc[11] = 0.44*dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+3] += 2 * dqdc[k];
J[13*k+8] -= dqdc[k];
}
J[159] += 2 * dqdT;// dwdot[OH]/dT
J[164] -= dqdT;// dwdot[H2O2]/dT

// reaction 30: H + O2 (+M) <=> HO2 (+M)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 0 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 0 - 1)*sc[10] + ( 0 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[5];
k_f = 1025000
            * exp(0.604 * logT - (-121.325538125518) * invT);
dlnkfdT = 0.604 * invT + (-121.325538125518) * invT2;
// pressure-fall-off
k_0 = 2.72e+19 * exp(-1.29 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.29 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.458)*exp(-T/1e-10);
Fcent2 = 0.458 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac + k_f*sc[5];
dqdc[1] = 0.0;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = 0.0;
dqdc[5] = dcdc_fac + k_f*sc[0];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r;
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = 0.0;
dqdc[11] = 0.0;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 31: H + O2 (+AR) <=> HO2 (+AR)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[5]*sc[10];
k_f = 564000
            * exp(0.691 * logT - (-359.497986050892) * invT);
dlnkfdT = 0.691 * invT + (-359.497986050892) * invT2;
// pressure-fall-off
k_0 = 1.19e+18 * exp(-1.04 * logT - (0.0553538332385192) * invT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.04 * invT + (0.0553538332385192) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.435)*exp(-T/1e-10);
Fcent2 = 0.435 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7]*sc[10];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5] + h_RT[10]) + (h_RT[7] + h_RT[10]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac + k_f*sc[5]*sc[10];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[10];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[10];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[4]*sc[5];
k_f = 2150000
            * exp(0.292 * logT - (-297.904266156394) * invT);
dlnkfdT = 0.292 * invT + (-297.904266156394) * invT2;
// pressure-fall-off
k_0 = 2.28e+19 * exp(-0.977 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -0.977 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-1)*exp(-T/1e-10);
Fcent2 = 1 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[4]*sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[4] + h_RT[5]) + (h_RT[4] + h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac + k_f*sc[4]*sc[5];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[4];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[4];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 33: H + O2 (+H2) <=> HO2 (+H2)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[1]*sc[5];
k_f = 4230000
            * exp(0.646 * logT - (702.339500463647) * invT);
dlnkfdT = 0.646 * invT + (702.339500463647) * invT2;
// pressure-fall-off
k_0 = 8.28e+20 * exp(-1.72 * logT - (-1.36371716433079) * invT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.72 * invT + (-1.36371716433079) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.353)*exp(-T/1e-10);
Fcent2 = 0.353 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[1]*sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[1] + h_RT[5]) + (h_RT[1] + h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac + k_f*sc[1]*sc[5];
dqdc[1] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[1];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[1];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 34: H + O2 (+HE) <=> HO2 (+HE)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[5]*sc[11];
k_f = 1025000
            * exp(0.604 * logT - (-121.325538125518) * invT);
dlnkfdT = 0.604 * invT + (-121.325538125518) * invT2;
// pressure-fall-off
k_0 = 2.72e+19 * exp(-1.29 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.29 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.458)*exp(-T/1e-10);
Fcent2 = 0.458 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7]*sc[11];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5] + h_RT[11]) + (h_RT[7] + h_RT[11]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = 0.0;
dqdc[0] = dcdc_fac + k_f*sc[5]*sc[11];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[11];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[11];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 0: H2 + M <=> 2 H + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.83 - 1)*sc[11];
// forward
phi_f = sc[1];
k_f = 45800000000000
            * exp(-1.4 * logT - (52535.8199100127) * invT);
dlnkfdT = -1.4 * invT + (52535.8199100127) * invT2;
// reverse
phi_r = (sc[0] * sc[0]);
Kc = refC * exp(-2.000000*g_RT[0] + g_RT[1]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[0]) - 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += 2 * q;// H
wdot[1] -= q;// H2
// for convenience
k_f *= alpha;
k_r *= alpha;
dqdc[0] = q_nocor - k_r*2.000000*sc[0];
dqdc[1] = 2.5*q_nocor + k_f;
dqdc[2] = q_nocor;
dqdc[3] = q_nocor;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = q_nocor;
dqdc[11] = 0.83*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] += 2 * dqdc[k];
J[13*k+1] -= dqdc[k];
}
J[156] += 2 * dqdT;// dwdot[H]/dT
J[157] -= dqdT;// dwdot[H2]/dT

// reaction 1: 2 O + M <=> O2 + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.83 - 1)*sc[10] + ( 0.83 - 1)*sc[11];
// forward
phi_f = (sc[2] * sc[2]);
k_f = 6160
            * exp(-0.5 * logT);
dlnkfdT = -0.5 * invT;
// reverse
phi_r = sc[5];
Kc = refCinv * exp(2.000000*g_RT[2] - g_RT[5]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= 2 * q;// O
wdot[5] += q;// O2
// for convenience
k_f *= alpha;
k_r *= alpha;
dqdc[0] = q_nocor;
dqdc[1] = 2.5*q_nocor;
dqdc[2] = q_nocor + k_f*2.000000*sc[2];
dqdc[3] = q_nocor;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor - k_r;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.83*q_nocor;
dqdc[11] = 0.83*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+2] += -2 * dqdc[k];
J[13*k+5] += dqdc[k];
}
J[158] += -2 * dqdT;// dwdot[O]/dT
J[161] += dqdT;// dwdot[O2]/dT

// reaction 5: H + OH + M <=> H2O + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 0.73 - 1)*sc[1] + ( 3.65 - 1)*sc[4] + ( 0.38 - 1)*sc[10] + ( 0 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[3];
k_f = 35000000000
            * exp(-2 * logT);
dlnkfdT = -2 * invT;
// reverse
phi_r = sc[4];
Kc = refCinv * exp(g_RT[0] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// for convenience
k_f *= alpha;
k_r *= alpha;
dqdc[0] = q_nocor + k_f*sc[3];
dqdc[1] = 0.73*q_nocor;
dqdc[2] = q_nocor;
dqdc[3] = q_nocor + k_f*sc[0];
dqdc[4] = 3.65*q_nocor - k_r;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.38*q_nocor;
dqdc[11] = 0.0;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+3] -= dqdc[k];
J[13*k+4] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[159] -= dqdT;// dwdot[OH]/dT
J[160] += dqdT;// dwdot[H2O]/dT

// reaction 6: H + OH + HE <=> H2O + HE
// a third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[3]*sc[11];
k_f = 404000000000000
            * exp(-3.234 * logT);
dlnkfdT = -3.234 * invT;
// reverse
phi_r = sc[4]*sc[11];
Kc = refCinv * exp(g_RT[0] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[3] + h_RT[11]) + (h_RT[4] + h_RT[11]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// d()/d[H]
dqdci =  + k_f*sc[3]*sc[11];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] -= dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
// d()/d[OH]
dqdci =  + k_f*sc[0]*sc[11];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[11];
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
// d()/d[HE]
dqdci =  + k_f*sc[0]*sc[3] - k_r*sc[4];
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[146] -= dqdci;              // dwdot[OH]/d[HE]
J[147] += dqdci;              // dwdot[H2O]/d[HE]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT

// reaction 8: H + O + M <=> OH + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[10] + ( 0.75 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[2];
k_f = 4710000
            * exp(-1 * logT);
dlnkfdT = -1 * invT;
// reverse
phi_r = sc[3];
Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] -= q;// O
wdot[3] += q;// OH
// for convenience
k_f *= alpha;
k_r *= alpha;
dqdc[0] = q_nocor + k_f*sc[2];
dqdc[1] = 2.5*q_nocor;
dqdc[2] = q_nocor + k_f*sc[0];
dqdc[3] = q_nocor - k_r;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.75*q_nocor;
dqdc[11] = 0.75*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+2] -= dqdc[k];
J[13*k+3] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[158] -= dqdT;// dwdot[O]/dT
J[159] += dqdT;// dwdot[OH]/dT

// reaction 9: H + O + M <=> OHV + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 6.5 - 1)*sc[4] + ( 0.4 - 1)*sc[5] + ( 0.4 - 1)*sc[6] + ( 0.35 - 1)*sc[10];
// forward
phi_f = sc[0]*sc[2];
k_f = 15
            * exp(- (3006.7195781832) * invT);
dlnkfdT = (3006.7195781832) * invT2;
// reverse
phi_r = sc[9];
Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[9]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] -= q;// O
wdot[9] += q;// OHV
// for convenience
k_f *= alpha;
k_r *= alpha;
dqdc[0] = q_nocor + k_f*sc[2];
dqdc[1] = q_nocor;
dqdc[2] = q_nocor + k_f*sc[0];
dqdc[3] = q_nocor;
dqdc[4] = 6.5*q_nocor;
dqdc[5] = 0.4*q_nocor;
dqdc[6] = 0.4*q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor - k_r;
dqdc[10] = 0.35*q_nocor;
dqdc[11] = q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+2] -= dqdc[k];
J[13*k+9] += dqdc[k];
}
J[156] -= dqdT;// dwdot[H]/dT
J[158] -= dqdT;// dwdot[O]/dT
J[165] += dqdT;// dwdot[OHV]/dT

// reaction 2: H2 + O <=> H + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[2];
k_f = 0.0508
            * exp(2.7 * logT - (3150.13632793755) * invT);
dlnkfdT = 2.7 * invT + (3150.13632793755) * invT2;
// reverse
phi_r = sc[0]*sc[3];
Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += q;// H
wdot[1] -= q;// H2
wdot[2] -= q;// O
wdot[3] += q;// OH
// d()/d[H]
dqdci =  - k_r*sc[3];
J[0] += dqdci;                // dwdot[H]/d[H]
J[1] -= dqdci;                // dwdot[H2]/d[H]
J[2] -= dqdci;                // dwdot[O]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
// d()/d[H2]
dqdci =  + k_f*sc[2];
J[13] += dqdci;               // dwdot[H]/d[H2]
J[14] -= dqdci;               // dwdot[H2]/d[H2]
J[15] -= dqdci;               // dwdot[O]/d[H2]
J[16] += dqdci;               // dwdot[OH]/d[H2]
// d()/d[O]
dqdci =  + k_f*sc[1];
J[26] += dqdci;               // dwdot[H]/d[O]
J[27] -= dqdci;               // dwdot[H2]/d[O]
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[0];
J[39] += dqdci;               // dwdot[H]/d[OH]
J[40] -= dqdci;               // dwdot[H2]/d[OH]
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
// d()/dT
J[156] += dqdT;               // dwdot[H]/dT
J[157] -= dqdT;               // dwdot[H2]/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT

// reaction 3: H2 + OH <=> H + H2O
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[3];
k_f = 43800000
            * exp(- (3517.48449397499) * invT);
dlnkfdT = (3517.48449397499) * invT2;
// reverse
phi_r = sc[0]*sc[4];
Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += q;// H
wdot[1] -= q;// H2
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// d()/d[H]
dqdci =  - k_r*sc[4];
J[0] += dqdci;                // dwdot[H]/d[H]
J[1] -= dqdci;                // dwdot[H2]/d[H]
J[3] -= dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
// d()/d[H2]
dqdci =  + k_f*sc[3];
J[13] += dqdci;               // dwdot[H]/d[H2]
J[14] -= dqdci;               // dwdot[H2]/d[H2]
J[16] -= dqdci;               // dwdot[OH]/d[H2]
J[17] += dqdci;               // dwdot[H2O]/d[H2]
// d()/d[OH]
dqdci =  + k_f*sc[1];
J[39] += dqdci;               // dwdot[H]/d[OH]
J[40] -= dqdci;               // dwdot[H2]/d[OH]
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[0];
J[52] += dqdci;               // dwdot[H]/d[H2O]
J[53] -= dqdci;               // dwdot[H2]/d[H2O]
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
// d()/dT
J[156] += dqdT;               // dwdot[H]/dT
J[157] -= dqdT;               // dwdot[H2]/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT

// reaction 4: H + O2 <=> O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[5];
k_f = 104000000
            * exp(- (7803.88405329959) * invT);
dlnkfdT = (7803.88405329959) * invT2;
// reverse
phi_r = sc[2]*sc[3];
Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] += q;// O
wdot[3] += q;// OH
wdot[5] -= q;// O2
// d()/d[H]
dqdci =  + k_f*sc[5];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[2] += dqdci;                // dwdot[O]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
// d()/d[O]
dqdci =  - k_r*sc[3];
J[26] -= dqdci;               // dwdot[H]/d[O]
J[28] += dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[31] -= dqdci;               // dwdot[O2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[2];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[41] += dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[44] -= dqdci;               // dwdot[O2]/d[OH]
// d()/d[O2]
dqdci =  + k_f*sc[0];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[67] += dqdci;               // dwdot[O]/d[O2]
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[158] += dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[161] -= dqdT;               // dwdot[O2]/dT

// reaction 7: H2O + O <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[4];
k_f = 2.97
            * exp(2.02 * logT - (6743.10332178324) * invT);
dlnkfdT = 2.02 * invT + (6743.10332178324) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[2] - 2.000000*g_RT[3] + g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += 2 * q;// OH
wdot[4] -= q;// H2O
// d()/d[O]
dqdci =  + k_f*sc[4];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += 2 * dqdci;           // dwdot[OH]/d[O]
J[30] -= dqdci;               // dwdot[H2O]/d[O]
// d()/d[OH]
dqdci =  - k_r*2.000000*sc[3];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += 2 * dqdci;           // dwdot[OH]/d[OH]
J[43] -= dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  + k_f*sc[2];
J[54] -= dqdci;               // dwdot[O]/d[H2O]
J[55] += 2 * dqdci;           // dwdot[OH]/d[H2O]
J[56] -= dqdci;               // dwdot[H2O]/d[H2O]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += 2 * dqdT;           // dwdot[OH]/dT
J[160] -= dqdT;               // dwdot[H2O]/dT

// reaction 10: H2O + OHV <=> H2O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[4]*sc[9];
k_f = 5930000
            * exp(0.5 * logT - (-432.766332592059) * invT);
dlnkfdT = 0.5 * invT + (-432.766332592059) * invT2;
// reverse
phi_r = sc[3]*sc[4];
Kc = exp(-g_RT[3] + g_RT[4] - g_RT[4] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[3] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[4];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[H2O]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[55] += dqdci;               // dwdot[OH]/d[H2O]
J[61] -= dqdci;               // dwdot[OHV]/d[H2O]
// d()/d[OHV]
dqdci =  + k_f*sc[4];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 11: H2 + OHV <=> H2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[9];
k_f = 2950000
            * exp(0.5 * logT - (-223.428199617296) * invT);
dlnkfdT = 0.5 * invT + (-223.428199617296) * invT2;
// reverse
phi_r = sc[1]*sc[3];
Kc = exp(g_RT[1] - g_RT[1] - g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[1] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[H2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[16] += dqdci;               // dwdot[OH]/d[H2]
J[22] -= dqdci;               // dwdot[OHV]/d[H2]
// d()/d[OH]
dqdci =  - k_r*sc[1];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[1];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 12: N2 + OHV <=> N2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[6]*sc[9];
k_f = 108000
            * exp(0.5 * logT - (-624.995098929462) * invT);
dlnkfdT = 0.5 * invT + (-624.995098929462) * invT2;
// reverse
phi_r = sc[3]*sc[6];
Kc = exp(-g_RT[3] + g_RT[6] - g_RT[6] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[6]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[6];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[N2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[81] += dqdci;               // dwdot[OH]/d[N2]
J[87] -= dqdci;               // dwdot[OHV]/d[N2]
// d()/d[OHV]
dqdci =  + k_f*sc[6];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 13: OH + OHV <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[9];
k_f = 6010000
            * exp(0.5 * logT - (-384.457532674806) * invT);
dlnkfdT = 0.5 * invT + (-384.457532674806) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[3] - 2.000000*g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  + k_f*sc[9] - k_r*2.000000*sc[3];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[3];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 14: H + OHV <=> H + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[9];
k_f = 1310000
            * exp(0.5 * logT - (-84.0371831893882) * invT);
dlnkfdT = 0.5 * invT + (-84.0371831893882) * invT2;
// reverse
phi_r = sc[0]*sc[3];
Kc = exp(g_RT[0] - g_RT[0] - g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[0] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[H]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[3] += dqdci;                // dwdot[OH]/d[H]
J[9] -= dqdci;                // dwdot[OHV]/d[H]
// d()/d[OH]
dqdci =  - k_r*sc[0];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[0];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 15: AR + OHV <=> AR + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[9]*sc[10];
k_f = 1690000
            * exp(- (2080.80091310252) * invT);
dlnkfdT = (2080.80091310252) * invT2;
// reverse
phi_r = sc[3]*sc[10];
Kc = exp(-g_RT[3] + g_RT[9] + g_RT[10] - g_RT[10]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[9] + h_RT[10]) + (h_RT[3] + h_RT[10]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[10];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[10];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/d[AR]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[133] += dqdci;              // dwdot[OH]/d[AR]
J[139] -= dqdci;              // dwdot[OHV]/d[AR]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 16: OHV <=> OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[9];
k_f = 1450000
           ;
dlnkfdT = 0.0;
// reverse
phi_r = sc[3];
Kc = exp(-g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[9]) + (h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r;
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f;
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 17: O2 + OHV <=> O2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[5]*sc[9];
k_f = 2100000
            * exp(0.5 * logT - (-240.537566254656) * invT);
dlnkfdT = 0.5 * invT + (-240.537566254656) * invT2;
// reverse
phi_r = sc[3]*sc[5];
Kc = exp(-g_RT[3] + g_RT[5] - g_RT[5] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[3] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[5];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[O2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[74] -= dqdci;               // dwdot[OHV]/d[O2]
// d()/d[OHV]
dqdci =  + k_f*sc[5];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 19: H + H2O2 <=> H2O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[8];
k_f = 24100000
            * exp(- (1997.77016324474) * invT);
dlnkfdT = (1997.77016324474) * invT2;
// reverse
phi_r = sc[3]*sc[4];
Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[3] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] += q;// OH
wdot[4] += q;// H2O
wdot[8] -= q;// H2O2
// d()/d[H]
dqdci =  + k_f*sc[8];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
J[8] -= dqdci;                // dwdot[H2O2]/d[H]
// d()/d[OH]
dqdci =  - k_r*sc[4];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[3];
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[55] += dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[H2O2]
dqdci =  + k_f*sc[0];
J[104] -= dqdci;              // dwdot[H]/d[H2O2]
J[107] += dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 20: H + H2O2 <=> H2 + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[8];
k_f = 21500
            * exp(1 * logT - (3019.29999482832) * invT);
dlnkfdT = 1 * invT + (3019.29999482832) * invT2;
// reverse
phi_r = sc[1]*sc[7];
Kc = exp(g_RT[0] - g_RT[1] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[1] += q;// H2
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[H]
dqdci =  + k_f*sc[8];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[1] += dqdci;                // dwdot[H2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
J[8] -= dqdci;                // dwdot[H2O2]/d[H]
// d()/d[H2]
dqdci =  - k_r*sc[7];
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[14] += dqdci;               // dwdot[H2]/d[H2]
J[20] += dqdci;               // dwdot[HO2]/d[H2]
J[21] -= dqdci;               // dwdot[H2O2]/d[H2]
// d()/d[HO2]
dqdci =  - k_r*sc[1];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[92] += dqdci;               // dwdot[H2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[0];
J[104] -= dqdci;              // dwdot[H]/d[H2O2]
J[105] += dqdci;              // dwdot[H2]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[157] += dqdT;               // dwdot[H2]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 21: H2O2 + O <=> HO2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[8];
k_f = 9.55
            * exp(2 * logT - (1997.77016324474) * invT);
dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
// reverse
phi_r = sc[3]*sc[7];
Kc = exp(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[3] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += q;// OH
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[O]
dqdci =  + k_f*sc[8];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[33] += dqdci;               // dwdot[HO2]/d[O]
J[34] -= dqdci;               // dwdot[H2O2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[7];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[HO2]
dqdci =  - k_r*sc[3];
J[93] -= dqdci;               // dwdot[O]/d[HO2]
J[94] += dqdci;               // dwdot[OH]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[2];
J[106] -= dqdci;              // dwdot[O]/d[H2O2]
J[107] += dqdci;              // dwdot[OH]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 22: H2O2 + OH <=> H2O + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[8];
k_f = 1740000
            * exp(- (160.022899725901) * invT);
dlnkfdT = (160.022899725901) * invT2;
// reverse
phi_r = sc[4]*sc[7];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[OH]
dqdci =  + k_f*sc[8];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[7];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[HO2]
dqdci =  - k_r*sc[4];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[3];
J[107] -= dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 23: H2O2 + OH <=> H2O + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[8];
k_f = 75900000
            * exp(- (3657.88194373451) * invT);
dlnkfdT = (3657.88194373451) * invT2;
// reverse
phi_r = sc[4]*sc[7];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[OH]
dqdci =  + k_f*sc[8];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[7];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[HO2]
dqdci =  - k_r*sc[4];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[3];
J[107] -= dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 24: H + HO2 <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[7];
k_f = 70800000
            * exp(- (148.448916412392) * invT);
dlnkfdT = (148.448916412392) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[0] - 2.000000*g_RT[3] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] += 2 * q;// OH
wdot[7] -= q;// HO2
// d()/d[H]
dqdci =  + k_f*sc[7];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] += 2 * dqdci;            // dwdot[OH]/d[H]
J[7] -= dqdci;                // dwdot[HO2]/d[H]
// d()/d[OH]
dqdci =  - k_r*2.000000*sc[3];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] += 2 * dqdci;           // dwdot[OH]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[HO2]
dqdci =  + k_f*sc[0];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[94] += 2 * dqdci;           // dwdot[OH]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] += 2 * dqdT;           // dwdot[OH]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 25: H + HO2 <=> H2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[7];
k_f = 11400
            * exp(1.083 * logT - (278.782032855815) * invT);
dlnkfdT = 1.083 * invT + (278.782032855815) * invT2;
// reverse
phi_r = sc[1]*sc[5];
Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[1] += q;// H2
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[H]
dqdci =  + k_f*sc[7];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[1] += dqdci;                // dwdot[H2]/d[H]
J[5] += dqdci;                // dwdot[O2]/d[H]
J[7] -= dqdci;                // dwdot[HO2]/d[H]
// d()/d[H2]
dqdci =  - k_r*sc[5];
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[14] += dqdci;               // dwdot[H2]/d[H2]
J[18] += dqdci;               // dwdot[O2]/d[H2]
J[20] -= dqdci;               // dwdot[HO2]/d[H2]
// d()/d[O2]
dqdci =  - k_r*sc[1];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[66] += dqdci;               // dwdot[H2]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[0];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[92] += dqdci;               // dwdot[H2]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[157] += dqdT;               // dwdot[H2]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 26: HO2 + O <=> O2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[7];
k_f = 32500000
           ;
dlnkfdT = 0.0;
// reverse
phi_r = sc[3]*sc[5];
Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += q;// OH
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[O]
dqdci =  + k_f*sc[7];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[31] += dqdci;               // dwdot[O2]/d[O]
J[33] -= dqdci;               // dwdot[HO2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[5];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[44] += dqdci;               // dwdot[O2]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[O2]
dqdci =  - k_r*sc[3];
J[67] -= dqdci;               // dwdot[O]/d[O2]
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[2];
J[93] -= dqdci;               // dwdot[O]/d[HO2]
J[94] += dqdci;               // dwdot[OH]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 27: HO2 + OH <=> H2O + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[7];
k_f = 31648300
            * exp(1.380154522 * logT - (15400.9555526519) * invT);
dlnkfdT = 1.380154522 * invT + (15400.9555526519) * invT2;
// reverse
phi_r = sc[4]*sc[5];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[OH]
dqdci =  + k_f*sc[7];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[44] += dqdci;               // dwdot[O2]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[5];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[57] += dqdci;               // dwdot[O2]/d[H2O]
J[59] -= dqdci;               // dwdot[HO2]/d[H2O]
// d()/d[O2]
dqdci =  - k_r*sc[4];
J[68] -= dqdci;               // dwdot[OH]/d[O2]
J[69] += dqdci;               // dwdot[H2O]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[3];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 28: 2 HO2 <=> H2O2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = (sc[7] * sc[7]);
k_f = 12100
            * exp(0.422 * logT - (-745.26388205679) * invT);
dlnkfdT = 0.422 * invT + (-745.26388205679) * invT2;
// reverse
phi_r = sc[5]*sc[8];
Kc = exp(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[7]) + (h_RT[5] + h_RT[8]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[5] += q;// O2
wdot[7] -= 2 * q;// HO2
wdot[8] += q;// H2O2
// d()/d[O2]
dqdci =  - k_r*sc[8];
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] += -2 * dqdci;          // dwdot[HO2]/d[O2]
J[73] += dqdci;               // dwdot[H2O2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*2.000000*sc[7];
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
J[99] += dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  - k_r*sc[5];
J[109] += dqdci;              // dwdot[O2]/d[H2O2]
J[111] += -2 * dqdci;         // dwdot[HO2]/d[H2O2]
J[112] += dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] += -2 * dqdT;          // dwdot[HO2]/dT
J[164] += dqdT;               // dwdot[H2O2]/dT

// reaction 29: 2 HO2 <=> H2O2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = (sc[7] * sc[7]);
k_f = 16900000000
            * exp(-0.681 * logT - (6507.59792218663) * invT);
dlnkfdT = -0.681 * invT + (6507.59792218663) * invT2;
// reverse
phi_r = sc[5]*sc[8];
Kc = exp(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[7]) + (h_RT[5] + h_RT[8]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[5] += q;// O2
wdot[7] -= 2 * q;// HO2
wdot[8] += q;// H2O2
// d()/d[O2]
dqdci =  - k_r*sc[8];
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] += -2 * dqdci;          // dwdot[HO2]/d[O2]
J[73] += dqdci;               // dwdot[H2O2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*2.000000*sc[7];
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
J[99] += dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  - k_r*sc[5];
J[109] += dqdci;              // dwdot[O2]/d[H2O2]
J[111] += -2 * dqdci;         // dwdot[HO2]/d[H2O2]
J[112] += dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] += -2 * dqdT;          // dwdot[HO2]/dT
J[164] += dqdT;               // dwdot[H2O2]/dT

amrex::Real c_R[12], dcRdT[12], e_RT[12];
amrex::Real * eh_RT;
if (HP == 1) {
cp_R(c_R, T);
dcvpRdT(dcRdT, T);
eh_RT = &h_RT[0];
}
else {
cv_R(c_R, T);
dcvpRdT(dcRdT, T);
speciesInternalEnergy(e_RT, T);
eh_RT = &e_RT[0];
}

amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
for (int k = 0; k < 12; ++k) {
cmix += c_R[k]*sc[k];
dcmixdT += dcRdT[k]*sc[k];
ehmix += eh_RT[k]*wdot[k];
dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[156+k];
}

amrex::Real cmixinv = 1.0/cmix;
amrex::Real tmp1 = ehmix*cmixinv;
amrex::Real tmp3 = cmixinv*T;
amrex::Real tmp2 = tmp1*tmp3;
amrex::Real dehmixdc;
// dTdot/d[X]
for (int k = 0; k < 12; ++k) {
dehmixdc = 0.0;
for (int m = 0; m < 12; ++m) {
dehmixdc += eh_RT[m]*J[k*13+m];
}
J[k*13+12] = tmp2*c_R[k] - tmp3*dehmixdc;
}
// dTdot/dT
J[168] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
#else
amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, const amrex::Real *  sc, const amrex::Real *  Tp, const int * HP)
{
amrex::Real c[12];

for (int k=0; k<12; k++) {
c[k] = 1.e6 * sc[k];
}

aJacobian_precond(J, c, *Tp, *HP);

// dwdot[k]/dT
// dTdot/d[X]
for (int k=0; k<12; k++) {
J[156+k] *= 1.e-6;
J[k*13+12] *= 1.e6;
}

}


// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, const amrex::Real * sc, const amrex::Real T, const int consP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
for (int i=0; i<169; i++) {
J[i] = 0.0;
}

amrex::Real wdot[12];
for (auto& val : wdot) {
val = 0.0;
}

const amrex::Real invT = 1.0 / T;
const amrex::Real invT2 = invT * invT;
const amrex::Real logT = log(T);

// reference concentration: P_atm / (RT) in inverse mol/m^3
amrex::Real refC = 101325 / 8.31446 / T;
amrex::Real refCinv = 1.0 / refC;

// compute the mixture concentration
amrex::Real mixture = 0.0;
for (int k = 0; k < 12; ++k) {
mixture += sc[k];
}

// compute the Gibbs free energy
amrex::Real g_RT[12];
gibbs(g_RT, T);

// compute the species enthalpy
amrex::Real h_RT[12];
speciesEnthalpy(h_RT, T);

amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
amrex::Real dqdci, dcdc_fac, dqdc[12];
amrex::Real Pr, fPr, F, k_0, logPr;
amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
const amrex::Real ln10 = log(10.0);
const amrex::Real log10e = 1.0/log(10.0);
// reaction 18: H2O2 (+M) <=> 2 OH (+M)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 3.7 - 1)*sc[1] + ( 5.1 - 1)*sc[4] + ( 0.79 - 1)*sc[5] + ( 5.2 - 1)*sc[8] + ( 0.68 - 1)*sc[10] + ( 0.44 - 1)*sc[11];
// forward
phi_f = sc[8];
k_f = 2000000000000
            * exp(0.9 * logT - (24531.3092413143) * invT);
dlnkfdT = 0.9 * invT + (24531.3092413143) * invT2;
// pressure-fall-off
k_0 = 2.35e+24 * exp(-2.293 * logT - (24528.2899413195) * invT);
Pr = 1e-6 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -2.293 * invT + (24528.2899413195) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.43)*exp(-T/1e-30);
Fcent2 = 0.43 * exp(-T/1e+30);
Fcent3 = 0.;
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-30
    -Fcent2/1e+30
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = refC * exp(-2.000000*g_RT[3] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[8]) + (2.000000*h_RT[3]) - 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[3] += 2 * q;// OH
wdot[8] -= q;// H2O2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H2]
dqdci = (3.7 - 1)*dcdc_fac;
J[16] += 2 * dqdci;           // dwdot[OH]/d[H2]
J[21] -= dqdci;               // dwdot[H2O2]/d[H2]
// d()/d[OH]
dqdci =  - k_r*2.000000*sc[3];
J[42] += 2 * dqdci;           // dwdot[OH]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci = (5.1 - 1)*dcdc_fac;
J[55] += 2 * dqdci;           // dwdot[OH]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[O2]
dqdci = (0.79 - 1)*dcdc_fac;
J[68] += 2 * dqdci;           // dwdot[OH]/d[O2]
J[73] -= dqdci;               // dwdot[H2O2]/d[O2]
// d()/d[H2O2]
dqdci = (5.2 - 1)*dcdc_fac + k_f;
J[107] += 2 * dqdci;          // dwdot[OH]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/d[AR]
dqdci = (0.68 - 1)*dcdc_fac;
J[133] += 2 * dqdci;          // dwdot[OH]/d[AR]
J[138] -= dqdci;              // dwdot[H2O2]/d[AR]
// d()/d[HE]
dqdci = (0.44 - 1)*dcdc_fac;
J[146] += 2 * dqdci;          // dwdot[OH]/d[HE]
J[151] -= dqdci;              // dwdot[H2O2]/d[HE]
}
else {
dqdc[0] = dcdc_fac;
dqdc[1] = 3.7*dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
dqdc[4] = 5.1*dcdc_fac;
dqdc[5] = 0.79*dcdc_fac;
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac;
dqdc[8] = 5.2*dcdc_fac + k_f;
dqdc[9] = dcdc_fac;
dqdc[10] = 0.68*dcdc_fac;
dqdc[11] = 0.44*dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+3] += 2 * dqdc[k];
J[13*k+8] -= dqdc[k];
}
}
J[159] += 2 * dqdT;// dwdot[OH]/dT
J[164] -= dqdT;// dwdot[H2O2]/dT

// reaction 30: H + O2 (+M) <=> HO2 (+M)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 0 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 0 - 1)*sc[10] + ( 0 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[5];
k_f = 1025000
            * exp(0.604 * logT - (-121.325538125518) * invT);
dlnkfdT = 0.604 * invT + (-121.325538125518) * invT2;
// pressure-fall-off
k_0 = 2.72e+19 * exp(-1.29 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.29 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.458)*exp(-T/1e-10);
Fcent2 = 0.458 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[5];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
// d()/d[H2]
dqdci = (0 - 1)*dcdc_fac;
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[18] -= dqdci;               // dwdot[O2]/d[H2]
J[20] += dqdci;               // dwdot[HO2]/d[H2]
// d()/d[H2O]
dqdci = (0 - 1)*dcdc_fac;
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[57] -= dqdci;               // dwdot[O2]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
// d()/d[O2]
dqdci =  + k_f*sc[0];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
J[72] += dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  - k_r;
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[96] -= dqdci;               // dwdot[O2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
// d()/d[AR]
dqdci = (0 - 1)*dcdc_fac;
J[130] -= dqdci;              // dwdot[H]/d[AR]
J[135] -= dqdci;              // dwdot[O2]/d[AR]
J[137] += dqdci;              // dwdot[HO2]/d[AR]
// d()/d[HE]
dqdci = (0 - 1)*dcdc_fac;
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[148] -= dqdci;              // dwdot[O2]/d[HE]
J[150] += dqdci;              // dwdot[HO2]/d[HE]
}
else {
dqdc[0] = dcdc_fac + k_f*sc[5];
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r;
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 31: H + O2 (+AR) <=> HO2 (+AR)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[5]*sc[10];
k_f = 564000
            * exp(0.691 * logT - (-359.497986050892) * invT);
dlnkfdT = 0.691 * invT + (-359.497986050892) * invT2;
// pressure-fall-off
k_0 = 1.19e+18 * exp(-1.04 * logT - (0.0553538332385192) * invT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.04 * invT + (0.0553538332385192) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.435)*exp(-T/1e-10);
Fcent2 = 0.435 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7]*sc[10];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5] + h_RT[10]) + (h_RT[7] + h_RT[10]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[5]*sc[10];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
// d()/d[O2]
dqdci =  + k_f*sc[0]*sc[10];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
J[72] += dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  - k_r*sc[10];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[96] -= dqdci;               // dwdot[O2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
// d()/d[AR]
dqdci = (1 - 1)*dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
J[130] -= dqdci;              // dwdot[H]/d[AR]
J[135] -= dqdci;              // dwdot[O2]/d[AR]
J[137] += dqdci;              // dwdot[HO2]/d[AR]
}
else {
dqdc[0] = dcdc_fac + k_f*sc[5]*sc[10];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[10];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[10];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[4]*sc[5];
k_f = 2150000
            * exp(0.292 * logT - (-297.904266156394) * invT);
dlnkfdT = 0.292 * invT + (-297.904266156394) * invT2;
// pressure-fall-off
k_0 = 2.28e+19 * exp(-0.977 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -0.977 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-1)*exp(-T/1e-10);
Fcent2 = 1 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[4]*sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[4] + h_RT[5]) + (h_RT[4] + h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[4]*sc[5];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
// d()/d[H2O]
dqdci = (1 - 1)*dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[57] -= dqdci;               // dwdot[O2]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
// d()/d[O2]
dqdci =  + k_f*sc[0]*sc[4];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
J[72] += dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  - k_r*sc[4];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[96] -= dqdci;               // dwdot[O2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
}
else {
dqdc[0] = dcdc_fac + k_f*sc[4]*sc[5];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[4];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[4];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 33: H + O2 (+H2) <=> HO2 (+H2)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[1]*sc[5];
k_f = 4230000
            * exp(0.646 * logT - (702.339500463647) * invT);
dlnkfdT = 0.646 * invT + (702.339500463647) * invT2;
// pressure-fall-off
k_0 = 8.28e+20 * exp(-1.72 * logT - (-1.36371716433079) * invT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.72 * invT + (-1.36371716433079) * invT2;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.353)*exp(-T/1e-10);
Fcent2 = 0.353 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[1]*sc[7];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[1] + h_RT[5]) + (h_RT[1] + h_RT[7]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[1]*sc[5];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
// d()/d[H2]
dqdci = (1 - 1)*dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[18] -= dqdci;               // dwdot[O2]/d[H2]
J[20] += dqdci;               // dwdot[HO2]/d[H2]
// d()/d[O2]
dqdci =  + k_f*sc[0]*sc[1];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
J[72] += dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  - k_r*sc[1];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[96] -= dqdci;               // dwdot[O2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
}
else {
dqdc[0] = dcdc_fac + k_f*sc[1]*sc[5];
dqdc[1] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[1];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[1];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 34: H + O2 (+HE) <=> HO2 (+HE)
// a third-body and pressure-fall-off reaction
// 3-body correction factor
alpha = mixture;
// forward
phi_f = sc[0]*sc[5]*sc[11];
k_f = 1025000
            * exp(0.604 * logT - (-121.325538125518) * invT);
dlnkfdT = 0.604 * invT + (-121.325538125518) * invT2;
// pressure-fall-off
k_0 = 2.72e+19 * exp(-1.29 * logT);
Pr = 1e-12 * alpha / k_f * k_0;
fPr = Pr / (1.0+Pr);
dlnk0dT = -1.29 * invT;
dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
dlogfPrdT = dlogPrdT / (1.0+Pr);
// Troe form
logPr = log10(Pr);
Fcent1 = (1.-0.458)*exp(-T/1e-10);
Fcent2 = 0.458 * exp(-T/1e+30);
Fcent3 = exp(-1e+30 * invT);
Fcent = Fcent1 + Fcent2 + Fcent3;
logFcent = log10(Fcent);
troe_c = -.4 - .67 * logFcent;
troe_n = .75 - 1.27 * logFcent;
troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
troePr = (troe_c + logPr) * troePr_den;
troe = 1.0 / (1.0 + troePr*troePr);
F = exp(M_LN10 * logFcent * troe);
dlogFcentdT = log10e/Fcent*( 
    -Fcent1/1e-10
    -Fcent2/1e+30
    + Fcent3*1e+30*invT2
);
dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
dlogFdn = dlogFdcn_fac * troePr;
dlogFdlogPr = dlogFdc;
dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
// reverse
phi_r = sc[7]*sc[11];
Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5] + h_RT[11]) + (h_RT[7] + h_RT[11]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
Corr = fPr * F;
q = Corr * q_nocor;
dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
// update wdot
wdot[0] -= q;// H
wdot[5] -= q;// O2
wdot[7] += q;// HO2
// for convenience
k_f *= Corr;
k_r *= Corr;
dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[5]*sc[11];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
// d()/d[O2]
dqdci =  + k_f*sc[0]*sc[11];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
J[72] += dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  - k_r*sc[11];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[96] -= dqdci;               // dwdot[O2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
// d()/d[HE]
dqdci = (1 - 1)*dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[148] -= dqdci;              // dwdot[O2]/d[HE]
J[150] += dqdci;              // dwdot[HO2]/d[HE]
}
else {
dqdc[0] = dcdc_fac + k_f*sc[5]*sc[11];
dqdc[1] = dcdc_fac;
dqdc[2] = dcdc_fac;
dqdc[3] = dcdc_fac;
dqdc[4] = dcdc_fac;
dqdc[5] = dcdc_fac + k_f*sc[0]*sc[11];
dqdc[6] = dcdc_fac;
dqdc[7] = dcdc_fac - k_r*sc[11];
dqdc[8] = dcdc_fac;
dqdc[9] = dcdc_fac;
dqdc[10] = dcdc_fac;
dqdc[11] = dcdc_fac + k_f*sc[0]*sc[5] - k_r*sc[7];
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+5] -= dqdc[k];
J[13*k+7] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[161] -= dqdT;// dwdot[O2]/dT
J[163] += dqdT;// dwdot[HO2]/dT

// reaction 0: H2 + M <=> 2 H + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.83 - 1)*sc[11];
// forward
phi_f = sc[1];
k_f = 45800000000000
            * exp(-1.4 * logT - (52535.8199100127) * invT);
dlnkfdT = -1.4 * invT + (52535.8199100127) * invT2;
// reverse
phi_r = (sc[0] * sc[0]);
Kc = refC * exp(-2.000000*g_RT[0] + g_RT[1]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[0]) - 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += 2 * q;// H
wdot[1] -= q;// H2
// for convenience
k_f *= alpha;
k_r *= alpha;
if (consP == 1) {
// d()/d[H]
dqdci =  - k_r*2.000000*sc[0];
J[0] += 2 * dqdci;            // dwdot[H]/d[H]
J[1] -= dqdci;                // dwdot[H2]/d[H]
// d()/d[H2]
dqdci = (2.5 - 1)*q_nocor + k_f;
J[13] += 2 * dqdci;           // dwdot[H]/d[H2]
J[14] -= dqdci;               // dwdot[H2]/d[H2]
// d()/d[H2O]
dqdci = (12 - 1)*q_nocor;
J[52] += 2 * dqdci;           // dwdot[H]/d[H2O]
J[53] -= dqdci;               // dwdot[H2]/d[H2O]
// d()/d[HE]
dqdci = (0.83 - 1)*q_nocor;
J[143] += 2 * dqdci;          // dwdot[H]/d[HE]
J[144] -= dqdci;              // dwdot[H2]/d[HE]
}
else {
dqdc[0] = q_nocor - k_r*2.000000*sc[0];
dqdc[1] = 2.5*q_nocor + k_f;
dqdc[2] = q_nocor;
dqdc[3] = q_nocor;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = q_nocor;
dqdc[11] = 0.83*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] += 2 * dqdc[k];
J[13*k+1] -= dqdc[k];
}
}
J[156] += 2 * dqdT;// dwdot[H]/dT
J[157] -= dqdT;// dwdot[H2]/dT

// reaction 1: 2 O + M <=> O2 + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.83 - 1)*sc[10] + ( 0.83 - 1)*sc[11];
// forward
phi_f = (sc[2] * sc[2]);
k_f = 6160
            * exp(-0.5 * logT);
dlnkfdT = -0.5 * invT;
// reverse
phi_r = sc[5];
Kc = refCinv * exp(2.000000*g_RT[2] - g_RT[5]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= 2 * q;// O
wdot[5] += q;// O2
// for convenience
k_f *= alpha;
k_r *= alpha;
if (consP == 1) {
// d()/d[H2]
dqdci = (2.5 - 1)*q_nocor;
J[15] += -2 * dqdci;          // dwdot[O]/d[H2]
J[18] += dqdci;               // dwdot[O2]/d[H2]
// d()/d[O]
dqdci =  + k_f*2.000000*sc[2];
J[28] += -2 * dqdci;          // dwdot[O]/d[O]
J[31] += dqdci;               // dwdot[O2]/d[O]
// d()/d[H2O]
dqdci = (12 - 1)*q_nocor;
J[54] += -2 * dqdci;          // dwdot[O]/d[H2O]
J[57] += dqdci;               // dwdot[O2]/d[H2O]
// d()/d[O2]
dqdci =  - k_r;
J[67] += -2 * dqdci;          // dwdot[O]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
// d()/d[AR]
dqdci = (0.83 - 1)*q_nocor;
J[132] += -2 * dqdci;         // dwdot[O]/d[AR]
J[135] += dqdci;              // dwdot[O2]/d[AR]
// d()/d[HE]
dqdci = (0.83 - 1)*q_nocor;
J[145] += -2 * dqdci;         // dwdot[O]/d[HE]
J[148] += dqdci;              // dwdot[O2]/d[HE]
}
else {
dqdc[0] = q_nocor;
dqdc[1] = 2.5*q_nocor;
dqdc[2] = q_nocor + k_f*2.000000*sc[2];
dqdc[3] = q_nocor;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor - k_r;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.83*q_nocor;
dqdc[11] = 0.83*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+2] += -2 * dqdc[k];
J[13*k+5] += dqdc[k];
}
}
J[158] += -2 * dqdT;// dwdot[O]/dT
J[161] += dqdT;// dwdot[O2]/dT

// reaction 5: H + OH + M <=> H2O + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 0.73 - 1)*sc[1] + ( 3.65 - 1)*sc[4] + ( 0.38 - 1)*sc[10] + ( 0 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[3];
k_f = 35000000000
            * exp(-2 * logT);
dlnkfdT = -2 * invT;
// reverse
phi_r = sc[4];
Kc = refCinv * exp(g_RT[0] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[3]) + (h_RT[4]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// for convenience
k_f *= alpha;
k_r *= alpha;
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[3];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] -= dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
// d()/d[H2]
dqdci = (0.73 - 1)*q_nocor;
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[16] -= dqdci;               // dwdot[OH]/d[H2]
J[17] += dqdci;               // dwdot[H2O]/d[H2]
// d()/d[OH]
dqdci =  + k_f*sc[0];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci = (3.65 - 1)*q_nocor - k_r;
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
// d()/d[AR]
dqdci = (0.38 - 1)*q_nocor;
J[130] -= dqdci;              // dwdot[H]/d[AR]
J[133] -= dqdci;              // dwdot[OH]/d[AR]
J[134] += dqdci;              // dwdot[H2O]/d[AR]
// d()/d[HE]
dqdci = (0 - 1)*q_nocor;
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[146] -= dqdci;              // dwdot[OH]/d[HE]
J[147] += dqdci;              // dwdot[H2O]/d[HE]
}
else {
dqdc[0] = q_nocor + k_f*sc[3];
dqdc[1] = 0.73*q_nocor;
dqdc[2] = q_nocor;
dqdc[3] = q_nocor + k_f*sc[0];
dqdc[4] = 3.65*q_nocor - k_r;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.38*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+3] -= dqdc[k];
J[13*k+4] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[159] -= dqdT;// dwdot[OH]/dT
J[160] += dqdT;// dwdot[H2O]/dT

// reaction 6: H + OH + HE <=> H2O + HE
// a third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[3]*sc[11];
k_f = 404000000000000
            * exp(-3.234 * logT);
dlnkfdT = -3.234 * invT;
// reverse
phi_r = sc[4]*sc[11];
Kc = refCinv * exp(g_RT[0] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[3] + h_RT[11]) + (h_RT[4] + h_RT[11]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// d()/d[H]
dqdci =  + k_f*sc[3]*sc[11];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] -= dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
// d()/d[OH]
dqdci =  + k_f*sc[0]*sc[11];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[11];
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
// d()/d[HE]
dqdci =  + k_f*sc[0]*sc[3] - k_r*sc[4];
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[146] -= dqdci;              // dwdot[OH]/d[HE]
J[147] += dqdci;              // dwdot[H2O]/d[HE]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT

// reaction 8: H + O + M <=> OH + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[10] + ( 0.75 - 1)*sc[11];
// forward
phi_f = sc[0]*sc[2];
k_f = 4710000
            * exp(-1 * logT);
dlnkfdT = -1 * invT;
// reverse
phi_r = sc[3];
Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] -= q;// O
wdot[3] += q;// OH
// for convenience
k_f *= alpha;
k_r *= alpha;
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[2];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[2] -= dqdci;                // dwdot[O]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
// d()/d[H2]
dqdci = (2.5 - 1)*q_nocor;
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[15] -= dqdci;               // dwdot[O]/d[H2]
J[16] += dqdci;               // dwdot[OH]/d[H2]
// d()/d[O]
dqdci =  + k_f*sc[0];
J[26] -= dqdci;               // dwdot[H]/d[O]
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
// d()/d[OH]
dqdci =  - k_r;
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
// d()/d[H2O]
dqdci = (12 - 1)*q_nocor;
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[54] -= dqdci;               // dwdot[O]/d[H2O]
J[55] += dqdci;               // dwdot[OH]/d[H2O]
// d()/d[AR]
dqdci = (0.75 - 1)*q_nocor;
J[130] -= dqdci;              // dwdot[H]/d[AR]
J[132] -= dqdci;              // dwdot[O]/d[AR]
J[133] += dqdci;              // dwdot[OH]/d[AR]
// d()/d[HE]
dqdci = (0.75 - 1)*q_nocor;
J[143] -= dqdci;              // dwdot[H]/d[HE]
J[145] -= dqdci;              // dwdot[O]/d[HE]
J[146] += dqdci;              // dwdot[OH]/d[HE]
}
else {
dqdc[0] = q_nocor + k_f*sc[2];
dqdc[1] = 2.5*q_nocor;
dqdc[2] = q_nocor + k_f*sc[0];
dqdc[3] = q_nocor - k_r;
dqdc[4] = 12*q_nocor;
dqdc[5] = q_nocor;
dqdc[6] = q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor;
dqdc[10] = 0.75*q_nocor;
dqdc[11] = 0.75*q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+2] -= dqdc[k];
J[13*k+3] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[158] -= dqdT;// dwdot[O]/dT
J[159] += dqdT;// dwdot[OH]/dT

// reaction 9: H + O + M <=> OHV + M
// a third-body and non-pressure-fall-off reaction
// 3-body correction factor
alpha = mixture + ( 6.5 - 1)*sc[4] + ( 0.4 - 1)*sc[5] + ( 0.4 - 1)*sc[6] + ( 0.35 - 1)*sc[10];
// forward
phi_f = sc[0]*sc[2];
k_f = 15
            * exp(- (3006.7195781832) * invT);
dlnkfdT = (3006.7195781832) * invT2;
// reverse
phi_r = sc[9];
Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[9]) + 1.000000);
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q_nocor = k_f*phi_f - k_r*phi_r;
q = alpha * q_nocor;
dqdT = alpha *(dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] -= q;// O
wdot[9] += q;// OHV
// for convenience
k_f *= alpha;
k_r *= alpha;
if (consP == 1) {
// d()/d[H]
dqdci =  + k_f*sc[2];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[2] -= dqdci;                // dwdot[O]/d[H]
J[9] += dqdci;                // dwdot[OHV]/d[H]
// d()/d[H2]
dqdci = (1 - 1)*q_nocor;
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[15] -= dqdci;               // dwdot[O]/d[H2]
J[22] += dqdci;               // dwdot[OHV]/d[H2]
// d()/d[O]
dqdci =  + k_f*sc[0];
J[26] -= dqdci;               // dwdot[H]/d[O]
J[28] -= dqdci;               // dwdot[O]/d[O]
J[35] += dqdci;               // dwdot[OHV]/d[O]
// d()/d[H2O]
dqdci = (6.5 - 1)*q_nocor;
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[54] -= dqdci;               // dwdot[O]/d[H2O]
J[61] += dqdci;               // dwdot[OHV]/d[H2O]
// d()/d[O2]
dqdci = (0.4 - 1)*q_nocor;
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[67] -= dqdci;               // dwdot[O]/d[O2]
J[74] += dqdci;               // dwdot[OHV]/d[O2]
// d()/d[N2]
dqdci = (0.4 - 1)*q_nocor;
J[78] -= dqdci;               // dwdot[H]/d[N2]
J[80] -= dqdci;               // dwdot[O]/d[N2]
J[87] += dqdci;               // dwdot[OHV]/d[N2]
// d()/d[OHV]
dqdci =  - k_r;
J[117] -= dqdci;              // dwdot[H]/d[OHV]
J[119] -= dqdci;              // dwdot[O]/d[OHV]
J[126] += dqdci;              // dwdot[OHV]/d[OHV]
// d()/d[AR]
dqdci = (0.35 - 1)*q_nocor;
J[130] -= dqdci;              // dwdot[H]/d[AR]
J[132] -= dqdci;              // dwdot[O]/d[AR]
J[139] += dqdci;              // dwdot[OHV]/d[AR]
}
else {
dqdc[0] = q_nocor + k_f*sc[2];
dqdc[1] = q_nocor;
dqdc[2] = q_nocor + k_f*sc[0];
dqdc[3] = q_nocor;
dqdc[4] = 6.5*q_nocor;
dqdc[5] = 0.4*q_nocor;
dqdc[6] = 0.4*q_nocor;
dqdc[7] = q_nocor;
dqdc[8] = q_nocor;
dqdc[9] = q_nocor - k_r;
dqdc[10] = 0.35*q_nocor;
dqdc[11] = q_nocor;
for (int k=0; k<12; k++) {
J[13*k+0] -= dqdc[k];
J[13*k+2] -= dqdc[k];
J[13*k+9] += dqdc[k];
}
}
J[156] -= dqdT;// dwdot[H]/dT
J[158] -= dqdT;// dwdot[O]/dT
J[165] += dqdT;// dwdot[OHV]/dT

// reaction 2: H2 + O <=> H + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[2];
k_f = 0.0508
            * exp(2.7 * logT - (3150.13632793755) * invT);
dlnkfdT = 2.7 * invT + (3150.13632793755) * invT2;
// reverse
phi_r = sc[0]*sc[3];
Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += q;// H
wdot[1] -= q;// H2
wdot[2] -= q;// O
wdot[3] += q;// OH
// d()/d[H]
dqdci =  - k_r*sc[3];
J[0] += dqdci;                // dwdot[H]/d[H]
J[1] -= dqdci;                // dwdot[H2]/d[H]
J[2] -= dqdci;                // dwdot[O]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
// d()/d[H2]
dqdci =  + k_f*sc[2];
J[13] += dqdci;               // dwdot[H]/d[H2]
J[14] -= dqdci;               // dwdot[H2]/d[H2]
J[15] -= dqdci;               // dwdot[O]/d[H2]
J[16] += dqdci;               // dwdot[OH]/d[H2]
// d()/d[O]
dqdci =  + k_f*sc[1];
J[26] += dqdci;               // dwdot[H]/d[O]
J[27] -= dqdci;               // dwdot[H2]/d[O]
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[0];
J[39] += dqdci;               // dwdot[H]/d[OH]
J[40] -= dqdci;               // dwdot[H2]/d[OH]
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
// d()/dT
J[156] += dqdT;               // dwdot[H]/dT
J[157] -= dqdT;               // dwdot[H2]/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT

// reaction 3: H2 + OH <=> H + H2O
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[3];
k_f = 43800000
            * exp(- (3517.48449397499) * invT);
dlnkfdT = (3517.48449397499) * invT2;
// reverse
phi_r = sc[0]*sc[4];
Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] += q;// H
wdot[1] -= q;// H2
wdot[3] -= q;// OH
wdot[4] += q;// H2O
// d()/d[H]
dqdci =  - k_r*sc[4];
J[0] += dqdci;                // dwdot[H]/d[H]
J[1] -= dqdci;                // dwdot[H2]/d[H]
J[3] -= dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
// d()/d[H2]
dqdci =  + k_f*sc[3];
J[13] += dqdci;               // dwdot[H]/d[H2]
J[14] -= dqdci;               // dwdot[H2]/d[H2]
J[16] -= dqdci;               // dwdot[OH]/d[H2]
J[17] += dqdci;               // dwdot[H2O]/d[H2]
// d()/d[OH]
dqdci =  + k_f*sc[1];
J[39] += dqdci;               // dwdot[H]/d[OH]
J[40] -= dqdci;               // dwdot[H2]/d[OH]
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[0];
J[52] += dqdci;               // dwdot[H]/d[H2O]
J[53] -= dqdci;               // dwdot[H2]/d[H2O]
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
// d()/dT
J[156] += dqdT;               // dwdot[H]/dT
J[157] -= dqdT;               // dwdot[H2]/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT

// reaction 4: H + O2 <=> O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[5];
k_f = 104000000
            * exp(- (7803.88405329959) * invT);
dlnkfdT = (7803.88405329959) * invT2;
// reverse
phi_r = sc[2]*sc[3];
Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[2] += q;// O
wdot[3] += q;// OH
wdot[5] -= q;// O2
// d()/d[H]
dqdci =  + k_f*sc[5];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[2] += dqdci;                // dwdot[O]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
J[5] -= dqdci;                // dwdot[O2]/d[H]
// d()/d[O]
dqdci =  - k_r*sc[3];
J[26] -= dqdci;               // dwdot[H]/d[O]
J[28] += dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[31] -= dqdci;               // dwdot[O2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[2];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[41] += dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[44] -= dqdci;               // dwdot[O2]/d[OH]
// d()/d[O2]
dqdci =  + k_f*sc[0];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[67] += dqdci;               // dwdot[O]/d[O2]
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[70] -= dqdci;               // dwdot[O2]/d[O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[158] += dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[161] -= dqdT;               // dwdot[O2]/dT

// reaction 7: H2O + O <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[4];
k_f = 2.97
            * exp(2.02 * logT - (6743.10332178324) * invT);
dlnkfdT = 2.02 * invT + (6743.10332178324) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[2] - 2.000000*g_RT[3] + g_RT[4]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += 2 * q;// OH
wdot[4] -= q;// H2O
// d()/d[O]
dqdci =  + k_f*sc[4];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += 2 * dqdci;           // dwdot[OH]/d[O]
J[30] -= dqdci;               // dwdot[H2O]/d[O]
// d()/d[OH]
dqdci =  - k_r*2.000000*sc[3];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += 2 * dqdci;           // dwdot[OH]/d[OH]
J[43] -= dqdci;               // dwdot[H2O]/d[OH]
// d()/d[H2O]
dqdci =  + k_f*sc[2];
J[54] -= dqdci;               // dwdot[O]/d[H2O]
J[55] += 2 * dqdci;           // dwdot[OH]/d[H2O]
J[56] -= dqdci;               // dwdot[H2O]/d[H2O]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += 2 * dqdT;           // dwdot[OH]/dT
J[160] -= dqdT;               // dwdot[H2O]/dT

// reaction 10: H2O + OHV <=> H2O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[4]*sc[9];
k_f = 5930000
            * exp(0.5 * logT - (-432.766332592059) * invT);
dlnkfdT = 0.5 * invT + (-432.766332592059) * invT2;
// reverse
phi_r = sc[3]*sc[4];
Kc = exp(-g_RT[3] + g_RT[4] - g_RT[4] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[3] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[4];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[H2O]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[55] += dqdci;               // dwdot[OH]/d[H2O]
J[61] -= dqdci;               // dwdot[OHV]/d[H2O]
// d()/d[OHV]
dqdci =  + k_f*sc[4];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 11: H2 + OHV <=> H2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[1]*sc[9];
k_f = 2950000
            * exp(0.5 * logT - (-223.428199617296) * invT);
dlnkfdT = 0.5 * invT + (-223.428199617296) * invT2;
// reverse
phi_r = sc[1]*sc[3];
Kc = exp(g_RT[1] - g_RT[1] - g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[1] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[H2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[16] += dqdci;               // dwdot[OH]/d[H2]
J[22] -= dqdci;               // dwdot[OHV]/d[H2]
// d()/d[OH]
dqdci =  - k_r*sc[1];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[1];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 12: N2 + OHV <=> N2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[6]*sc[9];
k_f = 108000
            * exp(0.5 * logT - (-624.995098929462) * invT);
dlnkfdT = 0.5 * invT + (-624.995098929462) * invT2;
// reverse
phi_r = sc[3]*sc[6];
Kc = exp(-g_RT[3] + g_RT[6] - g_RT[6] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[6]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[6];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[N2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[81] += dqdci;               // dwdot[OH]/d[N2]
J[87] -= dqdci;               // dwdot[OHV]/d[N2]
// d()/d[OHV]
dqdci =  + k_f*sc[6];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 13: OH + OHV <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[9];
k_f = 6010000
            * exp(0.5 * logT - (-384.457532674806) * invT);
dlnkfdT = 0.5 * invT + (-384.457532674806) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[3] - 2.000000*g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  + k_f*sc[9] - k_r*2.000000*sc[3];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[3];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 14: H + OHV <=> H + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[9];
k_f = 1310000
            * exp(0.5 * logT - (-84.0371831893882) * invT);
dlnkfdT = 0.5 * invT + (-84.0371831893882) * invT2;
// reverse
phi_r = sc[0]*sc[3];
Kc = exp(g_RT[0] - g_RT[0] - g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[0] + h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[H]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[3] += dqdci;                // dwdot[OH]/d[H]
J[9] -= dqdci;                // dwdot[OHV]/d[H]
// d()/d[OH]
dqdci =  - k_r*sc[0];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[0];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 15: AR + OHV <=> AR + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[9]*sc[10];
k_f = 1690000
            * exp(- (2080.80091310252) * invT);
dlnkfdT = (2080.80091310252) * invT2;
// reverse
phi_r = sc[3]*sc[10];
Kc = exp(-g_RT[3] + g_RT[9] + g_RT[10] - g_RT[10]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[9] + h_RT[10]) + (h_RT[3] + h_RT[10]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[10];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f*sc[10];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/d[AR]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[133] += dqdci;              // dwdot[OH]/d[AR]
J[139] -= dqdci;              // dwdot[OHV]/d[AR]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 16: OHV <=> OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[9];
k_f = 1450000
           ;
dlnkfdT = 0.0;
// reverse
phi_r = sc[3];
Kc = exp(-g_RT[3] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[9]) + (h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r;
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[OHV]
dqdci =  + k_f;
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 17: O2 + OHV <=> O2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[5]*sc[9];
k_f = 2100000
            * exp(0.5 * logT - (-240.537566254656) * invT);
dlnkfdT = 0.5 * invT + (-240.537566254656) * invT2;
// reverse
phi_r = sc[3]*sc[5];
Kc = exp(-g_RT[3] + g_RT[5] - g_RT[5] + g_RT[9]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[3] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] += q;// OH
wdot[9] -= q;// OHV
// d()/d[OH]
dqdci =  - k_r*sc[5];
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[48] -= dqdci;               // dwdot[OHV]/d[OH]
// d()/d[O2]
dqdci =  + k_f*sc[9] - k_r*sc[3];
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[74] -= dqdci;               // dwdot[OHV]/d[O2]
// d()/d[OHV]
dqdci =  + k_f*sc[5];
J[120] += dqdci;              // dwdot[OH]/d[OHV]
J[126] -= dqdci;              // dwdot[OHV]/d[OHV]
// d()/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[165] -= dqdT;               // dwdot[OHV]/dT

// reaction 19: H + H2O2 <=> H2O + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[8];
k_f = 24100000
            * exp(- (1997.77016324474) * invT);
dlnkfdT = (1997.77016324474) * invT2;
// reverse
phi_r = sc[3]*sc[4];
Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[3] + h_RT[4]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] += q;// OH
wdot[4] += q;// H2O
wdot[8] -= q;// H2O2
// d()/d[H]
dqdci =  + k_f*sc[8];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] += dqdci;                // dwdot[OH]/d[H]
J[4] += dqdci;                // dwdot[H2O]/d[H]
J[8] -= dqdci;                // dwdot[H2O2]/d[H]
// d()/d[OH]
dqdci =  - k_r*sc[4];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[3];
J[52] -= dqdci;               // dwdot[H]/d[H2O]
J[55] += dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[H2O2]
dqdci =  + k_f*sc[0];
J[104] -= dqdci;              // dwdot[H]/d[H2O2]
J[107] += dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 20: H + H2O2 <=> H2 + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[8];
k_f = 21500
            * exp(1 * logT - (3019.29999482832) * invT);
dlnkfdT = 1 * invT + (3019.29999482832) * invT2;
// reverse
phi_r = sc[1]*sc[7];
Kc = exp(g_RT[0] - g_RT[1] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[1] += q;// H2
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[H]
dqdci =  + k_f*sc[8];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[1] += dqdci;                // dwdot[H2]/d[H]
J[7] += dqdci;                // dwdot[HO2]/d[H]
J[8] -= dqdci;                // dwdot[H2O2]/d[H]
// d()/d[H2]
dqdci =  - k_r*sc[7];
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[14] += dqdci;               // dwdot[H2]/d[H2]
J[20] += dqdci;               // dwdot[HO2]/d[H2]
J[21] -= dqdci;               // dwdot[H2O2]/d[H2]
// d()/d[HO2]
dqdci =  - k_r*sc[1];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[92] += dqdci;               // dwdot[H2]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[0];
J[104] -= dqdci;              // dwdot[H]/d[H2O2]
J[105] += dqdci;              // dwdot[H2]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[157] += dqdT;               // dwdot[H2]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 21: H2O2 + O <=> HO2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[8];
k_f = 9.55
            * exp(2 * logT - (1997.77016324474) * invT);
dlnkfdT = 2 * invT + (1997.77016324474) * invT2;
// reverse
phi_r = sc[3]*sc[7];
Kc = exp(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[3] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += q;// OH
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[O]
dqdci =  + k_f*sc[8];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[33] += dqdci;               // dwdot[HO2]/d[O]
J[34] -= dqdci;               // dwdot[H2O2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[7];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[HO2]
dqdci =  - k_r*sc[3];
J[93] -= dqdci;               // dwdot[O]/d[HO2]
J[94] += dqdci;               // dwdot[OH]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[2];
J[106] -= dqdci;              // dwdot[O]/d[H2O2]
J[107] += dqdci;              // dwdot[OH]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 22: H2O2 + OH <=> H2O + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[8];
k_f = 1740000
            * exp(- (160.022899725901) * invT);
dlnkfdT = (160.022899725901) * invT2;
// reverse
phi_r = sc[4]*sc[7];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[OH]
dqdci =  + k_f*sc[8];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[7];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[HO2]
dqdci =  - k_r*sc[4];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[3];
J[107] -= dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 23: H2O2 + OH <=> H2O + HO2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[8];
k_f = 75900000
            * exp(- (3657.88194373451) * invT);
dlnkfdT = (3657.88194373451) * invT2;
// reverse
phi_r = sc[4]*sc[7];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[7]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[7] += q;// HO2
wdot[8] -= q;// H2O2
// d()/d[OH]
dqdci =  + k_f*sc[8];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[46] += dqdci;               // dwdot[HO2]/d[OH]
J[47] -= dqdci;               // dwdot[H2O2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[7];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[59] += dqdci;               // dwdot[HO2]/d[H2O]
J[60] -= dqdci;               // dwdot[H2O2]/d[H2O]
// d()/d[HO2]
dqdci =  - k_r*sc[4];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[98] += dqdci;               // dwdot[HO2]/d[HO2]
J[99] -= dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  + k_f*sc[3];
J[107] -= dqdci;              // dwdot[OH]/d[H2O2]
J[108] += dqdci;              // dwdot[H2O]/d[H2O2]
J[111] += dqdci;              // dwdot[HO2]/d[H2O2]
J[112] -= dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[163] += dqdT;               // dwdot[HO2]/dT
J[164] -= dqdT;               // dwdot[H2O2]/dT

// reaction 24: H + HO2 <=> 2 OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[7];
k_f = 70800000
            * exp(- (148.448916412392) * invT);
dlnkfdT = (148.448916412392) * invT2;
// reverse
phi_r = (sc[3] * sc[3]);
Kc = exp(g_RT[0] - 2.000000*g_RT[3] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (2.000000*h_RT[3]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[3] += 2 * q;// OH
wdot[7] -= q;// HO2
// d()/d[H]
dqdci =  + k_f*sc[7];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[3] += 2 * dqdci;            // dwdot[OH]/d[H]
J[7] -= dqdci;                // dwdot[HO2]/d[H]
// d()/d[OH]
dqdci =  - k_r*2.000000*sc[3];
J[39] -= dqdci;               // dwdot[H]/d[OH]
J[42] += 2 * dqdci;           // dwdot[OH]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[HO2]
dqdci =  + k_f*sc[0];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[94] += 2 * dqdci;           // dwdot[OH]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[159] += 2 * dqdT;           // dwdot[OH]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 25: H + HO2 <=> H2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[0]*sc[7];
k_f = 11400
            * exp(1.083 * logT - (278.782032855815) * invT);
dlnkfdT = 1.083 * invT + (278.782032855815) * invT2;
// reverse
phi_r = sc[1]*sc[5];
Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[0] -= q;// H
wdot[1] += q;// H2
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[H]
dqdci =  + k_f*sc[7];
J[0] -= dqdci;                // dwdot[H]/d[H]
J[1] += dqdci;                // dwdot[H2]/d[H]
J[5] += dqdci;                // dwdot[O2]/d[H]
J[7] -= dqdci;                // dwdot[HO2]/d[H]
// d()/d[H2]
dqdci =  - k_r*sc[5];
J[13] -= dqdci;               // dwdot[H]/d[H2]
J[14] += dqdci;               // dwdot[H2]/d[H2]
J[18] += dqdci;               // dwdot[O2]/d[H2]
J[20] -= dqdci;               // dwdot[HO2]/d[H2]
// d()/d[O2]
dqdci =  - k_r*sc[1];
J[65] -= dqdci;               // dwdot[H]/d[O2]
J[66] += dqdci;               // dwdot[H2]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[0];
J[91] -= dqdci;               // dwdot[H]/d[HO2]
J[92] += dqdci;               // dwdot[H2]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[156] -= dqdT;               // dwdot[H]/dT
J[157] += dqdT;               // dwdot[H2]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 26: HO2 + O <=> O2 + OH
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[2]*sc[7];
k_f = 32500000
           ;
dlnkfdT = 0.0;
// reverse
phi_r = sc[3]*sc[5];
Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[3] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[2] -= q;// O
wdot[3] += q;// OH
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[O]
dqdci =  + k_f*sc[7];
J[28] -= dqdci;               // dwdot[O]/d[O]
J[29] += dqdci;               // dwdot[OH]/d[O]
J[31] += dqdci;               // dwdot[O2]/d[O]
J[33] -= dqdci;               // dwdot[HO2]/d[O]
// d()/d[OH]
dqdci =  - k_r*sc[5];
J[41] -= dqdci;               // dwdot[O]/d[OH]
J[42] += dqdci;               // dwdot[OH]/d[OH]
J[44] += dqdci;               // dwdot[O2]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[O2]
dqdci =  - k_r*sc[3];
J[67] -= dqdci;               // dwdot[O]/d[O2]
J[68] += dqdci;               // dwdot[OH]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[2];
J[93] -= dqdci;               // dwdot[O]/d[HO2]
J[94] += dqdci;               // dwdot[OH]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[158] -= dqdT;               // dwdot[O]/dT
J[159] += dqdT;               // dwdot[OH]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 27: HO2 + OH <=> H2O + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = sc[3]*sc[7];
k_f = 31648300
            * exp(1.380154522 * logT - (15400.9555526519) * invT);
dlnkfdT = 1.380154522 * invT + (15400.9555526519) * invT2;
// reverse
phi_r = sc[4]*sc[5];
Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[5]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[3] -= q;// OH
wdot[4] += q;// H2O
wdot[5] += q;// O2
wdot[7] -= q;// HO2
// d()/d[OH]
dqdci =  + k_f*sc[7];
J[42] -= dqdci;               // dwdot[OH]/d[OH]
J[43] += dqdci;               // dwdot[H2O]/d[OH]
J[44] += dqdci;               // dwdot[O2]/d[OH]
J[46] -= dqdci;               // dwdot[HO2]/d[OH]
// d()/d[H2O]
dqdci =  - k_r*sc[5];
J[55] -= dqdci;               // dwdot[OH]/d[H2O]
J[56] += dqdci;               // dwdot[H2O]/d[H2O]
J[57] += dqdci;               // dwdot[O2]/d[H2O]
J[59] -= dqdci;               // dwdot[HO2]/d[H2O]
// d()/d[O2]
dqdci =  - k_r*sc[4];
J[68] -= dqdci;               // dwdot[OH]/d[O2]
J[69] += dqdci;               // dwdot[H2O]/d[O2]
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] -= dqdci;               // dwdot[HO2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*sc[3];
J[94] -= dqdci;               // dwdot[OH]/d[HO2]
J[95] += dqdci;               // dwdot[H2O]/d[HO2]
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] -= dqdci;               // dwdot[HO2]/d[HO2]
// d()/dT
J[159] -= dqdT;               // dwdot[OH]/dT
J[160] += dqdT;               // dwdot[H2O]/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] -= dqdT;               // dwdot[HO2]/dT

// reaction 28: 2 HO2 <=> H2O2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = (sc[7] * sc[7]);
k_f = 12100
            * exp(0.422 * logT - (-745.26388205679) * invT);
dlnkfdT = 0.422 * invT + (-745.26388205679) * invT2;
// reverse
phi_r = sc[5]*sc[8];
Kc = exp(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[7]) + (h_RT[5] + h_RT[8]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[5] += q;// O2
wdot[7] -= 2 * q;// HO2
wdot[8] += q;// H2O2
// d()/d[O2]
dqdci =  - k_r*sc[8];
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] += -2 * dqdci;          // dwdot[HO2]/d[O2]
J[73] += dqdci;               // dwdot[H2O2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*2.000000*sc[7];
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
J[99] += dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  - k_r*sc[5];
J[109] += dqdci;              // dwdot[O2]/d[H2O2]
J[111] += -2 * dqdci;         // dwdot[HO2]/d[H2O2]
J[112] += dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] += -2 * dqdT;          // dwdot[HO2]/dT
J[164] += dqdT;               // dwdot[H2O2]/dT

// reaction 29: 2 HO2 <=> H2O2 + O2
// a non-third-body and non-pressure-fall-off reaction
// forward
phi_f = (sc[7] * sc[7]);
k_f = 16900000000
            * exp(-0.681 * logT - (6507.59792218663) * invT);
dlnkfdT = -0.681 * invT + (6507.59792218663) * invT2;
// reverse
phi_r = sc[5]*sc[8];
Kc = exp(-g_RT[5] + 2.000000*g_RT[7] - g_RT[8]);
k_r = k_f / Kc;
dlnKcdT = invT * (-(2.000000*h_RT[7]) + (h_RT[5] + h_RT[8]));
dkrdT = (dlnkfdT - dlnKcdT)*k_r;
// rate of progress
q = k_f*phi_f - k_r*phi_r;
dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
// update wdot
wdot[5] += q;// O2
wdot[7] -= 2 * q;// HO2
wdot[8] += q;// H2O2
// d()/d[O2]
dqdci =  - k_r*sc[8];
J[70] += dqdci;               // dwdot[O2]/d[O2]
J[72] += -2 * dqdci;          // dwdot[HO2]/d[O2]
J[73] += dqdci;               // dwdot[H2O2]/d[O2]
// d()/d[HO2]
dqdci =  + k_f*2.000000*sc[7];
J[96] += dqdci;               // dwdot[O2]/d[HO2]
J[98] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
J[99] += dqdci;               // dwdot[H2O2]/d[HO2]
// d()/d[H2O2]
dqdci =  - k_r*sc[5];
J[109] += dqdci;              // dwdot[O2]/d[H2O2]
J[111] += -2 * dqdci;         // dwdot[HO2]/d[H2O2]
J[112] += dqdci;              // dwdot[H2O2]/d[H2O2]
// d()/dT
J[161] += dqdT;               // dwdot[O2]/dT
J[163] += -2 * dqdT;          // dwdot[HO2]/dT
J[164] += dqdT;               // dwdot[H2O2]/dT

amrex::Real c_R[12], dcRdT[12], e_RT[12];
amrex::Real * eh_RT;
if (consP == 1) {
cp_R(c_R, T);
dcvpRdT(dcRdT, T);
eh_RT = &h_RT[0];
}
else {
cv_R(c_R, T);
dcvpRdT(dcRdT, T);
speciesInternalEnergy(e_RT, T);
eh_RT = &e_RT[0];
}

amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
for (int k = 0; k < 12; ++k) {
cmix += c_R[k]*sc[k];
dcmixdT += dcRdT[k]*sc[k];
ehmix += eh_RT[k]*wdot[k];
dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[156+k];
}

amrex::Real cmixinv = 1.0/cmix;
amrex::Real tmp1 = ehmix*cmixinv;
amrex::Real tmp3 = cmixinv*T;
amrex::Real tmp2 = tmp1*tmp3;
amrex::Real dehmixdc;
// dTdot/d[X]
for (int k = 0; k < 12; ++k) {
dehmixdc = 0.0;
for (int m = 0; m < 12; ++m) {
dehmixdc += eh_RT[m]*J[k*13+m];
}
J[k*13+12] = tmp2*c_R[k] - tmp3*dehmixdc;
}
// dTdot/dT
J[168] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
#else
amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, const amrex::Real *  sc, const amrex::Real *  Tp, const int * consP)
{
amrex::Real c[12];

for (int k=0; k<12; k++) {
c[k] = 1.e6 * sc[k];
}

aJacobian(J, c, *Tp, *consP);

// dwdot[k]/dT
// dTdot/d[X]
for (int k=0; k<12; k++) {
J[156+k] *= 1.e-6;
J[k*13+12] *= 1.e6;
}

}


// Transport function declarations 

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
*LENIMC = 51;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
*LENRMC = 3204;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
*NO = 4;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
*KK = 12;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
*NLITE = 3;}


// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
*PATM =   0.1013250000000000E+07;}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
WT[0] = 1.00800000E+00;
WT[1] = 2.01600000E+00;
WT[2] = 1.59990000E+01;
WT[3] = 1.70070000E+01;
WT[4] = 1.80150000E+01;
WT[5] = 3.19980000E+01;
WT[6] = 2.80140000E+01;
WT[7] = 3.30060000E+01;
WT[8] = 3.40140000E+01;
WT[9] = 1.70070000E+01;
WT[10] = 3.99500000E+01;
WT[11] = 4.00260200E+00;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
EPS[0] = 1.45000000E+02;
EPS[1] = 3.80000000E+01;
EPS[2] = 8.00000000E+01;
EPS[3] = 8.00000000E+01;
EPS[4] = 5.72400000E+02;
EPS[5] = 1.07400000E+02;
EPS[6] = 9.75300000E+01;
EPS[7] = 1.07400000E+02;
EPS[8] = 1.07400000E+02;
EPS[9] = 8.00000000E+01;
EPS[10] = 1.36500000E+02;
EPS[11] = 1.02000000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
SIG[0] = 2.05000000E+00;
SIG[1] = 2.92000000E+00;
SIG[2] = 2.75000000E+00;
SIG[3] = 2.75000000E+00;
SIG[4] = 2.60500000E+00;
SIG[5] = 3.45800000E+00;
SIG[6] = 3.62100000E+00;
SIG[7] = 3.45800000E+00;
SIG[8] = 3.45800000E+00;
SIG[9] = 2.75000000E+00;
SIG[10] = 3.33000000E+00;
SIG[11] = 2.57600000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
DIP[0] = 0.00000000E+00;
DIP[1] = 0.00000000E+00;
DIP[2] = 0.00000000E+00;
DIP[3] = 0.00000000E+00;
DIP[4] = 1.84400000E+00;
DIP[5] = 0.00000000E+00;
DIP[6] = 0.00000000E+00;
DIP[7] = 0.00000000E+00;
DIP[8] = 0.00000000E+00;
DIP[9] = 0.00000000E+00;
DIP[10] = 0.00000000E+00;
DIP[11] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
POL[0] = 0.00000000E+00;
POL[1] = 7.90000000E-01;
POL[2] = 0.00000000E+00;
POL[3] = 0.00000000E+00;
POL[4] = 0.00000000E+00;
POL[5] = 1.60000000E+00;
POL[6] = 1.76000000E+00;
POL[7] = 0.00000000E+00;
POL[8] = 0.00000000E+00;
POL[9] = 0.00000000E+00;
POL[10] = 0.00000000E+00;
POL[11] = 0.00000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
ZROT[0] = 0.00000000E+00;
ZROT[1] = 2.80000000E+02;
ZROT[2] = 0.00000000E+00;
ZROT[3] = 0.00000000E+00;
ZROT[4] = 4.00000000E+00;
ZROT[5] = 3.80000000E+00;
ZROT[6] = 4.00000000E+00;
ZROT[7] = 1.00000000E+00;
ZROT[8] = 3.80000000E+00;
ZROT[9] = 0.00000000E+00;
ZROT[10] = 0.00000000E+00;
ZROT[11] = 0.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
NLIN[0] = 0;
NLIN[1] = 1;
NLIN[2] = 0;
NLIN[3] = 1;
NLIN[4] = 2;
NLIN[5] = 1;
NLIN[6] = 1;
NLIN[7] = 2;
NLIN[8] = 2;
NLIN[9] = 1;
NLIN[10] = 0;
NLIN[11] = 0;
}


// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
COFETA[0] = -1.89413239E+01;
COFETA[1] = 3.00307192E+00;
COFETA[2] = -3.02604779E-01;
COFETA[3] = 1.29293647E-02;
COFETA[4] = -1.35655309E+01;
COFETA[5] = 8.80923449E-01;
COFETA[6] = -3.20611225E-02;
COFETA[7] = 1.44621531E-03;
COFETA[8] = -1.40007723E+01;
COFETA[9] = 1.44311220E+00;
COFETA[10] = -1.02767741E-01;
COFETA[11] = 4.39662983E-03;
COFETA[12] = -1.39702229E+01;
COFETA[13] = 1.44311220E+00;
COFETA[14] = -1.02767741E-01;
COFETA[15] = 4.39662983E-03;
COFETA[16] = -1.48387839E+01;
COFETA[17] = 5.22429240E-01;
COFETA[18] = 1.42283123E-01;
COFETA[19] = -1.05730853E-02;
COFETA[20] = -1.60066453E+01;
COFETA[21] = 2.16753735E+00;
COFETA[22] = -1.97226850E-01;
COFETA[23] = 8.50065468E-03;
COFETA[24] = -1.55270224E+01;
COFETA[25] = 1.92766908E+00;
COFETA[26] = -1.66518287E-01;
COFETA[27] = 7.19100649E-03;
COFETA[28] = -1.59911374E+01;
COFETA[29] = 2.16753735E+00;
COFETA[30] = -1.97226850E-01;
COFETA[31] = 8.50065468E-03;
COFETA[32] = -1.59760959E+01;
COFETA[33] = 2.16753735E+00;
COFETA[34] = -1.97226850E-01;
COFETA[35] = 8.50065468E-03;
COFETA[36] = -1.39702229E+01;
COFETA[37] = 1.44311220E+00;
COFETA[38] = -1.02767741E-01;
COFETA[39] = 4.39662983E-03;
COFETA[40] = -1.76496901E+01;
COFETA[41] = 2.85078722E+00;
COFETA[42] = -2.83910962E-01;
COFETA[43] = 1.21649252E-02;
COFETA[44] = -1.39753223E+01;
COFETA[45] = 1.44641326E+00;
COFETA[46] = -1.20775609E-01;
COFETA[47] = 6.06541463E-03;
}


// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
COFLAM[0] = 6.08556807E-01;
COFLAM[1] = 3.00307192E+00;
COFLAM[2] = -3.02604779E-01;
COFLAM[3] = 1.29293647E-02;
COFLAM[4] = 1.31196501E+01;
COFLAM[5] = -2.21802328E+00;
COFLAM[6] = 3.76732980E-01;
COFLAM[7] = -1.54805563E-02;
COFLAM[8] = 2.78455034E+00;
COFLAM[9] = 1.44311220E+00;
COFLAM[10] = -1.02767741E-01;
COFLAM[11] = 4.39662983E-03;
COFLAM[12] = 1.48281611E+01;
COFLAM[13] = -3.56586994E+00;
COFLAM[14] = 5.84572541E-01;
COFLAM[15] = -2.59352423E-02;
COFLAM[16] = 1.72715245E+01;
COFLAM[17] = -6.29249329E+00;
COFLAM[18] = 1.14065197E+00;
COFLAM[19] = -5.72398931E-02;
COFLAM[20] = -2.22893370E+00;
COFLAM[21] = 3.04473937E+00;
COFLAM[22] = -2.97176396E-01;
COFLAM[23] = 1.29718512E-02;
COFLAM[24] = 8.03267097E+00;
COFLAM[25] = -1.36302697E+00;
COFLAM[26] = 3.28530730E-01;
COFLAM[27] = -1.65680815E-02;
COFLAM[28] = 5.13326062E+00;
COFLAM[29] = -4.21994736E-01;
COFLAM[30] = 2.40199676E-01;
COFLAM[31] = -1.37156124E-02;
COFLAM[32] = -5.97955209E-01;
COFLAM[33] = 1.97435083E+00;
COFLAM[34] = -7.77356813E-02;
COFLAM[35] = 2.98921365E-04;
COFLAM[36] = 1.65431190E+01;
COFLAM[37] = -4.32534876E+00;
COFLAM[38] = 6.96441548E-01;
COFLAM[39] = -3.14178907E-02;
COFLAM[40] = -1.77946992E+00;
COFLAM[41] = 2.85078722E+00;
COFLAM[42] = -2.83910962E-01;
COFLAM[43] = 1.21649252E-02;
COFLAM[44] = 4.19558188E+00;
COFLAM[45] = 1.44641326E+00;
COFLAM[46] = -1.20775609E-01;
COFLAM[47] = 6.06541463E-03;
}


// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
COFD[0] = -1.37059892E+01;
COFD[1] = 3.74572091E+00;
COFD[2] = -2.64900591E-01;
COFD[3] = 1.11831118E-02;
COFD[4] = -1.05322296E+01;
COFD[5] = 2.38287287E+00;
COFD[6] = -9.25451677E-02;
COFD[7] = 3.90477283E-03;
COFD[8] = -1.25811019E+01;
COFD[9] = 3.11263114E+00;
COFD[10] = -1.86688345E-01;
COFD[11] = 7.95911073E-03;
COFD[12] = -1.25828615E+01;
COFD[13] = 3.11263114E+00;
COFD[14] = -1.86688345E-01;
COFD[15] = 7.95911073E-03;
COFD[16] = -1.88399011E+01;
COFD[17] = 5.28766000E+00;
COFD[18] = -4.43001094E-01;
COFD[19] = 1.80463788E-02;
COFD[20] = -1.37328308E+01;
COFD[21] = 3.42674072E+00;
COFD[22] = -2.25944624E-01;
COFD[23] = 9.59602387E-03;
COFD[24] = -1.34540667E+01;
COFD[25] = 3.30341833E+00;
COFD[26] = -2.10450242E-01;
COFD[27] = 8.94670534E-03;
COFD[28] = -1.37332973E+01;
COFD[29] = 3.42674072E+00;
COFD[30] = -2.25944624E-01;
COFD[31] = 9.59602387E-03;
COFD[32] = -1.37337367E+01;
COFD[33] = 3.42674072E+00;
COFD[34] = -2.25944624E-01;
COFD[35] = 9.59602387E-03;
COFD[36] = -1.25828615E+01;
COFD[37] = 3.11263114E+00;
COFD[38] = -1.86688345E-01;
COFD[39] = 7.95911073E-03;
COFD[40] = -1.44320427E+01;
COFD[41] = 3.69390556E+00;
COFD[42] = -2.58816176E-01;
COFD[43] = 1.09455115E-02;
COFD[44] = -9.01870427E+00;
COFD[45] = 1.86885237E+00;
COFD[46] = -2.79898003E-02;
COFD[47] = 1.21327131E-03;
COFD[48] = -1.05322296E+01;
COFD[49] = 2.38287287E+00;
COFD[50] = -9.25451677E-02;
COFD[51] = 3.90477283E-03;
COFD[52] = -9.61324263E+00;
COFD[53] = 1.87738363E+00;
COFD[54] = -2.92985344E-02;
COFD[55] = 1.27857667E-03;
COFD[56] = -1.09715308E+01;
COFD[57] = 2.31040527E+00;
COFD[58] = -8.74442763E-02;
COFD[59] = 3.87525190E-03;
COFD[60] = -1.09748582E+01;
COFD[61] = 2.31040527E+00;
COFD[62] = -8.74442763E-02;
COFD[63] = 3.87525190E-03;
COFD[64] = -1.60490154E+01;
COFD[65] = 4.12894506E+00;
COFD[66] = -3.11280347E-01;
COFD[67] = 1.30569870E-02;
COFD[68] = -1.15530097E+01;
COFD[69] = 2.41720750E+00;
COFD[70] = -1.00182262E-01;
COFD[71] = 4.37376763E-03;
COFD[72] = -1.15279181E+01;
COFD[73] = 2.39753083E+00;
COFD[74] = -9.83466917E-02;
COFD[75] = 4.32694147E-03;
COFD[76] = -1.15539155E+01;
COFD[77] = 2.41720750E+00;
COFD[78] = -1.00182262E-01;
COFD[79] = 4.37376763E-03;
COFD[80] = -1.15547692E+01;
COFD[81] = 2.41720750E+00;
COFD[82] = -1.00182262E-01;
COFD[83] = 4.37376763E-03;
COFD[84] = -1.09748582E+01;
COFD[85] = 2.31040527E+00;
COFD[86] = -8.74442763E-02;
COFD[87] = 3.87525190E-03;
COFD[88] = -1.14507066E+01;
COFD[89] = 2.36127381E+00;
COFD[90] = -8.99526135E-02;
COFD[91] = 3.80113640E-03;
COFD[92] = -9.64985378E+00;
COFD[93] = 1.96102826E+00;
COFD[94] = -4.58416954E-02;
COFD[95] = 2.27372333E-03;
COFD[96] = -1.25811019E+01;
COFD[97] = 3.11263114E+00;
COFD[98] = -1.86688345E-01;
COFD[99] = 7.95911073E-03;
COFD[100] = -1.09715308E+01;
COFD[101] = 2.31040527E+00;
COFD[102] = -8.74442763E-02;
COFD[103] = 3.87525190E-03;
COFD[104] = -1.22676388E+01;
COFD[105] = 2.49197560E+00;
COFD[106] = -1.06603364E-01;
COFD[107] = 4.50971238E-03;
COFD[108] = -1.22826802E+01;
COFD[109] = 2.49197560E+00;
COFD[110] = -1.06603364E-01;
COFD[111] = 4.50971238E-03;
COFD[112] = -1.82875611E+01;
COFD[113] = 4.66435737E+00;
COFD[114] = -3.74296088E-01;
COFD[115] = 1.55332609E-02;
COFD[116] = -1.34772511E+01;
COFD[117] = 2.80526779E+00;
COFD[118] = -1.47524446E-01;
COFD[119] = 6.29376912E-03;
COFD[120] = -1.32254582E+01;
COFD[121] = 2.69833584E+00;
COFD[122] = -1.33541002E-01;
COFD[123] = 5.68346517E-03;
COFD[124] = -1.34823671E+01;
COFD[125] = 2.80526779E+00;
COFD[126] = -1.47524446E-01;
COFD[127] = 6.29376912E-03;
COFD[128] = -1.34872282E+01;
COFD[129] = 2.80526779E+00;
COFD[130] = -1.47524446E-01;
COFD[131] = 6.29376912E-03;
COFD[132] = -1.22826802E+01;
COFD[133] = 2.49197560E+00;
COFD[134] = -1.06603364E-01;
COFD[135] = 4.50971238E-03;
COFD[136] = -1.41177570E+01;
COFD[137] = 3.04580143E+00;
COFD[138] = -1.78216259E-01;
COFD[139] = 7.60087614E-03;
COFD[140] = -1.02432185E+01;
COFD[141] = 2.01382318E+00;
COFD[142] = -5.07824773E-02;
COFD[143] = 2.38045050E-03;
COFD[144] = -1.25828615E+01;
COFD[145] = 3.11263114E+00;
COFD[146] = -1.86688345E-01;
COFD[147] = 7.95911073E-03;
COFD[148] = -1.09748582E+01;
COFD[149] = 2.31040527E+00;
COFD[150] = -8.74442763E-02;
COFD[151] = 3.87525190E-03;
COFD[152] = -1.22826802E+01;
COFD[153] = 2.49197560E+00;
COFD[154] = -1.06603364E-01;
COFD[155] = 4.50971238E-03;
COFD[156] = -1.22981882E+01;
COFD[157] = 2.49197560E+00;
COFD[158] = -1.06603364E-01;
COFD[159] = 4.50971238E-03;
COFD[160] = -1.83035084E+01;
COFD[161] = 4.66435737E+00;
COFD[162] = -3.74296088E-01;
COFD[163] = 1.55332609E-02;
COFD[164] = -1.34974086E+01;
COFD[165] = 2.80526779E+00;
COFD[166] = -1.47524446E-01;
COFD[167] = 6.29376912E-03;
COFD[168] = -1.32446856E+01;
COFD[169] = 2.69833584E+00;
COFD[170] = -1.33541002E-01;
COFD[171] = 5.68346517E-03;
COFD[172] = -1.35027362E+01;
COFD[173] = 2.80526779E+00;
COFD[174] = -1.47524446E-01;
COFD[175] = 6.29376912E-03;
COFD[176] = -1.35078005E+01;
COFD[177] = 2.80526779E+00;
COFD[178] = -1.47524446E-01;
COFD[179] = 6.29376912E-03;
COFD[180] = -1.22981882E+01;
COFD[181] = 2.49197560E+00;
COFD[182] = -1.06603364E-01;
COFD[183] = 4.50971238E-03;
COFD[184] = -1.41393784E+01;
COFD[185] = 3.04580143E+00;
COFD[186] = -1.78216259E-01;
COFD[187] = 7.60087614E-03;
COFD[188] = -1.02491843E+01;
COFD[189] = 2.01382318E+00;
COFD[190] = -5.07824773E-02;
COFD[191] = 2.38045050E-03;
COFD[192] = -1.88399011E+01;
COFD[193] = 5.28766000E+00;
COFD[194] = -4.43001094E-01;
COFD[195] = 1.80463788E-02;
COFD[196] = -1.60490154E+01;
COFD[197] = 4.12894506E+00;
COFD[198] = -3.11280347E-01;
COFD[199] = 1.30569870E-02;
COFD[200] = -1.82875611E+01;
COFD[201] = 4.66435737E+00;
COFD[202] = -3.74296088E-01;
COFD[203] = 1.55332609E-02;
COFD[204] = -1.83035084E+01;
COFD[205] = 4.66435737E+00;
COFD[206] = -3.74296088E-01;
COFD[207] = 1.55332609E-02;
COFD[208] = -1.58086165E+01;
COFD[209] = 2.65423839E+00;
COFD[210] = -1.10240041E-02;
COFD[211] = -3.68459754E-03;
COFD[212] = -2.04397518E+01;
COFD[213] = 5.23624914E+00;
COFD[214] = -4.38651690E-01;
COFD[215] = 1.79454254E-02;
COFD[216] = -2.01962837E+01;
COFD[217] = 5.15333330E+00;
COFD[218] = -4.30146511E-01;
COFD[219] = 1.76625034E-02;
COFD[220] = -1.97212886E+01;
COFD[221] = 5.00923653E+00;
COFD[222] = -4.13897836E-01;
COFD[223] = 1.70524216E-02;
COFD[224] = -1.97265481E+01;
COFD[225] = 5.00923653E+00;
COFD[226] = -4.13897836E-01;
COFD[227] = 1.70524216E-02;
COFD[228] = -1.83035084E+01;
COFD[229] = 4.66435737E+00;
COFD[230] = -3.74296088E-01;
COFD[231] = 1.55332609E-02;
COFD[232] = -2.04725537E+01;
COFD[233] = 5.25060921E+00;
COFD[234] = -4.39775466E-01;
COFD[235] = 1.79662640E-02;
COFD[236] = -1.15042109E+01;
COFD[237] = 2.41805816E+00;
COFD[238] = -9.70163830E-02;
COFD[239] = 4.09454835E-03;
COFD[240] = -1.37328308E+01;
COFD[241] = 3.42674072E+00;
COFD[242] = -2.25944624E-01;
COFD[243] = 9.59602387E-03;
COFD[244] = -1.15530097E+01;
COFD[245] = 2.41720750E+00;
COFD[246] = -1.00182262E-01;
COFD[247] = 4.37376763E-03;
COFD[248] = -1.34772511E+01;
COFD[249] = 2.80526779E+00;
COFD[250] = -1.47524446E-01;
COFD[251] = 6.29376912E-03;
COFD[252] = -1.34974086E+01;
COFD[253] = 2.80526779E+00;
COFD[254] = -1.47524446E-01;
COFD[255] = 6.29376912E-03;
COFD[256] = -2.04397518E+01;
COFD[257] = 5.23624914E+00;
COFD[258] = -4.38651690E-01;
COFD[259] = 1.79454254E-02;
COFD[260] = -1.47079535E+01;
COFD[261] = 3.10657376E+00;
COFD[262] = -1.85922460E-01;
COFD[263] = 7.92680827E-03;
COFD[264] = -1.44285962E+01;
COFD[265] = 2.99858376E+00;
COFD[266] = -1.72232643E-01;
COFD[267] = 7.34804765E-03;
COFD[268] = -1.47156473E+01;
COFD[269] = 3.10657376E+00;
COFD[270] = -1.85922460E-01;
COFD[271] = 7.92680827E-03;
COFD[272] = -1.47229949E+01;
COFD[273] = 3.10657376E+00;
COFD[274] = -1.85922460E-01;
COFD[275] = 7.92680827E-03;
COFD[276] = -1.34974086E+01;
COFD[277] = 2.80526779E+00;
COFD[278] = -1.47524446E-01;
COFD[279] = 6.29376912E-03;
COFD[280] = -1.53783351E+01;
COFD[281] = 3.34436607E+00;
COFD[282] = -2.15557310E-01;
COFD[283] = 9.15916913E-03;
COFD[284] = -1.03186358E+01;
COFD[285] = 1.89854108E+00;
COFD[286] = -3.33173364E-02;
COFD[287] = 1.51065205E-03;
COFD[288] = -1.34540667E+01;
COFD[289] = 3.30341833E+00;
COFD[290] = -2.10450242E-01;
COFD[291] = 8.94670534E-03;
COFD[292] = -1.15279181E+01;
COFD[293] = 2.39753083E+00;
COFD[294] = -9.83466917E-02;
COFD[295] = 4.32694147E-03;
COFD[296] = -1.32254582E+01;
COFD[297] = 2.69833584E+00;
COFD[298] = -1.33541002E-01;
COFD[299] = 5.68346517E-03;
COFD[300] = -1.32446856E+01;
COFD[301] = 2.69833584E+00;
COFD[302] = -1.33541002E-01;
COFD[303] = 5.68346517E-03;
COFD[304] = -2.01962837E+01;
COFD[305] = 5.15333330E+00;
COFD[306] = -4.30146511E-01;
COFD[307] = 1.76625034E-02;
COFD[308] = -1.44285962E+01;
COFD[309] = 2.99858376E+00;
COFD[310] = -1.72232643E-01;
COFD[311] = 7.34804765E-03;
COFD[312] = -1.42056777E+01;
COFD[313] = 2.91297621E+00;
COFD[314] = -1.61544771E-01;
COFD[315] = 6.90271324E-03;
COFD[316] = -1.44357756E+01;
COFD[317] = 2.99858376E+00;
COFD[318] = -1.72232643E-01;
COFD[319] = 7.34804765E-03;
COFD[320] = -1.44426249E+01;
COFD[321] = 2.99858376E+00;
COFD[322] = -1.72232643E-01;
COFD[323] = 7.34804765E-03;
COFD[324] = -1.32446856E+01;
COFD[325] = 2.69833584E+00;
COFD[326] = -1.33541002E-01;
COFD[327] = 5.68346517E-03;
COFD[328] = -1.51258950E+01;
COFD[329] = 3.25013576E+00;
COFD[330] = -2.03913392E-01;
COFD[331] = 8.67915586E-03;
COFD[332] = -1.04392040E+01;
COFD[333] = 1.93661815E+00;
COFD[334] = -3.90724483E-02;
COFD[335] = 1.79631123E-03;
COFD[336] = -1.37332973E+01;
COFD[337] = 3.42674072E+00;
COFD[338] = -2.25944624E-01;
COFD[339] = 9.59602387E-03;
COFD[340] = -1.15539155E+01;
COFD[341] = 2.41720750E+00;
COFD[342] = -1.00182262E-01;
COFD[343] = 4.37376763E-03;
COFD[344] = -1.34823671E+01;
COFD[345] = 2.80526779E+00;
COFD[346] = -1.47524446E-01;
COFD[347] = 6.29376912E-03;
COFD[348] = -1.35027362E+01;
COFD[349] = 2.80526779E+00;
COFD[350] = -1.47524446E-01;
COFD[351] = 6.29376912E-03;
COFD[352] = -1.97212886E+01;
COFD[353] = 5.00923653E+00;
COFD[354] = -4.13897836E-01;
COFD[355] = 1.70524216E-02;
COFD[356] = -1.47156473E+01;
COFD[357] = 3.10657376E+00;
COFD[358] = -1.85922460E-01;
COFD[359] = 7.92680827E-03;
COFD[360] = -1.44357756E+01;
COFD[361] = 2.99858376E+00;
COFD[362] = -1.72232643E-01;
COFD[363] = 7.34804765E-03;
COFD[364] = -1.47234615E+01;
COFD[365] = 3.10657376E+00;
COFD[366] = -1.85922460E-01;
COFD[367] = 7.92680827E-03;
COFD[368] = -1.47309256E+01;
COFD[369] = 3.10657376E+00;
COFD[370] = -1.85922460E-01;
COFD[371] = 7.92680827E-03;
COFD[372] = -1.35027362E+01;
COFD[373] = 2.80526779E+00;
COFD[374] = -1.47524446E-01;
COFD[375] = 6.29376912E-03;
COFD[376] = -1.53868866E+01;
COFD[377] = 3.34436607E+00;
COFD[378] = -2.15557310E-01;
COFD[379] = 9.15916913E-03;
COFD[380] = -1.03203364E+01;
COFD[381] = 1.89854108E+00;
COFD[382] = -3.33173364E-02;
COFD[383] = 1.51065205E-03;
COFD[384] = -1.37337367E+01;
COFD[385] = 3.42674072E+00;
COFD[386] = -2.25944624E-01;
COFD[387] = 9.59602387E-03;
COFD[388] = -1.15547692E+01;
COFD[389] = 2.41720750E+00;
COFD[390] = -1.00182262E-01;
COFD[391] = 4.37376763E-03;
COFD[392] = -1.34872282E+01;
COFD[393] = 2.80526779E+00;
COFD[394] = -1.47524446E-01;
COFD[395] = 6.29376912E-03;
COFD[396] = -1.35078005E+01;
COFD[397] = 2.80526779E+00;
COFD[398] = -1.47524446E-01;
COFD[399] = 6.29376912E-03;
COFD[400] = -1.97265481E+01;
COFD[401] = 5.00923653E+00;
COFD[402] = -4.13897836E-01;
COFD[403] = 1.70524216E-02;
COFD[404] = -1.47229949E+01;
COFD[405] = 3.10657376E+00;
COFD[406] = -1.85922460E-01;
COFD[407] = 7.92680827E-03;
COFD[408] = -1.44426249E+01;
COFD[409] = 2.99858376E+00;
COFD[410] = -1.72232643E-01;
COFD[411] = 7.34804765E-03;
COFD[412] = -1.47309256E+01;
COFD[413] = 3.10657376E+00;
COFD[414] = -1.85922460E-01;
COFD[415] = 7.92680827E-03;
COFD[416] = -1.47385029E+01;
COFD[417] = 3.10657376E+00;
COFD[418] = -1.85922460E-01;
COFD[419] = 7.92680827E-03;
COFD[420] = -1.35078005E+01;
COFD[421] = 2.80526779E+00;
COFD[422] = -1.47524446E-01;
COFD[423] = 6.29376912E-03;
COFD[424] = -1.53950670E+01;
COFD[425] = 3.34436607E+00;
COFD[426] = -2.15557310E-01;
COFD[427] = 9.15916913E-03;
COFD[428] = -1.03219415E+01;
COFD[429] = 1.89854108E+00;
COFD[430] = -3.33173364E-02;
COFD[431] = 1.51065205E-03;
COFD[432] = -1.25828615E+01;
COFD[433] = 3.11263114E+00;
COFD[434] = -1.86688345E-01;
COFD[435] = 7.95911073E-03;
COFD[436] = -1.09748582E+01;
COFD[437] = 2.31040527E+00;
COFD[438] = -8.74442763E-02;
COFD[439] = 3.87525190E-03;
COFD[440] = -1.22826802E+01;
COFD[441] = 2.49197560E+00;
COFD[442] = -1.06603364E-01;
COFD[443] = 4.50971238E-03;
COFD[444] = -1.22981882E+01;
COFD[445] = 2.49197560E+00;
COFD[446] = -1.06603364E-01;
COFD[447] = 4.50971238E-03;
COFD[448] = -1.83035084E+01;
COFD[449] = 4.66435737E+00;
COFD[450] = -3.74296088E-01;
COFD[451] = 1.55332609E-02;
COFD[452] = -1.34974086E+01;
COFD[453] = 2.80526779E+00;
COFD[454] = -1.47524446E-01;
COFD[455] = 6.29376912E-03;
COFD[456] = -1.32446856E+01;
COFD[457] = 2.69833584E+00;
COFD[458] = -1.33541002E-01;
COFD[459] = 5.68346517E-03;
COFD[460] = -1.35027362E+01;
COFD[461] = 2.80526779E+00;
COFD[462] = -1.47524446E-01;
COFD[463] = 6.29376912E-03;
COFD[464] = -1.35078005E+01;
COFD[465] = 2.80526779E+00;
COFD[466] = -1.47524446E-01;
COFD[467] = 6.29376912E-03;
COFD[468] = -1.22981882E+01;
COFD[469] = 2.49197560E+00;
COFD[470] = -1.06603364E-01;
COFD[471] = 4.50971238E-03;
COFD[472] = -1.41393784E+01;
COFD[473] = 3.04580143E+00;
COFD[474] = -1.78216259E-01;
COFD[475] = 7.60087614E-03;
COFD[476] = -1.02491843E+01;
COFD[477] = 2.01382318E+00;
COFD[478] = -5.07824773E-02;
COFD[479] = 2.38045050E-03;
COFD[480] = -1.44320427E+01;
COFD[481] = 3.69390556E+00;
COFD[482] = -2.58816176E-01;
COFD[483] = 1.09455115E-02;
COFD[484] = -1.14507066E+01;
COFD[485] = 2.36127381E+00;
COFD[486] = -8.99526135E-02;
COFD[487] = 3.80113640E-03;
COFD[488] = -1.41177570E+01;
COFD[489] = 3.04580143E+00;
COFD[490] = -1.78216259E-01;
COFD[491] = 7.60087614E-03;
COFD[492] = -1.41393784E+01;
COFD[493] = 3.04580143E+00;
COFD[494] = -1.78216259E-01;
COFD[495] = 7.60087614E-03;
COFD[496] = -2.04725537E+01;
COFD[497] = 5.25060921E+00;
COFD[498] = -4.39775466E-01;
COFD[499] = 1.79662640E-02;
COFD[500] = -1.53783351E+01;
COFD[501] = 3.34436607E+00;
COFD[502] = -2.15557310E-01;
COFD[503] = 9.15916913E-03;
COFD[504] = -1.51258950E+01;
COFD[505] = 3.25013576E+00;
COFD[506] = -2.03913392E-01;
COFD[507] = 8.67915586E-03;
COFD[508] = -1.53868866E+01;
COFD[509] = 3.34436607E+00;
COFD[510] = -2.15557310E-01;
COFD[511] = 9.15916913E-03;
COFD[512] = -1.53950670E+01;
COFD[513] = 3.34436607E+00;
COFD[514] = -2.15557310E-01;
COFD[515] = 9.15916913E-03;
COFD[516] = -1.41393784E+01;
COFD[517] = 3.04580143E+00;
COFD[518] = -1.78216259E-01;
COFD[519] = 7.60087614E-03;
COFD[520] = -1.62043224E+01;
COFD[521] = 3.63818936E+00;
COFD[522] = -2.52168495E-01;
COFD[523] = 1.06812174E-02;
COFD[524] = -1.03014237E+01;
COFD[525] = 1.89006574E+00;
COFD[526] = -3.12440542E-02;
COFD[527] = 1.37570633E-03;
COFD[528] = -9.01870427E+00;
COFD[529] = 1.86885237E+00;
COFD[530] = -2.79898003E-02;
COFD[531] = 1.21327131E-03;
COFD[532] = -9.64985378E+00;
COFD[533] = 1.96102826E+00;
COFD[534] = -4.58416954E-02;
COFD[535] = 2.27372333E-03;
COFD[536] = -1.02432185E+01;
COFD[537] = 2.01382318E+00;
COFD[538] = -5.07824773E-02;
COFD[539] = 2.38045050E-03;
COFD[540] = -1.02491843E+01;
COFD[541] = 2.01382318E+00;
COFD[542] = -5.07824773E-02;
COFD[543] = 2.38045050E-03;
COFD[544] = -1.15042109E+01;
COFD[545] = 2.41805816E+00;
COFD[546] = -9.70163830E-02;
COFD[547] = 4.09454835E-03;
COFD[548] = -1.03186358E+01;
COFD[549] = 1.89854108E+00;
COFD[550] = -3.33173364E-02;
COFD[551] = 1.51065205E-03;
COFD[552] = -1.04392040E+01;
COFD[553] = 1.93661815E+00;
COFD[554] = -3.90724483E-02;
COFD[555] = 1.79631123E-03;
COFD[556] = -1.03203364E+01;
COFD[557] = 1.89854108E+00;
COFD[558] = -3.33173364E-02;
COFD[559] = 1.51065205E-03;
COFD[560] = -1.03219415E+01;
COFD[561] = 1.89854108E+00;
COFD[562] = -3.33173364E-02;
COFD[563] = 1.51065205E-03;
COFD[564] = -1.02491843E+01;
COFD[565] = 2.01382318E+00;
COFD[566] = -5.07824773E-02;
COFD[567] = 2.38045050E-03;
COFD[568] = -1.03014237E+01;
COFD[569] = 1.89006574E+00;
COFD[570] = -3.12440542E-02;
COFD[571] = 1.37570633E-03;
COFD[572] = -1.19059959E+01;
COFD[573] = 2.96079381E+00;
COFD[574] = -1.90702335E-01;
COFD[575] = 9.23708173E-03;
}


// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
KTDIF[0] = 0;
KTDIF[1] = 1;
KTDIF[2] = 11;
}


// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
COFTD[0] = 0.00000000E+00;
COFTD[1] = 0.00000000E+00;
COFTD[2] = 0.00000000E+00;
COFTD[3] = 0.00000000E+00;
COFTD[4] = 1.66147177E-01;
COFTD[5] = 1.85683391E-05;
COFTD[6] = -6.10633845E-09;
COFTD[7] = 7.25879508E-13;
COFTD[8] = 3.38241781E-01;
COFTD[9] = 1.85223136E-04;
COFTD[10] = -6.89052158E-08;
COFTD[11] = 7.81498433E-12;
COFTD[12] = 3.40786931E-01;
COFTD[13] = 1.86616874E-04;
COFTD[14] = -6.94237033E-08;
COFTD[15] = 7.87378934E-12;
COFTD[16] = -5.62837587E-02;
COFTD[17] = 5.52048103E-04;
COFTD[18] = -1.73183201E-07;
COFTD[19] = 1.71537861E-11;
COFTD[20] = 3.04088380E-01;
COFTD[21] = 2.65270116E-04;
COFTD[22] = -9.70943091E-08;
COFTD[23] = 1.07635184E-11;
COFTD[24] = 3.20290657E-01;
COFTD[25] = 2.40548803E-04;
COFTD[26] = -8.85885478E-08;
COFTD[27] = 9.89228612E-12;
COFTD[28] = 3.04674614E-01;
COFTD[29] = 2.65781515E-04;
COFTD[30] = -9.72814917E-08;
COFTD[31] = 1.07842688E-11;
COFTD[32] = 3.05227103E-01;
COFTD[33] = 2.66263476E-04;
COFTD[34] = -9.74578994E-08;
COFTD[35] = 1.08038247E-11;
COFTD[36] = 3.40786931E-01;
COFTD[37] = 1.86616874E-04;
COFTD[38] = -6.94237033E-08;
COFTD[39] = 7.87378934E-12;
COFTD[40] = 2.57370633E-01;
COFTD[41] = 3.25808299E-04;
COFTD[42] = -1.17191299E-07;
COFTD[43] = 1.27685314E-11;
COFTD[44] = 3.26769943E-01;
COFTD[45] = -7.04642731E-06;
COFTD[46] = 1.27278127E-08;
COFTD[47] = -2.88388137E-12;
COFTD[48] = -1.66147177E-01;
COFTD[49] = -1.85683391E-05;
COFTD[50] = 6.10633845E-09;
COFTD[51] = -7.25879508E-13;
COFTD[52] = 0.00000000E+00;
COFTD[53] = 0.00000000E+00;
COFTD[54] = 0.00000000E+00;
COFTD[55] = 0.00000000E+00;
COFTD[56] = 4.22959636E-01;
COFTD[57] = -1.12798929E-05;
COFTD[58] = 1.14156290E-08;
COFTD[59] = -1.73915067E-12;
COFTD[60] = 4.29422137E-01;
COFTD[61] = -1.14522411E-05;
COFTD[62] = 1.15900511E-08;
COFTD[63] = -1.76572357E-12;
COFTD[64] = 1.43934529E-01;
COFTD[65] = 3.48301797E-04;
COFTD[66] = -1.21363557E-07;
COFTD[67] = 1.28854938E-11;
COFTD[68] = 4.64675479E-01;
COFTD[69] = 1.14029890E-05;
COFTD[70] = 6.58530388E-10;
COFTD[71] = -2.03731242E-13;
COFTD[72] = 4.62401082E-01;
COFTD[73] = 1.96551071E-06;
COFTD[74] = 5.10373122E-09;
COFTD[75] = -8.08459096E-13;
COFTD[76] = 4.66474057E-01;
COFTD[77] = 1.14471255E-05;
COFTD[78] = 6.61079303E-10;
COFTD[79] = -2.04519806E-13;
COFTD[80] = 4.68171998E-01;
COFTD[81] = 1.14887925E-05;
COFTD[82] = 6.63485597E-10;
COFTD[83] = -2.05264248E-13;
COFTD[84] = 4.29422137E-01;
COFTD[85] = -1.14522411E-05;
COFTD[86] = 1.15900511E-08;
COFTD[87] = -1.76572357E-12;
COFTD[88] = 4.56551121E-01;
COFTD[89] = 4.15832915E-05;
COFTD[90] = -1.27925209E-08;
COFTD[91] = 1.50962262E-12;
COFTD[92] = 1.48732326E-01;
COFTD[93] = 8.04670424E-05;
COFTD[94] = -3.76755350E-08;
COFTD[95] = 1.81426636E-12;
COFTD[96] = -3.26769943E-01;
COFTD[97] = 7.04642731E-06;
COFTD[98] = -1.27278127E-08;
COFTD[99] = 2.88388137E-12;
COFTD[100] = -1.48732326E-01;
COFTD[101] = -8.04670424E-05;
COFTD[102] = 3.76755350E-08;
COFTD[103] = -1.81426636E-12;
COFTD[104] = 3.07129180E-01;
COFTD[105] = 4.11800129E-05;
COFTD[106] = -7.42984800E-09;
COFTD[107] = -1.92962221E-12;
COFTD[108] = 3.16962164E-01;
COFTD[109] = 4.24984236E-05;
COFTD[110] = -7.66772046E-09;
COFTD[111] = -1.99140060E-12;
COFTD[112] = 3.12914759E-01;
COFTD[113] = 4.16762509E-05;
COFTD[114] = -1.42870137E-08;
COFTD[115] = 1.70110251E-12;
COFTD[116] = 4.14043165E-01;
COFTD[117] = 1.57960962E-05;
COFTD[118] = 7.26518368E-09;
COFTD[119] = -3.52385131E-12;
COFTD[120] = 3.94752941E-01;
COFTD[121] = 2.57857550E-05;
COFTD[122] = 2.52233543E-09;
COFTD[123] = -3.16571792E-12;
COFTD[124] = 4.17267864E-01;
COFTD[125] = 1.59191212E-05;
COFTD[126] = 7.32176723E-09;
COFTD[127] = -3.55129618E-12;
COFTD[128] = 4.20321560E-01;
COFTD[129] = 1.60356223E-05;
COFTD[130] = 7.37535019E-09;
COFTD[131] = -3.57728566E-12;
COFTD[132] = 3.16962164E-01;
COFTD[133] = 4.24984236E-05;
COFTD[134] = -7.66772046E-09;
COFTD[135] = -1.99140060E-12;
COFTD[136] = 4.45321260E-01;
COFTD[137] = -5.55115003E-06;
COFTD[138] = 1.60873845E-08;
COFTD[139] = -3.96076426E-12;
COFTD[140] = 0.00000000E+00;
COFTD[141] = 0.00000000E+00;
COFTD[142] = 0.00000000E+00;
COFTD[143] = 0.00000000E+00;
}


// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i)
{

amrex::Real   EPS[12];
amrex::Real   SIG[12];
amrex::Real    wt[12];
amrex::Real Rcst = 83.144598;// in bar [CGS] !
amrex::Real avogadro = 6.02214199e23;
amrex::Real boltzmann = 1.3806503e-16;// we work in CGS

egtransetEPS(EPS);
egtransetSIG(SIG);
get_mw(wt);

// species 0: H
Tci[0] = 1.316 * EPS[0] ; 
ai[0] = (5.55 * avogadro * avogadro * EPS[0]*boltzmann * 1e-24 * SIG[0] * SIG[0] * SIG[0] ) / (wt[0] * wt[0]); 
bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]); 
acentric_i[0] = 0.0 ;

// species 1: H2
// Imported from NIST
Tci[1] = 33.145000 ; 
ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] / (2.015880 * 2.015880 * 12.964000); 
bi[1] = 0.08664 * Rcst * Tci[1] / (2.015880 * 12.964000); 
acentric_i[1] = -0.219000 ;

// species 2: O
Tci[2] = 1.316 * EPS[2] ; 
ai[2] = (5.55 * avogadro * avogadro * EPS[2]*boltzmann * 1e-24 * SIG[2] * SIG[2] * SIG[2] ) / (wt[2] * wt[2]); 
bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]); 
acentric_i[2] = 0.0 ;

// species 3: OH
Tci[3] = 1.316 * EPS[3] ; 
ai[3] = (5.55 * avogadro * avogadro * EPS[3]*boltzmann * 1e-24 * SIG[3] * SIG[3] * SIG[3] ) / (wt[3] * wt[3]); 
bi[3] = 0.855 * avogadro * 1e-24 * SIG[3] * SIG[3] * SIG[3] / (wt[3]); 
acentric_i[3] = 0.0 ;

// species 4: H2O
// Imported from NIST
Tci[4] = 647.096000 ; 
ai[4] = 1e6 * 0.42748 * Rcst * Rcst * Tci[4] * Tci[4] / (18.015340 * 18.015340 * 220.640000); 
bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000); 
acentric_i[4] = 0.344300 ;

// species 5: O2
// Imported from NIST
Tci[5] = 154.581000 ; 
ai[5] = 1e6 * 0.42748 * Rcst * Rcst * Tci[5] * Tci[5] / (31.998800 * 31.998800 * 50.430466); 
bi[5] = 0.08664 * Rcst * Tci[5] / (31.998800 * 50.430466); 
acentric_i[5] = 0.022200 ;

// species 6: N2
// Imported from NIST
Tci[6] = 126.192000 ; 
ai[6] = 1e6 * 0.42748 * Rcst * Rcst * Tci[6] * Tci[6] / (28.013400 * 28.013400 * 33.958000); 
bi[6] = 0.08664 * Rcst * Tci[6] / (28.013400 * 33.958000); 
acentric_i[6] = 0.037200 ;

// species 7: HO2
Tci[7] = 1.316 * EPS[7] ; 
ai[7] = (5.55 * avogadro * avogadro * EPS[7]*boltzmann * 1e-24 * SIG[7] * SIG[7] * SIG[7] ) / (wt[7] * wt[7]); 
bi[7] = 0.855 * avogadro * 1e-24 * SIG[7] * SIG[7] * SIG[7] / (wt[7]); 
acentric_i[7] = 0.0 ;

// species 8: H2O2
Tci[8] = 1.316 * EPS[8] ; 
ai[8] = (5.55 * avogadro * avogadro * EPS[8]*boltzmann * 1e-24 * SIG[8] * SIG[8] * SIG[8] ) / (wt[8] * wt[8]); 
bi[8] = 0.855 * avogadro * 1e-24 * SIG[8] * SIG[8] * SIG[8] / (wt[8]); 
acentric_i[8] = 0.0 ;

// species 9: OHV
Tci[9] = 1.316 * EPS[9] ; 
ai[9] = (5.55 * avogadro * avogadro * EPS[9]*boltzmann * 1e-24 * SIG[9] * SIG[9] * SIG[9] ) / (wt[9] * wt[9]); 
bi[9] = 0.855 * avogadro * 1e-24 * SIG[9] * SIG[9] * SIG[9] / (wt[9]); 
acentric_i[9] = 0.0 ;

// species 10: AR
// Imported from NIST
Tci[10] = 150.860000 ; 
ai[10] = 1e6 * 0.42748 * Rcst * Rcst * Tci[10] * Tci[10] / (39.948000 * 39.948000 * 48.980000); 
bi[10] = 0.08664 * Rcst * Tci[10] / (39.948000 * 48.980000); 
acentric_i[10] = -0.002000 ;

// species 11: HE
Tci[11] = 1.316 * EPS[11] ; 
ai[11] = (5.55 * avogadro * avogadro * EPS[11]*boltzmann * 1e-24 * SIG[11] * SIG[11] * SIG[11] ) / (wt[11] * wt[11]); 
bi[11] = 0.855 * avogadro * 1e-24 * SIG[11] * SIG[11] * SIG[11] / (wt[11]); 
acentric_i[11] = 0.0 ;

}


// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS_SRK(amrex::Real *  sqrtOneOverTc, amrex::Real *  sqrtAsti, amrex::Real *  Bi, amrex::Real *  Fomega)
{

// species 0: H
// Computed from Lennard-Jones
sqrtOneOverTc[0] = 7.2391579785311e-02;
sqrtAsti[0] = 5.8451221185793e+05;
Bi[0] = 4.4006634563410e+00;
Fomega[0] = 4.8508000000000e-01;

// species 1: H2
// Imported from NIST
sqrtOneOverTc[1] = 1.7369646834812e-01;
sqrtAsti[1] = 2.4824229305911e+05;
Bi[1] = 9.1362254169923e+00;
Fomega[1] = 1.3798618890700e-01;

// species 2: O
// Computed from Lennard-Jones
sqrtOneOverTc[2] = 9.7460146948181e-02;
sqrtAsti[2] = 4.2500156012386e+04;
Bi[2] = 6.6930234851885e-01;
Fomega[2] = 4.8508000000000e-01;

// species 3: OH
// Computed from Lennard-Jones
sqrtOneOverTc[3] = 9.7460146948181e-02;
sqrtAsti[3] = 3.9981183985545e+04;
Bi[3] = 6.2963299076575e-01;
Fomega[3] = 4.8508000000000e-01;

// species 4: H2O
// Imported from NIST
sqrtOneOverTc[4] = 3.9311140369715e-02;
sqrtAsti[4] = 1.3145468419592e+05;
Bi[4] = 1.1727204136223e+00;
Fomega[4] = 1.0013577274636e+00;

// species 5: O2
// Imported from NIST
sqrtOneOverTc[5] = 8.0430717653241e-02;
sqrtAsti[5] = 3.6980081924740e+04;
Bi[5] = 6.9005216444999e-01;
Fomega[5] = 5.1945301904908e-01;

// species 6: N2
// Imported from NIST
sqrtOneOverTc[6] = 8.9019282240563e-02;
sqrtAsti[6] = 4.2022983929562e+04;
Bi[6] = 9.5560052221563e-01;
Fomega[6] = 5.4259343186608e-01;

// species 7: HO2
// Computed from Lennard-Jones
sqrtOneOverTc[7] = 8.4114315417892e-02;
sqrtAsti[7] = 3.3657825710980e+04;
Bi[7] = 6.4505861132042e-01;
Fomega[7] = 4.8508000000000e-01;

// species 8: H2O2
// Computed from Lennard-Jones
sqrtOneOverTc[8] = 8.4114315417892e-02;
sqrtAsti[8] = 3.2660380884830e+04;
Bi[8] = 6.2594239211036e-01;
Fomega[8] = 4.8508000000000e-01;

// species 9: OHV
// Computed from Lennard-Jones
sqrtOneOverTc[9] = 9.7460146948181e-02;
sqrtAsti[9] = 3.9981183985545e+04;
Bi[9] = 6.2963299076575e-01;
Fomega[9] = 4.8508000000000e-01;

// species 10: AR
// Imported from NIST
sqrtOneOverTc[10] = 8.1416597423011e-02;
sqrtAsti[10] = 2.9333346340060e+04;
Bi[10] = 5.5540880812616e-01;
Fomega[10] = 4.8197599354800e-01;

// species 11: HE
// Computed from Lennard-Jones
sqrtOneOverTc[11] = 2.7294300551596e-01;
sqrtAsti[11] = 5.4994061413731e+04;
Bi[11] = 2.1989344021800e+00;
Fomega[11] = 4.8508000000000e-01;

}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
