#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {21,33,34,35,36,37,38,45,66,91,111,113,114,115,151,152,166,176,177,39,0,1,5,6,8,9,81,164,2,3,4,7,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,32,40,41,42,43,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,67,68,69,70,71,72,73,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,90,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,153,154,155,156,157,158,159,160,161,162,163,165,167,168,169,170,171,172,173,174,175,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200};

// Returns 0-based map of reaction order
void GET_RMAP
(int * _rmap)
{
for (int j=0; j<NUM_REACTIONS; ++j)
{
_rmap[j] = rmap[j];
}
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int& nspec, int ki[], int nu[])
{
const int ns[NUM_GAS_REACTIONS] =
     {2,2,4,4,4,3,3,3,3,3,4,4,4,3,4,4,2,4,4,4,4,2,4,4,4,4,4,3,4,4,4,3,3,3,3,3,3,3,3,3,4,4,4,3,3,3,4,4,4,4,4,4,3,4,4,4,4,4,5,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,2,4,3,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,3,4,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,3,5,4,4,5,4,4,3,3,4,4,4,4,4,4,4,4,4,4,5,3,5,3,4,4,4,4,4,3,4,5,3,3,3,3,3,3,2,4,3,3,3,3,4,4,4,4,4,4,4,4,5,4,4,3,4,4};
const int kiv[NUM_GAS_REACTIONS*5] =
     {1,0,0,0,0,2,5,0,0,0,1,2,0,3,0,1,3,0,4,0,0,5,2,3,0,0,3,4,0,0,0,3,4,0,0,4,2,3,0,0,0,2,3,0,0,0,2,36,0,0,4,36,4,3,0,1,36,1,3,0,6,36,6,3,0,3,36,3,0,0,0,36,0,3,0,37,36,37,3,0,36,3,0,0,0,5,36,5,3,0,11,36,11,3,0,10,36,10,3,0,33,36,33,3,0,8,3,0,0,0,0,8,4,3,0,0,8,1,7,0,8,2,7,3,0,8,3,4,7,0,8,3,4,7,0,0,7,3,0,0,0,7,1,5,0,7,2,5,3,0,7,3,4,5,0,7,8,5,0,0,7,8,5,0,0,0,5,7,0,0,0,5,7,0,0,0,5,7,0,0,0,5,7,0,0,0,5,7,0,0,0,5,7,0,0,10,2,11,0,0,10,3,11,0,0,10,7,11,3,0,10,5,11,2,0,15,10,3,0,0,15,11,0,0,0,29,0,33,0,0,29,7,33,5,0,33,0,29,1,0,33,2,29,3,0,33,3,29,4,0,33,7,29,8,0,27,33,29,12,0,21,33,29,0,0,25,6,21,6,0,37,25,37,21,0,25,4,21,4,0,25,10,21,10,0,25,11,21,11,0,25,5,10,0,3,25,5,10,4,0,25,2,10,1,0,25,2,0,14,0,25,1,29,0,0,25,0,13,1,0,25,3,18,0,0,25,11,18,10,0,21,0,29,0,0,21,5,11,0,0,21,5,18,2,0,21,2,10,0,0,21,0,13,1,0,21,3,13,4,0,37,17,37,13,0,17,4,13,4,0,17,10,13,10,0,17,11,13,11,0,17,5,13,5,0,17,1,13,1,0,33,17,13,33,0,17,13,0,0,0,17,6,13,6,0,9,0,17,0,0,13,5,10,36,0,13,5,14,2,0,13,2,10,0,0,13,0,9,1,0,13,3,0,14,0,13,4,18,0,0,13,11,10,14,0,9,3,10,0,0,9,5,10,2,0,29,5,27,0,0,29,5,26,2,0,29,5,18,3,0,29,2,18,0,0,29,3,25,4,0,29,3,18,1,0,29,3,30,0,0,29,3,26,0,0,29,3,1,22,0,29,3,21,4,0,29,7,26,3,0,27,2,26,5,0,27,0,26,3,0,27,3,34,5,0,27,7,12,5,0,27,8,12,7,0,29,27,26,0,0,27,18,34,5,0,27,26,5,0,0,27,1,12,0,0,12,26,3,0,0,31,18,3,0,0,34,29,3,0,0,34,25,4,0,0,34,30,0,0,0,34,5,30,7,0,34,0,30,1,0,34,2,30,3,0,34,3,30,4,0,34,7,30,8,0,29,34,30,33,0,34,14,18,30,0,26,34,30,34,0,27,34,30,12,0,34,5,26,7,0,34,0,26,1,0,34,2,26,3,0,34,3,26,4,0,34,7,26,8,0,29,34,26,33,0,30,5,18,7,0,30,5,18,7,0,30,0,18,1,0,30,7,18,8,0,30,26,18,34,0,30,3,18,4,0,30,2,18,3,0,30,18,34,0,0,30,7,35,3,0,26,5,18,7,0,26,0,18,1,0,26,7,18,8,0,29,26,18,33,0,26,18,34,0,0,22,5,11,0,3,22,5,11,4,0,22,2,11,0,0,22,2,10,0,3,0,22,18,0,0,22,3,4,14,0,0,14,18,0,0,10,1,18,0,0,18,5,14,7,0,18,2,14,3,0,18,0,1,14,0,18,3,4,14,0,18,7,8,14,0,18,29,33,14,0,18,16,14,20,0,18,19,14,23,0,18,26,34,14,0,18,27,12,14,0,18,3,10,0,4,18,10,1,0,0,18,2,10,0,3,14,10,0,0,0,14,5,10,7,0,14,2,10,3,0,0,14,10,1,0,14,3,10,4,0,29,14,33,10,0,14,18,10,0,0,14,2,11,0,0,14,7,11,0,3,14,10,1,0,0,18,0,30,0,0,26,18,0,0,0,18,3,35,0,0,35,0,23,0,0,18,7,28,0,0,28,24,0,0,0,7,24,32,5,0,35,3,32,0,0,14,5,16,0,0,23,10,4,0,0,23,11,1,0,0,7,15,23,5,0,23,5,7,19,0,0,23,1,15,0,0,23,1,19,0,23,2,15,3,0,23,2,19,3,0,23,3,4,15,0,23,3,4,19,0,29,23,33,10,3,7,23,8,15,0,7,23,8,19,0,19,11,0,0,0,5,19,11,7,0,26,14,34,10,0};
const int nuv[NUM_GAS_REACTIONS*5] =
     {-1,2,0,0,0,-2,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,1,0,-2,2,1,0,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,1,1,0,0,-1,-1,1,1,1,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,1,-2,2,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0};
if (i < 1) {
// Return max num species per reaction
nspec = 5;
} else {
if (i > NUM_GAS_REACTIONS) {
nspec = -1;
} else {
nspec = ns[i-1];
for (int j=0; j<nspec; ++j) {
ki[j] = kiv[(i-1)*5 + j] + 1;
nu[j] = nuv[(i-1)*5 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[], amrex::Real q_f[], amrex::Real q_r[])
{
amrex::Real c[39]; // temporary storage
amrex::Real PORT = 1e6 * P/(8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (int id = 0; id < 39; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, T);

// convert to chemkin units
for (int id = 0; id < 201; ++id) {
q_f[id] *= 1.0e-6;
q_r[id] *= 1.0e-6;
}
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
const amrex::Real invT = 1.0 / T;
const amrex::Real logT = log(T);
// compute the Gibbs free energy
amrex::Real g_RT[39];
gibbs(g_RT, T);

amrex::Real sc_qss[1];
comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);

}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 12.011000; // C
awt[1] = 1.008000; // H
awt[2] = 14.007000; // N
awt[3] = 15.999000; // O
awt[4] = 39.950000; // Ar
awt[5] = 4.002602; // He
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int kd = 6; 
// Zero ncf
for (int id = 0; id < kd * 39; ++ id) {
 ncf[id] = 0; 
}

// H
ncf[ 0 * kd + 1 ] = 1; // H

// H2
ncf[ 1 * kd + 1 ] = 2; // H

// O
ncf[ 2 * kd + 3 ] = 1; // O

// OH
ncf[ 3 * kd + 1 ] = 1; // H
ncf[ 3 * kd + 3 ] = 1; // O

// H2O
ncf[ 4 * kd + 1 ] = 2; // H
ncf[ 4 * kd + 3 ] = 1; // O

// O2
ncf[ 5 * kd + 3 ] = 2; // O

// N2
ncf[ 6 * kd + 2 ] = 2; // N

// HO2
ncf[ 7 * kd + 1 ] = 1; // H
ncf[ 7 * kd + 3 ] = 2; // O

// H2O2
ncf[ 8 * kd + 1 ] = 2; // H
ncf[ 8 * kd + 3 ] = 2; // O

// C
ncf[ 9 * kd + 0 ] = 1; // C

// CO
ncf[ 10 * kd + 0 ] = 1; // C
ncf[ 10 * kd + 3 ] = 1; // O

// CO2
ncf[ 11 * kd + 0 ] = 1; // C
ncf[ 11 * kd + 3 ] = 2; // O

// CH3O2H
ncf[ 12 * kd + 0 ] = 1; // C
ncf[ 12 * kd + 1 ] = 4; // H
ncf[ 12 * kd + 3 ] = 2; // O

// CH
ncf[ 13 * kd + 0 ] = 1; // C
ncf[ 13 * kd + 1 ] = 1; // H

// HCO
ncf[ 14 * kd + 0 ] = 1; // C
ncf[ 14 * kd + 1 ] = 1; // H
ncf[ 14 * kd + 3 ] = 1; // O

// HOCO
ncf[ 15 * kd + 0 ] = 1; // C
ncf[ 15 * kd + 1 ] = 1; // H
ncf[ 15 * kd + 3 ] = 2; // O

// O2CHO
ncf[ 16 * kd + 0 ] = 1; // C
ncf[ 16 * kd + 1 ] = 1; // H
ncf[ 16 * kd + 3 ] = 3; // O

// CHV
ncf[ 17 * kd + 0 ] = 1; // C
ncf[ 17 * kd + 1 ] = 1; // H

// CH2O
ncf[ 18 * kd + 0 ] = 1; // C
ncf[ 18 * kd + 1 ] = 2; // H
ncf[ 18 * kd + 3 ] = 1; // O

// OCHO
ncf[ 19 * kd + 0 ] = 1; // C
ncf[ 19 * kd + 1 ] = 1; // H
ncf[ 19 * kd + 3 ] = 2; // O

// HO2CHO
ncf[ 20 * kd + 0 ] = 1; // C
ncf[ 20 * kd + 1 ] = 2; // H
ncf[ 20 * kd + 3 ] = 3; // O

// CH2
ncf[ 21 * kd + 0 ] = 1; // C
ncf[ 21 * kd + 1 ] = 2; // H

// HCOH
ncf[ 22 * kd + 0 ] = 1; // C
ncf[ 22 * kd + 1 ] = 2; // H
ncf[ 22 * kd + 3 ] = 1; // O

// HOCHO
ncf[ 23 * kd + 0 ] = 1; // C
ncf[ 23 * kd + 1 ] = 2; // H
ncf[ 23 * kd + 3 ] = 2; // O

// HOCH2O2
ncf[ 24 * kd + 0 ] = 1; // C
ncf[ 24 * kd + 1 ] = 3; // H
ncf[ 24 * kd + 3 ] = 3; // O

// CH2(S)
ncf[ 25 * kd + 0 ] = 1; // C
ncf[ 25 * kd + 1 ] = 2; // H

// CH3O
ncf[ 26 * kd + 0 ] = 1; // C
ncf[ 26 * kd + 1 ] = 3; // H
ncf[ 26 * kd + 3 ] = 1; // O

// CH3O2
ncf[ 27 * kd + 0 ] = 1; // C
ncf[ 27 * kd + 1 ] = 3; // H
ncf[ 27 * kd + 3 ] = 2; // O

// OCH2O2H
ncf[ 28 * kd + 0 ] = 1; // C
ncf[ 28 * kd + 1 ] = 3; // H
ncf[ 28 * kd + 3 ] = 3; // O

// CH3
ncf[ 29 * kd + 0 ] = 1; // C
ncf[ 29 * kd + 1 ] = 3; // H

// CH2OH
ncf[ 30 * kd + 0 ] = 1; // C
ncf[ 30 * kd + 1 ] = 3; // H
ncf[ 30 * kd + 3 ] = 1; // O

// CH2O2H
ncf[ 31 * kd + 0 ] = 1; // C
ncf[ 31 * kd + 1 ] = 3; // H
ncf[ 31 * kd + 3 ] = 2; // O

// HOCH2O2H
ncf[ 32 * kd + 0 ] = 1; // C
ncf[ 32 * kd + 1 ] = 4; // H
ncf[ 32 * kd + 3 ] = 3; // O

// CH4
ncf[ 33 * kd + 0 ] = 1; // C
ncf[ 33 * kd + 1 ] = 4; // H

// CH3OH
ncf[ 34 * kd + 0 ] = 1; // C
ncf[ 34 * kd + 1 ] = 4; // H
ncf[ 34 * kd + 3 ] = 1; // O

// HOCH2O
ncf[ 35 * kd + 0 ] = 1; // C
ncf[ 35 * kd + 1 ] = 3; // H
ncf[ 35 * kd + 3 ] = 2; // O

// OHV
ncf[ 36 * kd + 1 ] = 1; // H
ncf[ 36 * kd + 3 ] = 1; // O

// AR
ncf[ 37 * kd + 4 ] = 1; // Ar

// HE
ncf[ 38 * kd + 5 ] = 1; // He

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(6);
ename[0] = "C";
ename[1] = "H";
ename[2] = "N";
ename[3] = "O";
ename[4] = "Ar";
ename[5] = "He";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(39);
kname[0] = "H";
kname[1] = "H2";
kname[2] = "O";
kname[3] = "OH";
kname[4] = "H2O";
kname[5] = "O2";
kname[6] = "N2";
kname[7] = "HO2";
kname[8] = "H2O2";
kname[9] = "C";
kname[10] = "CO";
kname[11] = "CO2";
kname[12] = "CH3O2H";
kname[13] = "CH";
kname[14] = "HCO";
kname[15] = "HOCO";
kname[16] = "O2CHO";
kname[17] = "CHV";
kname[18] = "CH2O";
kname[19] = "OCHO";
kname[20] = "HO2CHO";
kname[21] = "CH2";
kname[22] = "HCOH";
kname[23] = "HOCHO";
kname[24] = "HOCH2O2";
kname[25] = "CH2(S)";
kname[26] = "CH3O";
kname[27] = "CH3O2";
kname[28] = "OCH2O2H";
kname[29] = "CH3";
kname[30] = "CH2OH";
kname[31] = "CH2O2H";
kname[32] = "HOCH2O2H";
kname[33] = "CH4";
kname[34] = "CH3OH";
kname[35] = "HOCH2O";
kname[36] = "OHV";
kname[37] = "AR";
kname[38] = "HE";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<40; k++) {
for (int l=0; l<40; l++) {
if(Jac[ 40 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<40; k++) {
for (int l=0; l<40; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 40 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<40; k++) {
for (int l=0; l<40; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 40 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 40;
int offset_col = nc * 40;
for (int k=0; k<40; k++) {
for (int l=0; l<40; l++) {
if(Jac[40*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 40;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 40;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 40;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 40;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<40; k++) {
for (int l=0; l<40; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 40*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[40*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 40*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,1600> Jac = {0.0};
amrex::GpuArray<amrex::Real,39> conc = {0.0};
for (int n=0; n<39; n++) {
    conc[n] = 1.0/ 39.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<40; l++) {
for (int k=0; k<40; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[40*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
