
#include <boost/multiprecision/mpfr.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <cmath>
#include <vector>
#include <stack>
#include <stdexcept>
#include <Eigen/Dense>
#include <stdlib.h>
#include <iostream>

using namespace std;
using namespace Eigen;
namespace py=pybind11;

/**From kmnam MarkovDigraphs module
 * Compute the nullspace of A by performing a singular value decomposition.
 * 
 * This function returns the column of V in the SVD of A = USV corresponding
 * to the least singular value (recall that singular values are always
 * non-negative). It therefore effectively assumes that the A has a 
 * nullspace of dimension one.
 * 
 * @param A Input matrix to be decomposed. 
 * @returns The column of V in the singular value decomposition A = USV 
 *          corresponding to the least singular value. 
 */

template <typename T>
Matrix<T, Dynamic, 1> getOneDimNullspaceFromSVD(const Ref<const Matrix<T, Dynamic, Dynamic> >& A)
{
    // Perform a singular value decomposition of A, only computing V in full
    Eigen::BDCSVD<Matrix<T, Dynamic, Dynamic> > svd(A, Eigen::ComputeFullV);

    // Return the column of V corresponding to the least singular value of A
    // (always given last in the decomposition) 
    Matrix<T, Dynamic, 1> singular = svd.singularValues(); 
    Matrix<T, Dynamic, Dynamic> V = svd.matrixV();
    return V.col(singular.rows() - 1); 
}

typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50> > InternalType;
Matrix<InternalType, Dynamic, 1> getrhos(py::array_t<double> parsar, py::array_t<double> varvals){
    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    auto varsarbuf=varvals.request();
    double *vars=(double *) varsarbuf.ptr;
        long double x=vars[0];
    long double a1_0_1=pars[0];
    long double k_1_0=pars[1];
    long double kr_1_0=pars[2];
    long double a1_0_2=pars[3];
    long double k_2_0=pars[4];
    long double a1_0_3=pars[5];
    long double k_3_0=pars[6];
    long double b1_1_1=pars[7];
    long double k_1_1=pars[8];
    long double kr_1_1=pars[9];
    long double b1_1_2=pars[10];
    long double k_2_1=pars[11];
    long double b1_1_3=pars[12];
    long double k_3_1=pars[13];
    long double a2_0_1=pars[14];
    long double a2_0_2=pars[15];
    long double a2_0_3=pars[16];
    long double b2_2_1=pars[17];
    long double k_1_2=pars[18];
    long double kr_1_2=pars[19];
    long double b2_2_2=pars[20];
    long double k_2_2=pars[21];
    long double b2_2_3=pars[22];
    long double k_3_2=pars[23];
    long double a3_0_1=pars[24];
    long double a3_0_2=pars[25];
    long double a3_0_3=pars[26];
    long double b3_3_1=pars[27];
    long double k_1_3=pars[28];
    long double kr_1_3=pars[29];
    long double b3_3_2=pars[30];
    long double k_2_3=pars[31];
    long double b3_3_3=pars[32];
    long double k_3_3=pars[33];
    long double a4_0_1=pars[34];
    long double a4_0_2=pars[35];
    long double a4_0_3=pars[36];
    long double b4_4_1=pars[37];
    long double k_1_4=pars[38];
    long double kr_1_4=pars[39];
    long double b4_4_2=pars[40];
    long double k_2_4=pars[41];
    long double b4_4_3=pars[42];
    long double k_3_4=pars[43];
    long double a2_1_1=pars[44];
    long double a2_1_2=pars[45];
    long double a2_1_3=pars[46];
    long double b2_1U2_1=pars[47];
    long double k_1_1U2=pars[48];
    long double kr_1_1U2=pars[49];
    long double b2_1U2_2=pars[50];
    long double k_2_1U2=pars[51];
    long double b2_1U2_3=pars[52];
    long double k_3_1U2=pars[53];
    long double a3_1_1=pars[54];
    long double a3_1_2=pars[55];
    long double a3_1_3=pars[56];
    long double b3_1U3_1=pars[57];
    long double k_1_1U3=pars[58];
    long double kr_1_1U3=pars[59];
    long double b3_1U3_2=pars[60];
    long double k_2_1U3=pars[61];
    long double b3_1U3_3=pars[62];
    long double k_3_1U3=pars[63];
    long double a4_1_1=pars[64];
    long double a4_1_2=pars[65];
    long double a4_1_3=pars[66];
    long double b4_1U4_1=pars[67];
    long double k_1_1U4=pars[68];
    long double kr_1_1U4=pars[69];
    long double b4_1U4_2=pars[70];
    long double k_2_1U4=pars[71];
    long double b4_1U4_3=pars[72];
    long double k_3_1U4=pars[73];
    long double a1_2_1=pars[74];
    long double a1_2_2=pars[75];
    long double a1_2_3=pars[76];
    long double b1_1U2_1=pars[77];
    long double b1_1U2_2=pars[78];
    long double b1_1U2_3=pars[79];
    long double a3_2_1=pars[80];
    long double a3_2_2=pars[81];
    long double a3_2_3=pars[82];
    long double b3_2U3_1=pars[83];
    long double k_1_2U3=pars[84];
    long double kr_1_2U3=pars[85];
    long double b3_2U3_2=pars[86];
    long double k_2_2U3=pars[87];
    long double b3_2U3_3=pars[88];
    long double k_3_2U3=pars[89];
    long double a4_2_1=pars[90];
    long double a4_2_2=pars[91];
    long double a4_2_3=pars[92];
    long double b4_2U4_1=pars[93];
    long double k_1_2U4=pars[94];
    long double kr_1_2U4=pars[95];
    long double b4_2U4_2=pars[96];
    long double k_2_2U4=pars[97];
    long double b4_2U4_3=pars[98];
    long double k_3_2U4=pars[99];
    long double a1_3_1=pars[100];
    long double a1_3_2=pars[101];
    long double a1_3_3=pars[102];
    long double b1_1U3_1=pars[103];
    long double b1_1U3_2=pars[104];
    long double b1_1U3_3=pars[105];
    long double a2_3_1=pars[106];
    long double a2_3_2=pars[107];
    long double a2_3_3=pars[108];
    long double b2_2U3_1=pars[109];
    long double b2_2U3_2=pars[110];
    long double b2_2U3_3=pars[111];
    long double a4_3_1=pars[112];
    long double a4_3_2=pars[113];
    long double a4_3_3=pars[114];
    long double b4_3U4_1=pars[115];
    long double k_1_3U4=pars[116];
    long double kr_1_3U4=pars[117];
    long double b4_3U4_2=pars[118];
    long double k_2_3U4=pars[119];
    long double b4_3U4_3=pars[120];
    long double k_3_3U4=pars[121];
    long double a1_4_1=pars[122];
    long double a1_4_2=pars[123];
    long double a1_4_3=pars[124];
    long double b1_1U4_1=pars[125];
    long double b1_1U4_2=pars[126];
    long double b1_1U4_3=pars[127];
    long double a2_4_1=pars[128];
    long double a2_4_2=pars[129];
    long double a2_4_3=pars[130];
    long double b2_2U4_1=pars[131];
    long double b2_2U4_2=pars[132];
    long double b2_2U4_3=pars[133];
    long double a3_4_1=pars[134];
    long double a3_4_2=pars[135];
    long double a3_4_3=pars[136];
    long double b3_3U4_1=pars[137];
    long double b3_3U4_2=pars[138];
    long double b3_3U4_3=pars[139];
    long double a3_1U2_1=pars[140];
    long double a3_1U2_2=pars[141];
    long double a3_1U2_3=pars[142];
    long double b3_1U2U3_1=pars[143];
    long double k_1_1U2U3=pars[144];
    long double kr_1_1U2U3=pars[145];
    long double b3_1U2U3_2=pars[146];
    long double k_2_1U2U3=pars[147];
    long double b3_1U2U3_3=pars[148];
    long double k_3_1U2U3=pars[149];
    long double a4_1U2_1=pars[150];
    long double a4_1U2_2=pars[151];
    long double a4_1U2_3=pars[152];
    long double b4_1U2U4_1=pars[153];
    long double k_1_1U2U4=pars[154];
    long double kr_1_1U2U4=pars[155];
    long double b4_1U2U4_2=pars[156];
    long double k_2_1U2U4=pars[157];
    long double b4_1U2U4_3=pars[158];
    long double k_3_1U2U4=pars[159];
    long double a2_1U3_1=pars[160];
    long double a2_1U3_2=pars[161];
    long double a2_1U3_3=pars[162];
    long double b2_1U2U3_1=pars[163];
    long double b2_1U2U3_2=pars[164];
    long double b2_1U2U3_3=pars[165];
    long double a4_1U3_1=pars[166];
    long double a4_1U3_2=pars[167];
    long double a4_1U3_3=pars[168];
    long double b4_1U3U4_1=pars[169];
    long double k_1_1U3U4=pars[170];
    long double kr_1_1U3U4=pars[171];
    long double b4_1U3U4_2=pars[172];
    long double k_2_1U3U4=pars[173];
    long double b4_1U3U4_3=pars[174];
    long double k_3_1U3U4=pars[175];
    long double a2_1U4_1=pars[176];
    long double a2_1U4_2=pars[177];
    long double a2_1U4_3=pars[178];
    long double b2_1U2U4_1=pars[179];
    long double b2_1U2U4_2=pars[180];
    long double b2_1U2U4_3=pars[181];
    long double a3_1U4_1=pars[182];
    long double a3_1U4_2=pars[183];
    long double a3_1U4_3=pars[184];
    long double b3_1U3U4_1=pars[185];
    long double b3_1U3U4_2=pars[186];
    long double b3_1U3U4_3=pars[187];
    long double a1_2U3_1=pars[188];
    long double a1_2U3_2=pars[189];
    long double a1_2U3_3=pars[190];
    long double b1_1U2U3_1=pars[191];
    long double b1_1U2U3_2=pars[192];
    long double b1_1U2U3_3=pars[193];
    long double a4_2U3_1=pars[194];
    long double a4_2U3_2=pars[195];
    long double a4_2U3_3=pars[196];
    long double b4_2U3U4_1=pars[197];
    long double k_1_2U3U4=pars[198];
    long double kr_1_2U3U4=pars[199];
    long double b4_2U3U4_2=pars[200];
    long double k_2_2U3U4=pars[201];
    long double b4_2U3U4_3=pars[202];
    long double k_3_2U3U4=pars[203];
    long double a1_2U4_1=pars[204];
    long double a1_2U4_2=pars[205];
    long double a1_2U4_3=pars[206];
    long double b1_1U2U4_1=pars[207];
    long double b1_1U2U4_2=pars[208];
    long double b1_1U2U4_3=pars[209];
    long double a3_2U4_1=pars[210];
    long double a3_2U4_2=pars[211];
    long double a3_2U4_3=pars[212];
    long double b3_2U3U4_1=pars[213];
    long double b3_2U3U4_2=pars[214];
    long double b3_2U3U4_3=pars[215];
    long double a1_3U4_1=pars[216];
    long double a1_3U4_2=pars[217];
    long double a1_3U4_3=pars[218];
    long double b1_1U3U4_1=pars[219];
    long double b1_1U3U4_2=pars[220];
    long double b1_1U3U4_3=pars[221];
    long double a2_3U4_1=pars[222];
    long double a2_3U4_2=pars[223];
    long double a2_3U4_3=pars[224];
    long double b2_2U3U4_1=pars[225];
    long double b2_2U3U4_2=pars[226];
    long double b2_2U3U4_3=pars[227];
    long double a4_1U2U3_1=pars[228];
    long double a4_1U2U3_2=pars[229];
    long double a4_1U2U3_3=pars[230];
    long double b4_1U2U3U4_1=pars[231];
    long double k_1_1U2U3U4=pars[232];
    long double kr_1_1U2U3U4=pars[233];
    long double b4_1U2U3U4_2=pars[234];
    long double k_2_1U2U3U4=pars[235];
    long double b4_1U2U3U4_3=pars[236];
    long double k_3_1U2U3U4=pars[237];
    long double a3_1U2U4_1=pars[238];
    long double a3_1U2U4_2=pars[239];
    long double a3_1U2U4_3=pars[240];
    long double b3_1U2U3U4_1=pars[241];
    long double b3_1U2U3U4_2=pars[242];
    long double b3_1U2U3U4_3=pars[243];
    long double a2_1U3U4_1=pars[244];
    long double a2_1U3U4_2=pars[245];
    long double a2_1U3U4_3=pars[246];
    long double b2_1U2U3U4_1=pars[247];
    long double b2_1U2U3U4_2=pars[248];
    long double b2_1U2U3U4_3=pars[249];
    long double a1_2U3U4_1=pars[250];
    long double a1_2U3U4_2=pars[251];
    long double a1_2U3U4_3=pars[252];
    long double b1_1U2U3U4_1=pars[253];
    long double b1_1U2U3U4_2=pars[254];
    long double b1_1U2U3U4_3=pars[255];
    Matrix<InternalType, Dynamic, Dynamic> L = Matrix<InternalType, Dynamic, Dynamic>::Zero(48, 48);
    L<<0,b1_1_1,b2_2_1,b3_3_1,b4_4_1,0,0,0,0,0,0,0,0,0,0,0,kr_1_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a1_0_1*x,0,0,0,0,b2_1U2_1,b3_1U3_1,b4_1U4_1,0,0,0,0,0,0,0,0,0,kr_1_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a2_0_1*x,0,0,0,0,b1_1U2_1,0,0,b3_2U3_1,b4_2U4_1,0,0,0,0,0,0,0,0,kr_1_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2,0,0,0,0,0,0,0,0,0,0,0,0,0,
a3_0_1*x,0,0,0,0,0,b1_1U3_1,0,b2_2U3_1,0,b4_3U4_1,0,0,0,0,0,0,0,0,kr_1_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3,0,0,0,0,0,0,0,0,0,0,0,0,
a4_0_1*x,0,0,0,0,0,0,b1_1U4_1,0,b2_2U4_1,b3_3U4_1,0,0,0,0,0,0,0,0,0,kr_1_4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_4,0,0,0,0,0,0,0,0,0,0,0,
0,a2_1_1*x,a1_2_1*x,0,0,0,0,0,0,0,0,b3_1U2U3_1,b4_1U2U4_1,0,0,0,0,0,0,0,0,kr_1_1U2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2,0,0,0,0,0,0,0,0,0,0,
0,a3_1_1*x,0,a1_3_1*x,0,0,0,0,0,0,0,b2_1U2U3_1,0,b4_1U3U4_1,0,0,0,0,0,0,0,0,kr_1_1U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3,0,0,0,0,0,0,0,0,0,
0,a4_1_1*x,0,0,a1_4_1*x,0,0,0,0,0,0,0,b2_1U2U4_1,b3_1U3U4_1,0,0,0,0,0,0,0,0,0,kr_1_1U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U4,0,0,0,0,0,0,0,0,
0,0,a3_2_1*x,a2_3_1*x,0,0,0,0,0,0,0,b1_1U2U3_1,0,0,b4_2U3U4_1,0,0,0,0,0,0,0,0,0,kr_1_2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3,0,0,0,0,0,0,0,
0,0,a4_2_1*x,0,a2_4_1*x,0,0,0,0,0,0,0,b1_1U2U4_1,0,b3_2U3U4_1,0,0,0,0,0,0,0,0,0,0,kr_1_2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U4,0,0,0,0,0,0,
0,0,0,a4_3_1*x,a3_4_1*x,0,0,0,0,0,0,0,0,b1_1U3U4_1,b2_2U3U4_1,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U4,0,0,0,0,0,
0,0,0,0,0,a3_1U2_1*x,a2_1U3_1*x,0,a1_2U3_1*x,0,0,0,0,0,0,b4_1U2U3U4_1,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3,0,0,0,0,
0,0,0,0,0,a4_1U2_1*x,0,a2_1U4_1*x,0,a1_2U4_1*x,0,0,0,0,0,b3_1U2U3U4_1,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U4,0,0,0,
0,0,0,0,0,0,a4_1U3_1*x,a3_1U4_1*x,0,0,a1_3U4_1*x,0,0,0,0,b2_1U2U3U4_1,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U4,0,0,
0,0,0,0,0,0,0,0,a4_2U3_1*x,a3_2U4_1*x,a2_3U4_1*x,0,0,0,0,b1_1U2U3U4_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U4,0,
0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_1*x,a3_1U2U4_1*x,a2_1U3U4_1*x,a1_2U3U4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U4,
k_1_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1_2,b2_2_2,b3_3_2,b4_4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,k_1_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a1_0_2*x,0,0,0,0,b2_1U2_2,b3_1U3_2,b4_1U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,k_1_2,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_0_2*x,0,0,0,0,b1_1U2_2,0,0,b3_2U3_2,b4_2U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,k_1_3,0,0,0,0,0,0,0,0,0,0,0,0,a3_0_2*x,0,0,0,0,0,b1_1U3_2,0,b2_2U3_2,0,b4_3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,k_1_4,0,0,0,0,0,0,0,0,0,0,0,a4_0_2*x,0,0,0,0,0,0,b1_1U4_2,0,b2_2U4_2,b3_3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,k_1_1U2,0,0,0,0,0,0,0,0,0,0,0,a2_1_2*x,a1_2_2*x,0,0,0,0,0,0,0,0,b3_1U2U3_2,b4_1U2U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,k_1_1U3,0,0,0,0,0,0,0,0,0,0,a3_1_2*x,0,a1_3_2*x,0,0,0,0,0,0,0,b2_1U2U3_2,0,b4_1U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,k_1_1U4,0,0,0,0,0,0,0,0,0,a4_1_2*x,0,0,a1_4_2*x,0,0,0,0,0,0,0,b2_1U2U4_2,b3_1U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,k_1_2U3,0,0,0,0,0,0,0,0,0,a3_2_2*x,a2_3_2*x,0,0,0,0,0,0,0,b1_1U2U3_2,0,0,b4_2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,k_1_2U4,0,0,0,0,0,0,0,0,a4_2_2*x,0,a2_4_2*x,0,0,0,0,0,0,0,b1_1U2U4_2,0,b3_2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,k_1_3U4,0,0,0,0,0,0,0,0,a4_3_2*x,a3_4_2*x,0,0,0,0,0,0,0,0,b1_1U3U4_2,b2_2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3,0,0,0,0,0,0,0,0,0,a3_1U2_2*x,a2_1U3_2*x,0,a1_2U3_2*x,0,0,0,0,0,0,b4_1U2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U4,0,0,0,0,0,0,0,0,a4_1U2_2*x,0,a2_1U4_2*x,0,a1_2U4_2*x,0,0,0,0,0,b3_1U2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U4,0,0,0,0,0,0,0,0,a4_1U3_2*x,a3_1U4_2*x,0,0,a1_3U4_2*x,0,0,0,0,b2_1U2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U4,0,0,0,0,0,0,0,0,0,a4_2U3_2*x,a3_2U4_2*x,a2_3U4_2*x,0,0,0,0,b1_1U2U3U4_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_2*x,a3_1U2U4_2*x,a2_1U3U4_2*x,a1_2U3U4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1_3,b2_2_3,b3_3_3,b4_4_3,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a1_0_3*x,0,0,0,0,b2_1U2_3,b3_1U3_3,b4_1U4_3,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_0_3*x,0,0,0,0,b1_1U2_3,0,0,b3_2U3_3,b4_2U4_3,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3,0,0,0,0,0,0,0,0,0,0,0,0,a3_0_3*x,0,0,0,0,0,b1_1U3_3,0,b2_2U3_3,0,b4_3U4_3,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_4,0,0,0,0,0,0,0,0,0,0,0,a4_0_3*x,0,0,0,0,0,0,b1_1U4_3,0,b2_2U4_3,b3_3U4_3,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2,0,0,0,0,0,0,0,0,0,0,0,a2_1_3*x,a1_2_3*x,0,0,0,0,0,0,0,0,b3_1U2U3_3,b4_1U2U4_3,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3,0,0,0,0,0,0,0,0,0,0,a3_1_3*x,0,a1_3_3*x,0,0,0,0,0,0,0,b2_1U2U3_3,0,b4_1U3U4_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U4,0,0,0,0,0,0,0,0,0,a4_1_3*x,0,0,a1_4_3*x,0,0,0,0,0,0,0,b2_1U2U4_3,b3_1U3U4_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3,0,0,0,0,0,0,0,0,0,a3_2_3*x,a2_3_3*x,0,0,0,0,0,0,0,b1_1U2U3_3,0,0,b4_2U3U4_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U4,0,0,0,0,0,0,0,0,a4_2_3*x,0,a2_4_3*x,0,0,0,0,0,0,0,b1_1U2U4_3,0,b3_2U3U4_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U4,0,0,0,0,0,0,0,0,a4_3_3*x,a3_4_3*x,0,0,0,0,0,0,0,0,b1_1U3U4_3,b2_2U3U4_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3,0,0,0,0,0,0,0,0,0,a3_1U2_3*x,a2_1U3_3*x,0,a1_2U3_3*x,0,0,0,0,0,0,b4_1U2U3U4_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U4,0,0,0,0,0,0,0,0,a4_1U2_3*x,0,a2_1U4_3*x,0,a1_2U4_3*x,0,0,0,0,0,b3_1U2U3U4_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U4,0,0,0,0,0,0,0,0,a4_1U3_3*x,a3_1U4_3*x,0,0,a1_3U4_3*x,0,0,0,0,b2_1U2U3U4_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U4,0,0,0,0,0,0,0,0,0,a4_2U3_3*x,a3_2U4_3*x,a2_3U4_3*x,0,0,0,0,b1_1U2U3U4_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_3*x,a3_1U2U4_3*x,a2_1U3U4_3*x,a1_2U3U4_3*x,0;

    for (unsigned i = 0; i < 48; ++i){

        L(i, i) = -(L.col(i).sum());
    }
    
    Matrix<InternalType, Dynamic, 1> steady_state;
    try
    {
    steady_state = getOneDimNullspaceFromSVD<InternalType>(L);
    }
    catch (const std::runtime_error& e)
    {
        throw;
    }
    InternalType norm = steady_state.sum();
    for (int i=0; i<steady_state.size();i++){
        steady_state[i]=steady_state[i] / norm;
    }
    return steady_state;}
    
py::array_t<double> interfacerhos(py::array_t<double> parsar, py::array_t<double> varvals ) {
    Matrix<InternalType, Dynamic, 1> rhos;
    rhos=getrhos(parsar,varvals);
    py::array_t<double> resultpy = py::array_t<double>(48);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    
    for (int i=0;i<48;i++){
        ptrresultpy[i]=rhos[i].template convert_to<double>();
    }
    return resultpy;}
double interfacess(py::array_t<double> parsar, py::array_t<double> varvals){
    std::vector<int> indicesC={32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
    std::vector<int> coeffsC={6,13,23,33,43,53,63,73,89,99,121,149,159,175,203,237};

    Matrix<InternalType, Dynamic, 1> rhos;
    rhos=getrhos(parsar,varvals);
    double ss=0;
    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    for (int i=0;i<indicesC.size();i++){
            ss+=(rhos[indicesC[i]]*pars[coeffsC[i]]).template convert_to<double>();
        }
    return ss;
    }
    PYBIND11_MODULE(base_4sites_highprec,m){

    m.def("interfacerhos", &interfacerhos, "A function which returns the normalised ss rhos.",
    py::arg("parsar"), py::arg("varvals"));
    m.def("interfacess", &interfacess, "A function which returns the ss output, where appropriate rhos are multiplied by rates.",
    py::arg("parsar"), py::arg("varvals"));
}

    