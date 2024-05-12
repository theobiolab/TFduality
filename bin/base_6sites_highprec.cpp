
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
    long double a5_0_1=pars[44];
    long double a5_0_2=pars[45];
    long double a5_0_3=pars[46];
    long double b5_5_1=pars[47];
    long double k_1_5=pars[48];
    long double kr_1_5=pars[49];
    long double b5_5_2=pars[50];
    long double k_2_5=pars[51];
    long double b5_5_3=pars[52];
    long double k_3_5=pars[53];
    long double a6_0_1=pars[54];
    long double a6_0_2=pars[55];
    long double a6_0_3=pars[56];
    long double b6_6_1=pars[57];
    long double k_1_6=pars[58];
    long double kr_1_6=pars[59];
    long double b6_6_2=pars[60];
    long double k_2_6=pars[61];
    long double b6_6_3=pars[62];
    long double k_3_6=pars[63];
    long double a2_1_1=pars[64];
    long double a2_1_2=pars[65];
    long double a2_1_3=pars[66];
    long double b2_1U2_1=pars[67];
    long double k_1_1U2=pars[68];
    long double kr_1_1U2=pars[69];
    long double b2_1U2_2=pars[70];
    long double k_2_1U2=pars[71];
    long double b2_1U2_3=pars[72];
    long double k_3_1U2=pars[73];
    long double a3_1_1=pars[74];
    long double a3_1_2=pars[75];
    long double a3_1_3=pars[76];
    long double b3_1U3_1=pars[77];
    long double k_1_1U3=pars[78];
    long double kr_1_1U3=pars[79];
    long double b3_1U3_2=pars[80];
    long double k_2_1U3=pars[81];
    long double b3_1U3_3=pars[82];
    long double k_3_1U3=pars[83];
    long double a4_1_1=pars[84];
    long double a4_1_2=pars[85];
    long double a4_1_3=pars[86];
    long double b4_1U4_1=pars[87];
    long double k_1_1U4=pars[88];
    long double kr_1_1U4=pars[89];
    long double b4_1U4_2=pars[90];
    long double k_2_1U4=pars[91];
    long double b4_1U4_3=pars[92];
    long double k_3_1U4=pars[93];
    long double a5_1_1=pars[94];
    long double a5_1_2=pars[95];
    long double a5_1_3=pars[96];
    long double b5_1U5_1=pars[97];
    long double k_1_1U5=pars[98];
    long double kr_1_1U5=pars[99];
    long double b5_1U5_2=pars[100];
    long double k_2_1U5=pars[101];
    long double b5_1U5_3=pars[102];
    long double k_3_1U5=pars[103];
    long double a6_1_1=pars[104];
    long double a6_1_2=pars[105];
    long double a6_1_3=pars[106];
    long double b6_1U6_1=pars[107];
    long double k_1_1U6=pars[108];
    long double kr_1_1U6=pars[109];
    long double b6_1U6_2=pars[110];
    long double k_2_1U6=pars[111];
    long double b6_1U6_3=pars[112];
    long double k_3_1U6=pars[113];
    long double a1_2_1=pars[114];
    long double a1_2_2=pars[115];
    long double a1_2_3=pars[116];
    long double b1_1U2_1=pars[117];
    long double b1_1U2_2=pars[118];
    long double b1_1U2_3=pars[119];
    long double a3_2_1=pars[120];
    long double a3_2_2=pars[121];
    long double a3_2_3=pars[122];
    long double b3_2U3_1=pars[123];
    long double k_1_2U3=pars[124];
    long double kr_1_2U3=pars[125];
    long double b3_2U3_2=pars[126];
    long double k_2_2U3=pars[127];
    long double b3_2U3_3=pars[128];
    long double k_3_2U3=pars[129];
    long double a4_2_1=pars[130];
    long double a4_2_2=pars[131];
    long double a4_2_3=pars[132];
    long double b4_2U4_1=pars[133];
    long double k_1_2U4=pars[134];
    long double kr_1_2U4=pars[135];
    long double b4_2U4_2=pars[136];
    long double k_2_2U4=pars[137];
    long double b4_2U4_3=pars[138];
    long double k_3_2U4=pars[139];
    long double a5_2_1=pars[140];
    long double a5_2_2=pars[141];
    long double a5_2_3=pars[142];
    long double b5_2U5_1=pars[143];
    long double k_1_2U5=pars[144];
    long double kr_1_2U5=pars[145];
    long double b5_2U5_2=pars[146];
    long double k_2_2U5=pars[147];
    long double b5_2U5_3=pars[148];
    long double k_3_2U5=pars[149];
    long double a6_2_1=pars[150];
    long double a6_2_2=pars[151];
    long double a6_2_3=pars[152];
    long double b6_2U6_1=pars[153];
    long double k_1_2U6=pars[154];
    long double kr_1_2U6=pars[155];
    long double b6_2U6_2=pars[156];
    long double k_2_2U6=pars[157];
    long double b6_2U6_3=pars[158];
    long double k_3_2U6=pars[159];
    long double a1_3_1=pars[160];
    long double a1_3_2=pars[161];
    long double a1_3_3=pars[162];
    long double b1_1U3_1=pars[163];
    long double b1_1U3_2=pars[164];
    long double b1_1U3_3=pars[165];
    long double a2_3_1=pars[166];
    long double a2_3_2=pars[167];
    long double a2_3_3=pars[168];
    long double b2_2U3_1=pars[169];
    long double b2_2U3_2=pars[170];
    long double b2_2U3_3=pars[171];
    long double a4_3_1=pars[172];
    long double a4_3_2=pars[173];
    long double a4_3_3=pars[174];
    long double b4_3U4_1=pars[175];
    long double k_1_3U4=pars[176];
    long double kr_1_3U4=pars[177];
    long double b4_3U4_2=pars[178];
    long double k_2_3U4=pars[179];
    long double b4_3U4_3=pars[180];
    long double k_3_3U4=pars[181];
    long double a5_3_1=pars[182];
    long double a5_3_2=pars[183];
    long double a5_3_3=pars[184];
    long double b5_3U5_1=pars[185];
    long double k_1_3U5=pars[186];
    long double kr_1_3U5=pars[187];
    long double b5_3U5_2=pars[188];
    long double k_2_3U5=pars[189];
    long double b5_3U5_3=pars[190];
    long double k_3_3U5=pars[191];
    long double a6_3_1=pars[192];
    long double a6_3_2=pars[193];
    long double a6_3_3=pars[194];
    long double b6_3U6_1=pars[195];
    long double k_1_3U6=pars[196];
    long double kr_1_3U6=pars[197];
    long double b6_3U6_2=pars[198];
    long double k_2_3U6=pars[199];
    long double b6_3U6_3=pars[200];
    long double k_3_3U6=pars[201];
    long double a1_4_1=pars[202];
    long double a1_4_2=pars[203];
    long double a1_4_3=pars[204];
    long double b1_1U4_1=pars[205];
    long double b1_1U4_2=pars[206];
    long double b1_1U4_3=pars[207];
    long double a2_4_1=pars[208];
    long double a2_4_2=pars[209];
    long double a2_4_3=pars[210];
    long double b2_2U4_1=pars[211];
    long double b2_2U4_2=pars[212];
    long double b2_2U4_3=pars[213];
    long double a3_4_1=pars[214];
    long double a3_4_2=pars[215];
    long double a3_4_3=pars[216];
    long double b3_3U4_1=pars[217];
    long double b3_3U4_2=pars[218];
    long double b3_3U4_3=pars[219];
    long double a5_4_1=pars[220];
    long double a5_4_2=pars[221];
    long double a5_4_3=pars[222];
    long double b5_4U5_1=pars[223];
    long double k_1_4U5=pars[224];
    long double kr_1_4U5=pars[225];
    long double b5_4U5_2=pars[226];
    long double k_2_4U5=pars[227];
    long double b5_4U5_3=pars[228];
    long double k_3_4U5=pars[229];
    long double a6_4_1=pars[230];
    long double a6_4_2=pars[231];
    long double a6_4_3=pars[232];
    long double b6_4U6_1=pars[233];
    long double k_1_4U6=pars[234];
    long double kr_1_4U6=pars[235];
    long double b6_4U6_2=pars[236];
    long double k_2_4U6=pars[237];
    long double b6_4U6_3=pars[238];
    long double k_3_4U6=pars[239];
    long double a1_5_1=pars[240];
    long double a1_5_2=pars[241];
    long double a1_5_3=pars[242];
    long double b1_1U5_1=pars[243];
    long double b1_1U5_2=pars[244];
    long double b1_1U5_3=pars[245];
    long double a2_5_1=pars[246];
    long double a2_5_2=pars[247];
    long double a2_5_3=pars[248];
    long double b2_2U5_1=pars[249];
    long double b2_2U5_2=pars[250];
    long double b2_2U5_3=pars[251];
    long double a3_5_1=pars[252];
    long double a3_5_2=pars[253];
    long double a3_5_3=pars[254];
    long double b3_3U5_1=pars[255];
    long double b3_3U5_2=pars[256];
    long double b3_3U5_3=pars[257];
    long double a4_5_1=pars[258];
    long double a4_5_2=pars[259];
    long double a4_5_3=pars[260];
    long double b4_4U5_1=pars[261];
    long double b4_4U5_2=pars[262];
    long double b4_4U5_3=pars[263];
    long double a6_5_1=pars[264];
    long double a6_5_2=pars[265];
    long double a6_5_3=pars[266];
    long double b6_5U6_1=pars[267];
    long double k_1_5U6=pars[268];
    long double kr_1_5U6=pars[269];
    long double b6_5U6_2=pars[270];
    long double k_2_5U6=pars[271];
    long double b6_5U6_3=pars[272];
    long double k_3_5U6=pars[273];
    long double a1_6_1=pars[274];
    long double a1_6_2=pars[275];
    long double a1_6_3=pars[276];
    long double b1_1U6_1=pars[277];
    long double b1_1U6_2=pars[278];
    long double b1_1U6_3=pars[279];
    long double a2_6_1=pars[280];
    long double a2_6_2=pars[281];
    long double a2_6_3=pars[282];
    long double b2_2U6_1=pars[283];
    long double b2_2U6_2=pars[284];
    long double b2_2U6_3=pars[285];
    long double a3_6_1=pars[286];
    long double a3_6_2=pars[287];
    long double a3_6_3=pars[288];
    long double b3_3U6_1=pars[289];
    long double b3_3U6_2=pars[290];
    long double b3_3U6_3=pars[291];
    long double a4_6_1=pars[292];
    long double a4_6_2=pars[293];
    long double a4_6_3=pars[294];
    long double b4_4U6_1=pars[295];
    long double b4_4U6_2=pars[296];
    long double b4_4U6_3=pars[297];
    long double a5_6_1=pars[298];
    long double a5_6_2=pars[299];
    long double a5_6_3=pars[300];
    long double b5_5U6_1=pars[301];
    long double b5_5U6_2=pars[302];
    long double b5_5U6_3=pars[303];
    long double a3_1U2_1=pars[304];
    long double a3_1U2_2=pars[305];
    long double a3_1U2_3=pars[306];
    long double b3_1U2U3_1=pars[307];
    long double k_1_1U2U3=pars[308];
    long double kr_1_1U2U3=pars[309];
    long double b3_1U2U3_2=pars[310];
    long double k_2_1U2U3=pars[311];
    long double b3_1U2U3_3=pars[312];
    long double k_3_1U2U3=pars[313];
    long double a4_1U2_1=pars[314];
    long double a4_1U2_2=pars[315];
    long double a4_1U2_3=pars[316];
    long double b4_1U2U4_1=pars[317];
    long double k_1_1U2U4=pars[318];
    long double kr_1_1U2U4=pars[319];
    long double b4_1U2U4_2=pars[320];
    long double k_2_1U2U4=pars[321];
    long double b4_1U2U4_3=pars[322];
    long double k_3_1U2U4=pars[323];
    long double a5_1U2_1=pars[324];
    long double a5_1U2_2=pars[325];
    long double a5_1U2_3=pars[326];
    long double b5_1U2U5_1=pars[327];
    long double k_1_1U2U5=pars[328];
    long double kr_1_1U2U5=pars[329];
    long double b5_1U2U5_2=pars[330];
    long double k_2_1U2U5=pars[331];
    long double b5_1U2U5_3=pars[332];
    long double k_3_1U2U5=pars[333];
    long double a6_1U2_1=pars[334];
    long double a6_1U2_2=pars[335];
    long double a6_1U2_3=pars[336];
    long double b6_1U2U6_1=pars[337];
    long double k_1_1U2U6=pars[338];
    long double kr_1_1U2U6=pars[339];
    long double b6_1U2U6_2=pars[340];
    long double k_2_1U2U6=pars[341];
    long double b6_1U2U6_3=pars[342];
    long double k_3_1U2U6=pars[343];
    long double a2_1U3_1=pars[344];
    long double a2_1U3_2=pars[345];
    long double a2_1U3_3=pars[346];
    long double b2_1U2U3_1=pars[347];
    long double b2_1U2U3_2=pars[348];
    long double b2_1U2U3_3=pars[349];
    long double a4_1U3_1=pars[350];
    long double a4_1U3_2=pars[351];
    long double a4_1U3_3=pars[352];
    long double b4_1U3U4_1=pars[353];
    long double k_1_1U3U4=pars[354];
    long double kr_1_1U3U4=pars[355];
    long double b4_1U3U4_2=pars[356];
    long double k_2_1U3U4=pars[357];
    long double b4_1U3U4_3=pars[358];
    long double k_3_1U3U4=pars[359];
    long double a5_1U3_1=pars[360];
    long double a5_1U3_2=pars[361];
    long double a5_1U3_3=pars[362];
    long double b5_1U3U5_1=pars[363];
    long double k_1_1U3U5=pars[364];
    long double kr_1_1U3U5=pars[365];
    long double b5_1U3U5_2=pars[366];
    long double k_2_1U3U5=pars[367];
    long double b5_1U3U5_3=pars[368];
    long double k_3_1U3U5=pars[369];
    long double a6_1U3_1=pars[370];
    long double a6_1U3_2=pars[371];
    long double a6_1U3_3=pars[372];
    long double b6_1U3U6_1=pars[373];
    long double k_1_1U3U6=pars[374];
    long double kr_1_1U3U6=pars[375];
    long double b6_1U3U6_2=pars[376];
    long double k_2_1U3U6=pars[377];
    long double b6_1U3U6_3=pars[378];
    long double k_3_1U3U6=pars[379];
    long double a2_1U4_1=pars[380];
    long double a2_1U4_2=pars[381];
    long double a2_1U4_3=pars[382];
    long double b2_1U2U4_1=pars[383];
    long double b2_1U2U4_2=pars[384];
    long double b2_1U2U4_3=pars[385];
    long double a3_1U4_1=pars[386];
    long double a3_1U4_2=pars[387];
    long double a3_1U4_3=pars[388];
    long double b3_1U3U4_1=pars[389];
    long double b3_1U3U4_2=pars[390];
    long double b3_1U3U4_3=pars[391];
    long double a5_1U4_1=pars[392];
    long double a5_1U4_2=pars[393];
    long double a5_1U4_3=pars[394];
    long double b5_1U4U5_1=pars[395];
    long double k_1_1U4U5=pars[396];
    long double kr_1_1U4U5=pars[397];
    long double b5_1U4U5_2=pars[398];
    long double k_2_1U4U5=pars[399];
    long double b5_1U4U5_3=pars[400];
    long double k_3_1U4U5=pars[401];
    long double a6_1U4_1=pars[402];
    long double a6_1U4_2=pars[403];
    long double a6_1U4_3=pars[404];
    long double b6_1U4U6_1=pars[405];
    long double k_1_1U4U6=pars[406];
    long double kr_1_1U4U6=pars[407];
    long double b6_1U4U6_2=pars[408];
    long double k_2_1U4U6=pars[409];
    long double b6_1U4U6_3=pars[410];
    long double k_3_1U4U6=pars[411];
    long double a2_1U5_1=pars[412];
    long double a2_1U5_2=pars[413];
    long double a2_1U5_3=pars[414];
    long double b2_1U2U5_1=pars[415];
    long double b2_1U2U5_2=pars[416];
    long double b2_1U2U5_3=pars[417];
    long double a3_1U5_1=pars[418];
    long double a3_1U5_2=pars[419];
    long double a3_1U5_3=pars[420];
    long double b3_1U3U5_1=pars[421];
    long double b3_1U3U5_2=pars[422];
    long double b3_1U3U5_3=pars[423];
    long double a4_1U5_1=pars[424];
    long double a4_1U5_2=pars[425];
    long double a4_1U5_3=pars[426];
    long double b4_1U4U5_1=pars[427];
    long double b4_1U4U5_2=pars[428];
    long double b4_1U4U5_3=pars[429];
    long double a6_1U5_1=pars[430];
    long double a6_1U5_2=pars[431];
    long double a6_1U5_3=pars[432];
    long double b6_1U5U6_1=pars[433];
    long double k_1_1U5U6=pars[434];
    long double kr_1_1U5U6=pars[435];
    long double b6_1U5U6_2=pars[436];
    long double k_2_1U5U6=pars[437];
    long double b6_1U5U6_3=pars[438];
    long double k_3_1U5U6=pars[439];
    long double a2_1U6_1=pars[440];
    long double a2_1U6_2=pars[441];
    long double a2_1U6_3=pars[442];
    long double b2_1U2U6_1=pars[443];
    long double b2_1U2U6_2=pars[444];
    long double b2_1U2U6_3=pars[445];
    long double a3_1U6_1=pars[446];
    long double a3_1U6_2=pars[447];
    long double a3_1U6_3=pars[448];
    long double b3_1U3U6_1=pars[449];
    long double b3_1U3U6_2=pars[450];
    long double b3_1U3U6_3=pars[451];
    long double a4_1U6_1=pars[452];
    long double a4_1U6_2=pars[453];
    long double a4_1U6_3=pars[454];
    long double b4_1U4U6_1=pars[455];
    long double b4_1U4U6_2=pars[456];
    long double b4_1U4U6_3=pars[457];
    long double a5_1U6_1=pars[458];
    long double a5_1U6_2=pars[459];
    long double a5_1U6_3=pars[460];
    long double b5_1U5U6_1=pars[461];
    long double b5_1U5U6_2=pars[462];
    long double b5_1U5U6_3=pars[463];
    long double a1_2U3_1=pars[464];
    long double a1_2U3_2=pars[465];
    long double a1_2U3_3=pars[466];
    long double b1_1U2U3_1=pars[467];
    long double b1_1U2U3_2=pars[468];
    long double b1_1U2U3_3=pars[469];
    long double a4_2U3_1=pars[470];
    long double a4_2U3_2=pars[471];
    long double a4_2U3_3=pars[472];
    long double b4_2U3U4_1=pars[473];
    long double k_1_2U3U4=pars[474];
    long double kr_1_2U3U4=pars[475];
    long double b4_2U3U4_2=pars[476];
    long double k_2_2U3U4=pars[477];
    long double b4_2U3U4_3=pars[478];
    long double k_3_2U3U4=pars[479];
    long double a5_2U3_1=pars[480];
    long double a5_2U3_2=pars[481];
    long double a5_2U3_3=pars[482];
    long double b5_2U3U5_1=pars[483];
    long double k_1_2U3U5=pars[484];
    long double kr_1_2U3U5=pars[485];
    long double b5_2U3U5_2=pars[486];
    long double k_2_2U3U5=pars[487];
    long double b5_2U3U5_3=pars[488];
    long double k_3_2U3U5=pars[489];
    long double a6_2U3_1=pars[490];
    long double a6_2U3_2=pars[491];
    long double a6_2U3_3=pars[492];
    long double b6_2U3U6_1=pars[493];
    long double k_1_2U3U6=pars[494];
    long double kr_1_2U3U6=pars[495];
    long double b6_2U3U6_2=pars[496];
    long double k_2_2U3U6=pars[497];
    long double b6_2U3U6_3=pars[498];
    long double k_3_2U3U6=pars[499];
    long double a1_2U4_1=pars[500];
    long double a1_2U4_2=pars[501];
    long double a1_2U4_3=pars[502];
    long double b1_1U2U4_1=pars[503];
    long double b1_1U2U4_2=pars[504];
    long double b1_1U2U4_3=pars[505];
    long double a3_2U4_1=pars[506];
    long double a3_2U4_2=pars[507];
    long double a3_2U4_3=pars[508];
    long double b3_2U3U4_1=pars[509];
    long double b3_2U3U4_2=pars[510];
    long double b3_2U3U4_3=pars[511];
    long double a5_2U4_1=pars[512];
    long double a5_2U4_2=pars[513];
    long double a5_2U4_3=pars[514];
    long double b5_2U4U5_1=pars[515];
    long double k_1_2U4U5=pars[516];
    long double kr_1_2U4U5=pars[517];
    long double b5_2U4U5_2=pars[518];
    long double k_2_2U4U5=pars[519];
    long double b5_2U4U5_3=pars[520];
    long double k_3_2U4U5=pars[521];
    long double a6_2U4_1=pars[522];
    long double a6_2U4_2=pars[523];
    long double a6_2U4_3=pars[524];
    long double b6_2U4U6_1=pars[525];
    long double k_1_2U4U6=pars[526];
    long double kr_1_2U4U6=pars[527];
    long double b6_2U4U6_2=pars[528];
    long double k_2_2U4U6=pars[529];
    long double b6_2U4U6_3=pars[530];
    long double k_3_2U4U6=pars[531];
    long double a1_2U5_1=pars[532];
    long double a1_2U5_2=pars[533];
    long double a1_2U5_3=pars[534];
    long double b1_1U2U5_1=pars[535];
    long double b1_1U2U5_2=pars[536];
    long double b1_1U2U5_3=pars[537];
    long double a3_2U5_1=pars[538];
    long double a3_2U5_2=pars[539];
    long double a3_2U5_3=pars[540];
    long double b3_2U3U5_1=pars[541];
    long double b3_2U3U5_2=pars[542];
    long double b3_2U3U5_3=pars[543];
    long double a4_2U5_1=pars[544];
    long double a4_2U5_2=pars[545];
    long double a4_2U5_3=pars[546];
    long double b4_2U4U5_1=pars[547];
    long double b4_2U4U5_2=pars[548];
    long double b4_2U4U5_3=pars[549];
    long double a6_2U5_1=pars[550];
    long double a6_2U5_2=pars[551];
    long double a6_2U5_3=pars[552];
    long double b6_2U5U6_1=pars[553];
    long double k_1_2U5U6=pars[554];
    long double kr_1_2U5U6=pars[555];
    long double b6_2U5U6_2=pars[556];
    long double k_2_2U5U6=pars[557];
    long double b6_2U5U6_3=pars[558];
    long double k_3_2U5U6=pars[559];
    long double a1_2U6_1=pars[560];
    long double a1_2U6_2=pars[561];
    long double a1_2U6_3=pars[562];
    long double b1_1U2U6_1=pars[563];
    long double b1_1U2U6_2=pars[564];
    long double b1_1U2U6_3=pars[565];
    long double a3_2U6_1=pars[566];
    long double a3_2U6_2=pars[567];
    long double a3_2U6_3=pars[568];
    long double b3_2U3U6_1=pars[569];
    long double b3_2U3U6_2=pars[570];
    long double b3_2U3U6_3=pars[571];
    long double a4_2U6_1=pars[572];
    long double a4_2U6_2=pars[573];
    long double a4_2U6_3=pars[574];
    long double b4_2U4U6_1=pars[575];
    long double b4_2U4U6_2=pars[576];
    long double b4_2U4U6_3=pars[577];
    long double a5_2U6_1=pars[578];
    long double a5_2U6_2=pars[579];
    long double a5_2U6_3=pars[580];
    long double b5_2U5U6_1=pars[581];
    long double b5_2U5U6_2=pars[582];
    long double b5_2U5U6_3=pars[583];
    long double a1_3U4_1=pars[584];
    long double a1_3U4_2=pars[585];
    long double a1_3U4_3=pars[586];
    long double b1_1U3U4_1=pars[587];
    long double b1_1U3U4_2=pars[588];
    long double b1_1U3U4_3=pars[589];
    long double a2_3U4_1=pars[590];
    long double a2_3U4_2=pars[591];
    long double a2_3U4_3=pars[592];
    long double b2_2U3U4_1=pars[593];
    long double b2_2U3U4_2=pars[594];
    long double b2_2U3U4_3=pars[595];
    long double a5_3U4_1=pars[596];
    long double a5_3U4_2=pars[597];
    long double a5_3U4_3=pars[598];
    long double b5_3U4U5_1=pars[599];
    long double k_1_3U4U5=pars[600];
    long double kr_1_3U4U5=pars[601];
    long double b5_3U4U5_2=pars[602];
    long double k_2_3U4U5=pars[603];
    long double b5_3U4U5_3=pars[604];
    long double k_3_3U4U5=pars[605];
    long double a6_3U4_1=pars[606];
    long double a6_3U4_2=pars[607];
    long double a6_3U4_3=pars[608];
    long double b6_3U4U6_1=pars[609];
    long double k_1_3U4U6=pars[610];
    long double kr_1_3U4U6=pars[611];
    long double b6_3U4U6_2=pars[612];
    long double k_2_3U4U6=pars[613];
    long double b6_3U4U6_3=pars[614];
    long double k_3_3U4U6=pars[615];
    long double a1_3U5_1=pars[616];
    long double a1_3U5_2=pars[617];
    long double a1_3U5_3=pars[618];
    long double b1_1U3U5_1=pars[619];
    long double b1_1U3U5_2=pars[620];
    long double b1_1U3U5_3=pars[621];
    long double a2_3U5_1=pars[622];
    long double a2_3U5_2=pars[623];
    long double a2_3U5_3=pars[624];
    long double b2_2U3U5_1=pars[625];
    long double b2_2U3U5_2=pars[626];
    long double b2_2U3U5_3=pars[627];
    long double a4_3U5_1=pars[628];
    long double a4_3U5_2=pars[629];
    long double a4_3U5_3=pars[630];
    long double b4_3U4U5_1=pars[631];
    long double b4_3U4U5_2=pars[632];
    long double b4_3U4U5_3=pars[633];
    long double a6_3U5_1=pars[634];
    long double a6_3U5_2=pars[635];
    long double a6_3U5_3=pars[636];
    long double b6_3U5U6_1=pars[637];
    long double k_1_3U5U6=pars[638];
    long double kr_1_3U5U6=pars[639];
    long double b6_3U5U6_2=pars[640];
    long double k_2_3U5U6=pars[641];
    long double b6_3U5U6_3=pars[642];
    long double k_3_3U5U6=pars[643];
    long double a1_3U6_1=pars[644];
    long double a1_3U6_2=pars[645];
    long double a1_3U6_3=pars[646];
    long double b1_1U3U6_1=pars[647];
    long double b1_1U3U6_2=pars[648];
    long double b1_1U3U6_3=pars[649];
    long double a2_3U6_1=pars[650];
    long double a2_3U6_2=pars[651];
    long double a2_3U6_3=pars[652];
    long double b2_2U3U6_1=pars[653];
    long double b2_2U3U6_2=pars[654];
    long double b2_2U3U6_3=pars[655];
    long double a4_3U6_1=pars[656];
    long double a4_3U6_2=pars[657];
    long double a4_3U6_3=pars[658];
    long double b4_3U4U6_1=pars[659];
    long double b4_3U4U6_2=pars[660];
    long double b4_3U4U6_3=pars[661];
    long double a5_3U6_1=pars[662];
    long double a5_3U6_2=pars[663];
    long double a5_3U6_3=pars[664];
    long double b5_3U5U6_1=pars[665];
    long double b5_3U5U6_2=pars[666];
    long double b5_3U5U6_3=pars[667];
    long double a1_4U5_1=pars[668];
    long double a1_4U5_2=pars[669];
    long double a1_4U5_3=pars[670];
    long double b1_1U4U5_1=pars[671];
    long double b1_1U4U5_2=pars[672];
    long double b1_1U4U5_3=pars[673];
    long double a2_4U5_1=pars[674];
    long double a2_4U5_2=pars[675];
    long double a2_4U5_3=pars[676];
    long double b2_2U4U5_1=pars[677];
    long double b2_2U4U5_2=pars[678];
    long double b2_2U4U5_3=pars[679];
    long double a3_4U5_1=pars[680];
    long double a3_4U5_2=pars[681];
    long double a3_4U5_3=pars[682];
    long double b3_3U4U5_1=pars[683];
    long double b3_3U4U5_2=pars[684];
    long double b3_3U4U5_3=pars[685];
    long double a6_4U5_1=pars[686];
    long double a6_4U5_2=pars[687];
    long double a6_4U5_3=pars[688];
    long double b6_4U5U6_1=pars[689];
    long double k_1_4U5U6=pars[690];
    long double kr_1_4U5U6=pars[691];
    long double b6_4U5U6_2=pars[692];
    long double k_2_4U5U6=pars[693];
    long double b6_4U5U6_3=pars[694];
    long double k_3_4U5U6=pars[695];
    long double a1_4U6_1=pars[696];
    long double a1_4U6_2=pars[697];
    long double a1_4U6_3=pars[698];
    long double b1_1U4U6_1=pars[699];
    long double b1_1U4U6_2=pars[700];
    long double b1_1U4U6_3=pars[701];
    long double a2_4U6_1=pars[702];
    long double a2_4U6_2=pars[703];
    long double a2_4U6_3=pars[704];
    long double b2_2U4U6_1=pars[705];
    long double b2_2U4U6_2=pars[706];
    long double b2_2U4U6_3=pars[707];
    long double a3_4U6_1=pars[708];
    long double a3_4U6_2=pars[709];
    long double a3_4U6_3=pars[710];
    long double b3_3U4U6_1=pars[711];
    long double b3_3U4U6_2=pars[712];
    long double b3_3U4U6_3=pars[713];
    long double a5_4U6_1=pars[714];
    long double a5_4U6_2=pars[715];
    long double a5_4U6_3=pars[716];
    long double b5_4U5U6_1=pars[717];
    long double b5_4U5U6_2=pars[718];
    long double b5_4U5U6_3=pars[719];
    long double a1_5U6_1=pars[720];
    long double a1_5U6_2=pars[721];
    long double a1_5U6_3=pars[722];
    long double b1_1U5U6_1=pars[723];
    long double b1_1U5U6_2=pars[724];
    long double b1_1U5U6_3=pars[725];
    long double a2_5U6_1=pars[726];
    long double a2_5U6_2=pars[727];
    long double a2_5U6_3=pars[728];
    long double b2_2U5U6_1=pars[729];
    long double b2_2U5U6_2=pars[730];
    long double b2_2U5U6_3=pars[731];
    long double a3_5U6_1=pars[732];
    long double a3_5U6_2=pars[733];
    long double a3_5U6_3=pars[734];
    long double b3_3U5U6_1=pars[735];
    long double b3_3U5U6_2=pars[736];
    long double b3_3U5U6_3=pars[737];
    long double a4_5U6_1=pars[738];
    long double a4_5U6_2=pars[739];
    long double a4_5U6_3=pars[740];
    long double b4_4U5U6_1=pars[741];
    long double b4_4U5U6_2=pars[742];
    long double b4_4U5U6_3=pars[743];
    long double a4_1U2U3_1=pars[744];
    long double a4_1U2U3_2=pars[745];
    long double a4_1U2U3_3=pars[746];
    long double b4_1U2U3U4_1=pars[747];
    long double k_1_1U2U3U4=pars[748];
    long double kr_1_1U2U3U4=pars[749];
    long double b4_1U2U3U4_2=pars[750];
    long double k_2_1U2U3U4=pars[751];
    long double b4_1U2U3U4_3=pars[752];
    long double k_3_1U2U3U4=pars[753];
    long double a5_1U2U3_1=pars[754];
    long double a5_1U2U3_2=pars[755];
    long double a5_1U2U3_3=pars[756];
    long double b5_1U2U3U5_1=pars[757];
    long double k_1_1U2U3U5=pars[758];
    long double kr_1_1U2U3U5=pars[759];
    long double b5_1U2U3U5_2=pars[760];
    long double k_2_1U2U3U5=pars[761];
    long double b5_1U2U3U5_3=pars[762];
    long double k_3_1U2U3U5=pars[763];
    long double a6_1U2U3_1=pars[764];
    long double a6_1U2U3_2=pars[765];
    long double a6_1U2U3_3=pars[766];
    long double b6_1U2U3U6_1=pars[767];
    long double k_1_1U2U3U6=pars[768];
    long double kr_1_1U2U3U6=pars[769];
    long double b6_1U2U3U6_2=pars[770];
    long double k_2_1U2U3U6=pars[771];
    long double b6_1U2U3U6_3=pars[772];
    long double k_3_1U2U3U6=pars[773];
    long double a3_1U2U4_1=pars[774];
    long double a3_1U2U4_2=pars[775];
    long double a3_1U2U4_3=pars[776];
    long double b3_1U2U3U4_1=pars[777];
    long double b3_1U2U3U4_2=pars[778];
    long double b3_1U2U3U4_3=pars[779];
    long double a5_1U2U4_1=pars[780];
    long double a5_1U2U4_2=pars[781];
    long double a5_1U2U4_3=pars[782];
    long double b5_1U2U4U5_1=pars[783];
    long double k_1_1U2U4U5=pars[784];
    long double kr_1_1U2U4U5=pars[785];
    long double b5_1U2U4U5_2=pars[786];
    long double k_2_1U2U4U5=pars[787];
    long double b5_1U2U4U5_3=pars[788];
    long double k_3_1U2U4U5=pars[789];
    long double a6_1U2U4_1=pars[790];
    long double a6_1U2U4_2=pars[791];
    long double a6_1U2U4_3=pars[792];
    long double b6_1U2U4U6_1=pars[793];
    long double k_1_1U2U4U6=pars[794];
    long double kr_1_1U2U4U6=pars[795];
    long double b6_1U2U4U6_2=pars[796];
    long double k_2_1U2U4U6=pars[797];
    long double b6_1U2U4U6_3=pars[798];
    long double k_3_1U2U4U6=pars[799];
    long double a3_1U2U5_1=pars[800];
    long double a3_1U2U5_2=pars[801];
    long double a3_1U2U5_3=pars[802];
    long double b3_1U2U3U5_1=pars[803];
    long double b3_1U2U3U5_2=pars[804];
    long double b3_1U2U3U5_3=pars[805];
    long double a4_1U2U5_1=pars[806];
    long double a4_1U2U5_2=pars[807];
    long double a4_1U2U5_3=pars[808];
    long double b4_1U2U4U5_1=pars[809];
    long double b4_1U2U4U5_2=pars[810];
    long double b4_1U2U4U5_3=pars[811];
    long double a6_1U2U5_1=pars[812];
    long double a6_1U2U5_2=pars[813];
    long double a6_1U2U5_3=pars[814];
    long double b6_1U2U5U6_1=pars[815];
    long double k_1_1U2U5U6=pars[816];
    long double kr_1_1U2U5U6=pars[817];
    long double b6_1U2U5U6_2=pars[818];
    long double k_2_1U2U5U6=pars[819];
    long double b6_1U2U5U6_3=pars[820];
    long double k_3_1U2U5U6=pars[821];
    long double a3_1U2U6_1=pars[822];
    long double a3_1U2U6_2=pars[823];
    long double a3_1U2U6_3=pars[824];
    long double b3_1U2U3U6_1=pars[825];
    long double b3_1U2U3U6_2=pars[826];
    long double b3_1U2U3U6_3=pars[827];
    long double a4_1U2U6_1=pars[828];
    long double a4_1U2U6_2=pars[829];
    long double a4_1U2U6_3=pars[830];
    long double b4_1U2U4U6_1=pars[831];
    long double b4_1U2U4U6_2=pars[832];
    long double b4_1U2U4U6_3=pars[833];
    long double a5_1U2U6_1=pars[834];
    long double a5_1U2U6_2=pars[835];
    long double a5_1U2U6_3=pars[836];
    long double b5_1U2U5U6_1=pars[837];
    long double b5_1U2U5U6_2=pars[838];
    long double b5_1U2U5U6_3=pars[839];
    long double a2_1U3U4_1=pars[840];
    long double a2_1U3U4_2=pars[841];
    long double a2_1U3U4_3=pars[842];
    long double b2_1U2U3U4_1=pars[843];
    long double b2_1U2U3U4_2=pars[844];
    long double b2_1U2U3U4_3=pars[845];
    long double a5_1U3U4_1=pars[846];
    long double a5_1U3U4_2=pars[847];
    long double a5_1U3U4_3=pars[848];
    long double b5_1U3U4U5_1=pars[849];
    long double k_1_1U3U4U5=pars[850];
    long double kr_1_1U3U4U5=pars[851];
    long double b5_1U3U4U5_2=pars[852];
    long double k_2_1U3U4U5=pars[853];
    long double b5_1U3U4U5_3=pars[854];
    long double k_3_1U3U4U5=pars[855];
    long double a6_1U3U4_1=pars[856];
    long double a6_1U3U4_2=pars[857];
    long double a6_1U3U4_3=pars[858];
    long double b6_1U3U4U6_1=pars[859];
    long double k_1_1U3U4U6=pars[860];
    long double kr_1_1U3U4U6=pars[861];
    long double b6_1U3U4U6_2=pars[862];
    long double k_2_1U3U4U6=pars[863];
    long double b6_1U3U4U6_3=pars[864];
    long double k_3_1U3U4U6=pars[865];
    long double a2_1U3U5_1=pars[866];
    long double a2_1U3U5_2=pars[867];
    long double a2_1U3U5_3=pars[868];
    long double b2_1U2U3U5_1=pars[869];
    long double b2_1U2U3U5_2=pars[870];
    long double b2_1U2U3U5_3=pars[871];
    long double a4_1U3U5_1=pars[872];
    long double a4_1U3U5_2=pars[873];
    long double a4_1U3U5_3=pars[874];
    long double b4_1U3U4U5_1=pars[875];
    long double b4_1U3U4U5_2=pars[876];
    long double b4_1U3U4U5_3=pars[877];
    long double a6_1U3U5_1=pars[878];
    long double a6_1U3U5_2=pars[879];
    long double a6_1U3U5_3=pars[880];
    long double b6_1U3U5U6_1=pars[881];
    long double k_1_1U3U5U6=pars[882];
    long double kr_1_1U3U5U6=pars[883];
    long double b6_1U3U5U6_2=pars[884];
    long double k_2_1U3U5U6=pars[885];
    long double b6_1U3U5U6_3=pars[886];
    long double k_3_1U3U5U6=pars[887];
    long double a2_1U3U6_1=pars[888];
    long double a2_1U3U6_2=pars[889];
    long double a2_1U3U6_3=pars[890];
    long double b2_1U2U3U6_1=pars[891];
    long double b2_1U2U3U6_2=pars[892];
    long double b2_1U2U3U6_3=pars[893];
    long double a4_1U3U6_1=pars[894];
    long double a4_1U3U6_2=pars[895];
    long double a4_1U3U6_3=pars[896];
    long double b4_1U3U4U6_1=pars[897];
    long double b4_1U3U4U6_2=pars[898];
    long double b4_1U3U4U6_3=pars[899];
    long double a5_1U3U6_1=pars[900];
    long double a5_1U3U6_2=pars[901];
    long double a5_1U3U6_3=pars[902];
    long double b5_1U3U5U6_1=pars[903];
    long double b5_1U3U5U6_2=pars[904];
    long double b5_1U3U5U6_3=pars[905];
    long double a2_1U4U5_1=pars[906];
    long double a2_1U4U5_2=pars[907];
    long double a2_1U4U5_3=pars[908];
    long double b2_1U2U4U5_1=pars[909];
    long double b2_1U2U4U5_2=pars[910];
    long double b2_1U2U4U5_3=pars[911];
    long double a3_1U4U5_1=pars[912];
    long double a3_1U4U5_2=pars[913];
    long double a3_1U4U5_3=pars[914];
    long double b3_1U3U4U5_1=pars[915];
    long double b3_1U3U4U5_2=pars[916];
    long double b3_1U3U4U5_3=pars[917];
    long double a6_1U4U5_1=pars[918];
    long double a6_1U4U5_2=pars[919];
    long double a6_1U4U5_3=pars[920];
    long double b6_1U4U5U6_1=pars[921];
    long double k_1_1U4U5U6=pars[922];
    long double kr_1_1U4U5U6=pars[923];
    long double b6_1U4U5U6_2=pars[924];
    long double k_2_1U4U5U6=pars[925];
    long double b6_1U4U5U6_3=pars[926];
    long double k_3_1U4U5U6=pars[927];
    long double a2_1U4U6_1=pars[928];
    long double a2_1U4U6_2=pars[929];
    long double a2_1U4U6_3=pars[930];
    long double b2_1U2U4U6_1=pars[931];
    long double b2_1U2U4U6_2=pars[932];
    long double b2_1U2U4U6_3=pars[933];
    long double a3_1U4U6_1=pars[934];
    long double a3_1U4U6_2=pars[935];
    long double a3_1U4U6_3=pars[936];
    long double b3_1U3U4U6_1=pars[937];
    long double b3_1U3U4U6_2=pars[938];
    long double b3_1U3U4U6_3=pars[939];
    long double a5_1U4U6_1=pars[940];
    long double a5_1U4U6_2=pars[941];
    long double a5_1U4U6_3=pars[942];
    long double b5_1U4U5U6_1=pars[943];
    long double b5_1U4U5U6_2=pars[944];
    long double b5_1U4U5U6_3=pars[945];
    long double a2_1U5U6_1=pars[946];
    long double a2_1U5U6_2=pars[947];
    long double a2_1U5U6_3=pars[948];
    long double b2_1U2U5U6_1=pars[949];
    long double b2_1U2U5U6_2=pars[950];
    long double b2_1U2U5U6_3=pars[951];
    long double a3_1U5U6_1=pars[952];
    long double a3_1U5U6_2=pars[953];
    long double a3_1U5U6_3=pars[954];
    long double b3_1U3U5U6_1=pars[955];
    long double b3_1U3U5U6_2=pars[956];
    long double b3_1U3U5U6_3=pars[957];
    long double a4_1U5U6_1=pars[958];
    long double a4_1U5U6_2=pars[959];
    long double a4_1U5U6_3=pars[960];
    long double b4_1U4U5U6_1=pars[961];
    long double b4_1U4U5U6_2=pars[962];
    long double b4_1U4U5U6_3=pars[963];
    long double a1_2U3U4_1=pars[964];
    long double a1_2U3U4_2=pars[965];
    long double a1_2U3U4_3=pars[966];
    long double b1_1U2U3U4_1=pars[967];
    long double b1_1U2U3U4_2=pars[968];
    long double b1_1U2U3U4_3=pars[969];
    long double a5_2U3U4_1=pars[970];
    long double a5_2U3U4_2=pars[971];
    long double a5_2U3U4_3=pars[972];
    long double b5_2U3U4U5_1=pars[973];
    long double k_1_2U3U4U5=pars[974];
    long double kr_1_2U3U4U5=pars[975];
    long double b5_2U3U4U5_2=pars[976];
    long double k_2_2U3U4U5=pars[977];
    long double b5_2U3U4U5_3=pars[978];
    long double k_3_2U3U4U5=pars[979];
    long double a6_2U3U4_1=pars[980];
    long double a6_2U3U4_2=pars[981];
    long double a6_2U3U4_3=pars[982];
    long double b6_2U3U4U6_1=pars[983];
    long double k_1_2U3U4U6=pars[984];
    long double kr_1_2U3U4U6=pars[985];
    long double b6_2U3U4U6_2=pars[986];
    long double k_2_2U3U4U6=pars[987];
    long double b6_2U3U4U6_3=pars[988];
    long double k_3_2U3U4U6=pars[989];
    long double a1_2U3U5_1=pars[990];
    long double a1_2U3U5_2=pars[991];
    long double a1_2U3U5_3=pars[992];
    long double b1_1U2U3U5_1=pars[993];
    long double b1_1U2U3U5_2=pars[994];
    long double b1_1U2U3U5_3=pars[995];
    long double a4_2U3U5_1=pars[996];
    long double a4_2U3U5_2=pars[997];
    long double a4_2U3U5_3=pars[998];
    long double b4_2U3U4U5_1=pars[999];
    long double b4_2U3U4U5_2=pars[1000];
    long double b4_2U3U4U5_3=pars[1001];
    long double a6_2U3U5_1=pars[1002];
    long double a6_2U3U5_2=pars[1003];
    long double a6_2U3U5_3=pars[1004];
    long double b6_2U3U5U6_1=pars[1005];
    long double k_1_2U3U5U6=pars[1006];
    long double kr_1_2U3U5U6=pars[1007];
    long double b6_2U3U5U6_2=pars[1008];
    long double k_2_2U3U5U6=pars[1009];
    long double b6_2U3U5U6_3=pars[1010];
    long double k_3_2U3U5U6=pars[1011];
    long double a1_2U3U6_1=pars[1012];
    long double a1_2U3U6_2=pars[1013];
    long double a1_2U3U6_3=pars[1014];
    long double b1_1U2U3U6_1=pars[1015];
    long double b1_1U2U3U6_2=pars[1016];
    long double b1_1U2U3U6_3=pars[1017];
    long double a4_2U3U6_1=pars[1018];
    long double a4_2U3U6_2=pars[1019];
    long double a4_2U3U6_3=pars[1020];
    long double b4_2U3U4U6_1=pars[1021];
    long double b4_2U3U4U6_2=pars[1022];
    long double b4_2U3U4U6_3=pars[1023];
    long double a5_2U3U6_1=pars[1024];
    long double a5_2U3U6_2=pars[1025];
    long double a5_2U3U6_3=pars[1026];
    long double b5_2U3U5U6_1=pars[1027];
    long double b5_2U3U5U6_2=pars[1028];
    long double b5_2U3U5U6_3=pars[1029];
    long double a1_2U4U5_1=pars[1030];
    long double a1_2U4U5_2=pars[1031];
    long double a1_2U4U5_3=pars[1032];
    long double b1_1U2U4U5_1=pars[1033];
    long double b1_1U2U4U5_2=pars[1034];
    long double b1_1U2U4U5_3=pars[1035];
    long double a3_2U4U5_1=pars[1036];
    long double a3_2U4U5_2=pars[1037];
    long double a3_2U4U5_3=pars[1038];
    long double b3_2U3U4U5_1=pars[1039];
    long double b3_2U3U4U5_2=pars[1040];
    long double b3_2U3U4U5_3=pars[1041];
    long double a6_2U4U5_1=pars[1042];
    long double a6_2U4U5_2=pars[1043];
    long double a6_2U4U5_3=pars[1044];
    long double b6_2U4U5U6_1=pars[1045];
    long double k_1_2U4U5U6=pars[1046];
    long double kr_1_2U4U5U6=pars[1047];
    long double b6_2U4U5U6_2=pars[1048];
    long double k_2_2U4U5U6=pars[1049];
    long double b6_2U4U5U6_3=pars[1050];
    long double k_3_2U4U5U6=pars[1051];
    long double a1_2U4U6_1=pars[1052];
    long double a1_2U4U6_2=pars[1053];
    long double a1_2U4U6_3=pars[1054];
    long double b1_1U2U4U6_1=pars[1055];
    long double b1_1U2U4U6_2=pars[1056];
    long double b1_1U2U4U6_3=pars[1057];
    long double a3_2U4U6_1=pars[1058];
    long double a3_2U4U6_2=pars[1059];
    long double a3_2U4U6_3=pars[1060];
    long double b3_2U3U4U6_1=pars[1061];
    long double b3_2U3U4U6_2=pars[1062];
    long double b3_2U3U4U6_3=pars[1063];
    long double a5_2U4U6_1=pars[1064];
    long double a5_2U4U6_2=pars[1065];
    long double a5_2U4U6_3=pars[1066];
    long double b5_2U4U5U6_1=pars[1067];
    long double b5_2U4U5U6_2=pars[1068];
    long double b5_2U4U5U6_3=pars[1069];
    long double a1_2U5U6_1=pars[1070];
    long double a1_2U5U6_2=pars[1071];
    long double a1_2U5U6_3=pars[1072];
    long double b1_1U2U5U6_1=pars[1073];
    long double b1_1U2U5U6_2=pars[1074];
    long double b1_1U2U5U6_3=pars[1075];
    long double a3_2U5U6_1=pars[1076];
    long double a3_2U5U6_2=pars[1077];
    long double a3_2U5U6_3=pars[1078];
    long double b3_2U3U5U6_1=pars[1079];
    long double b3_2U3U5U6_2=pars[1080];
    long double b3_2U3U5U6_3=pars[1081];
    long double a4_2U5U6_1=pars[1082];
    long double a4_2U5U6_2=pars[1083];
    long double a4_2U5U6_3=pars[1084];
    long double b4_2U4U5U6_1=pars[1085];
    long double b4_2U4U5U6_2=pars[1086];
    long double b4_2U4U5U6_3=pars[1087];
    long double a1_3U4U5_1=pars[1088];
    long double a1_3U4U5_2=pars[1089];
    long double a1_3U4U5_3=pars[1090];
    long double b1_1U3U4U5_1=pars[1091];
    long double b1_1U3U4U5_2=pars[1092];
    long double b1_1U3U4U5_3=pars[1093];
    long double a2_3U4U5_1=pars[1094];
    long double a2_3U4U5_2=pars[1095];
    long double a2_3U4U5_3=pars[1096];
    long double b2_2U3U4U5_1=pars[1097];
    long double b2_2U3U4U5_2=pars[1098];
    long double b2_2U3U4U5_3=pars[1099];
    long double a6_3U4U5_1=pars[1100];
    long double a6_3U4U5_2=pars[1101];
    long double a6_3U4U5_3=pars[1102];
    long double b6_3U4U5U6_1=pars[1103];
    long double k_1_3U4U5U6=pars[1104];
    long double kr_1_3U4U5U6=pars[1105];
    long double b6_3U4U5U6_2=pars[1106];
    long double k_2_3U4U5U6=pars[1107];
    long double b6_3U4U5U6_3=pars[1108];
    long double k_3_3U4U5U6=pars[1109];
    long double a1_3U4U6_1=pars[1110];
    long double a1_3U4U6_2=pars[1111];
    long double a1_3U4U6_3=pars[1112];
    long double b1_1U3U4U6_1=pars[1113];
    long double b1_1U3U4U6_2=pars[1114];
    long double b1_1U3U4U6_3=pars[1115];
    long double a2_3U4U6_1=pars[1116];
    long double a2_3U4U6_2=pars[1117];
    long double a2_3U4U6_3=pars[1118];
    long double b2_2U3U4U6_1=pars[1119];
    long double b2_2U3U4U6_2=pars[1120];
    long double b2_2U3U4U6_3=pars[1121];
    long double a5_3U4U6_1=pars[1122];
    long double a5_3U4U6_2=pars[1123];
    long double a5_3U4U6_3=pars[1124];
    long double b5_3U4U5U6_1=pars[1125];
    long double b5_3U4U5U6_2=pars[1126];
    long double b5_3U4U5U6_3=pars[1127];
    long double a1_3U5U6_1=pars[1128];
    long double a1_3U5U6_2=pars[1129];
    long double a1_3U5U6_3=pars[1130];
    long double b1_1U3U5U6_1=pars[1131];
    long double b1_1U3U5U6_2=pars[1132];
    long double b1_1U3U5U6_3=pars[1133];
    long double a2_3U5U6_1=pars[1134];
    long double a2_3U5U6_2=pars[1135];
    long double a2_3U5U6_3=pars[1136];
    long double b2_2U3U5U6_1=pars[1137];
    long double b2_2U3U5U6_2=pars[1138];
    long double b2_2U3U5U6_3=pars[1139];
    long double a4_3U5U6_1=pars[1140];
    long double a4_3U5U6_2=pars[1141];
    long double a4_3U5U6_3=pars[1142];
    long double b4_3U4U5U6_1=pars[1143];
    long double b4_3U4U5U6_2=pars[1144];
    long double b4_3U4U5U6_3=pars[1145];
    long double a1_4U5U6_1=pars[1146];
    long double a1_4U5U6_2=pars[1147];
    long double a1_4U5U6_3=pars[1148];
    long double b1_1U4U5U6_1=pars[1149];
    long double b1_1U4U5U6_2=pars[1150];
    long double b1_1U4U5U6_3=pars[1151];
    long double a2_4U5U6_1=pars[1152];
    long double a2_4U5U6_2=pars[1153];
    long double a2_4U5U6_3=pars[1154];
    long double b2_2U4U5U6_1=pars[1155];
    long double b2_2U4U5U6_2=pars[1156];
    long double b2_2U4U5U6_3=pars[1157];
    long double a3_4U5U6_1=pars[1158];
    long double a3_4U5U6_2=pars[1159];
    long double a3_4U5U6_3=pars[1160];
    long double b3_3U4U5U6_1=pars[1161];
    long double b3_3U4U5U6_2=pars[1162];
    long double b3_3U4U5U6_3=pars[1163];
    long double a5_1U2U3U4_1=pars[1164];
    long double a5_1U2U3U4_2=pars[1165];
    long double a5_1U2U3U4_3=pars[1166];
    long double b5_1U2U3U4U5_1=pars[1167];
    long double k_1_1U2U3U4U5=pars[1168];
    long double kr_1_1U2U3U4U5=pars[1169];
    long double b5_1U2U3U4U5_2=pars[1170];
    long double k_2_1U2U3U4U5=pars[1171];
    long double b5_1U2U3U4U5_3=pars[1172];
    long double k_3_1U2U3U4U5=pars[1173];
    long double a6_1U2U3U4_1=pars[1174];
    long double a6_1U2U3U4_2=pars[1175];
    long double a6_1U2U3U4_3=pars[1176];
    long double b6_1U2U3U4U6_1=pars[1177];
    long double k_1_1U2U3U4U6=pars[1178];
    long double kr_1_1U2U3U4U6=pars[1179];
    long double b6_1U2U3U4U6_2=pars[1180];
    long double k_2_1U2U3U4U6=pars[1181];
    long double b6_1U2U3U4U6_3=pars[1182];
    long double k_3_1U2U3U4U6=pars[1183];
    long double a4_1U2U3U5_1=pars[1184];
    long double a4_1U2U3U5_2=pars[1185];
    long double a4_1U2U3U5_3=pars[1186];
    long double b4_1U2U3U4U5_1=pars[1187];
    long double b4_1U2U3U4U5_2=pars[1188];
    long double b4_1U2U3U4U5_3=pars[1189];
    long double a6_1U2U3U5_1=pars[1190];
    long double a6_1U2U3U5_2=pars[1191];
    long double a6_1U2U3U5_3=pars[1192];
    long double b6_1U2U3U5U6_1=pars[1193];
    long double k_1_1U2U3U5U6=pars[1194];
    long double kr_1_1U2U3U5U6=pars[1195];
    long double b6_1U2U3U5U6_2=pars[1196];
    long double k_2_1U2U3U5U6=pars[1197];
    long double b6_1U2U3U5U6_3=pars[1198];
    long double k_3_1U2U3U5U6=pars[1199];
    long double a4_1U2U3U6_1=pars[1200];
    long double a4_1U2U3U6_2=pars[1201];
    long double a4_1U2U3U6_3=pars[1202];
    long double b4_1U2U3U4U6_1=pars[1203];
    long double b4_1U2U3U4U6_2=pars[1204];
    long double b4_1U2U3U4U6_3=pars[1205];
    long double a5_1U2U3U6_1=pars[1206];
    long double a5_1U2U3U6_2=pars[1207];
    long double a5_1U2U3U6_3=pars[1208];
    long double b5_1U2U3U5U6_1=pars[1209];
    long double b5_1U2U3U5U6_2=pars[1210];
    long double b5_1U2U3U5U6_3=pars[1211];
    long double a3_1U2U4U5_1=pars[1212];
    long double a3_1U2U4U5_2=pars[1213];
    long double a3_1U2U4U5_3=pars[1214];
    long double b3_1U2U3U4U5_1=pars[1215];
    long double b3_1U2U3U4U5_2=pars[1216];
    long double b3_1U2U3U4U5_3=pars[1217];
    long double a6_1U2U4U5_1=pars[1218];
    long double a6_1U2U4U5_2=pars[1219];
    long double a6_1U2U4U5_3=pars[1220];
    long double b6_1U2U4U5U6_1=pars[1221];
    long double k_1_1U2U4U5U6=pars[1222];
    long double kr_1_1U2U4U5U6=pars[1223];
    long double b6_1U2U4U5U6_2=pars[1224];
    long double k_2_1U2U4U5U6=pars[1225];
    long double b6_1U2U4U5U6_3=pars[1226];
    long double k_3_1U2U4U5U6=pars[1227];
    long double a3_1U2U4U6_1=pars[1228];
    long double a3_1U2U4U6_2=pars[1229];
    long double a3_1U2U4U6_3=pars[1230];
    long double b3_1U2U3U4U6_1=pars[1231];
    long double b3_1U2U3U4U6_2=pars[1232];
    long double b3_1U2U3U4U6_3=pars[1233];
    long double a5_1U2U4U6_1=pars[1234];
    long double a5_1U2U4U6_2=pars[1235];
    long double a5_1U2U4U6_3=pars[1236];
    long double b5_1U2U4U5U6_1=pars[1237];
    long double b5_1U2U4U5U6_2=pars[1238];
    long double b5_1U2U4U5U6_3=pars[1239];
    long double a3_1U2U5U6_1=pars[1240];
    long double a3_1U2U5U6_2=pars[1241];
    long double a3_1U2U5U6_3=pars[1242];
    long double b3_1U2U3U5U6_1=pars[1243];
    long double b3_1U2U3U5U6_2=pars[1244];
    long double b3_1U2U3U5U6_3=pars[1245];
    long double a4_1U2U5U6_1=pars[1246];
    long double a4_1U2U5U6_2=pars[1247];
    long double a4_1U2U5U6_3=pars[1248];
    long double b4_1U2U4U5U6_1=pars[1249];
    long double b4_1U2U4U5U6_2=pars[1250];
    long double b4_1U2U4U5U6_3=pars[1251];
    long double a2_1U3U4U5_1=pars[1252];
    long double a2_1U3U4U5_2=pars[1253];
    long double a2_1U3U4U5_3=pars[1254];
    long double b2_1U2U3U4U5_1=pars[1255];
    long double b2_1U2U3U4U5_2=pars[1256];
    long double b2_1U2U3U4U5_3=pars[1257];
    long double a6_1U3U4U5_1=pars[1258];
    long double a6_1U3U4U5_2=pars[1259];
    long double a6_1U3U4U5_3=pars[1260];
    long double b6_1U3U4U5U6_1=pars[1261];
    long double k_1_1U3U4U5U6=pars[1262];
    long double kr_1_1U3U4U5U6=pars[1263];
    long double b6_1U3U4U5U6_2=pars[1264];
    long double k_2_1U3U4U5U6=pars[1265];
    long double b6_1U3U4U5U6_3=pars[1266];
    long double k_3_1U3U4U5U6=pars[1267];
    long double a2_1U3U4U6_1=pars[1268];
    long double a2_1U3U4U6_2=pars[1269];
    long double a2_1U3U4U6_3=pars[1270];
    long double b2_1U2U3U4U6_1=pars[1271];
    long double b2_1U2U3U4U6_2=pars[1272];
    long double b2_1U2U3U4U6_3=pars[1273];
    long double a5_1U3U4U6_1=pars[1274];
    long double a5_1U3U4U6_2=pars[1275];
    long double a5_1U3U4U6_3=pars[1276];
    long double b5_1U3U4U5U6_1=pars[1277];
    long double b5_1U3U4U5U6_2=pars[1278];
    long double b5_1U3U4U5U6_3=pars[1279];
    long double a2_1U3U5U6_1=pars[1280];
    long double a2_1U3U5U6_2=pars[1281];
    long double a2_1U3U5U6_3=pars[1282];
    long double b2_1U2U3U5U6_1=pars[1283];
    long double b2_1U2U3U5U6_2=pars[1284];
    long double b2_1U2U3U5U6_3=pars[1285];
    long double a4_1U3U5U6_1=pars[1286];
    long double a4_1U3U5U6_2=pars[1287];
    long double a4_1U3U5U6_3=pars[1288];
    long double b4_1U3U4U5U6_1=pars[1289];
    long double b4_1U3U4U5U6_2=pars[1290];
    long double b4_1U3U4U5U6_3=pars[1291];
    long double a2_1U4U5U6_1=pars[1292];
    long double a2_1U4U5U6_2=pars[1293];
    long double a2_1U4U5U6_3=pars[1294];
    long double b2_1U2U4U5U6_1=pars[1295];
    long double b2_1U2U4U5U6_2=pars[1296];
    long double b2_1U2U4U5U6_3=pars[1297];
    long double a3_1U4U5U6_1=pars[1298];
    long double a3_1U4U5U6_2=pars[1299];
    long double a3_1U4U5U6_3=pars[1300];
    long double b3_1U3U4U5U6_1=pars[1301];
    long double b3_1U3U4U5U6_2=pars[1302];
    long double b3_1U3U4U5U6_3=pars[1303];
    long double a1_2U3U4U5_1=pars[1304];
    long double a1_2U3U4U5_2=pars[1305];
    long double a1_2U3U4U5_3=pars[1306];
    long double b1_1U2U3U4U5_1=pars[1307];
    long double b1_1U2U3U4U5_2=pars[1308];
    long double b1_1U2U3U4U5_3=pars[1309];
    long double a6_2U3U4U5_1=pars[1310];
    long double a6_2U3U4U5_2=pars[1311];
    long double a6_2U3U4U5_3=pars[1312];
    long double b6_2U3U4U5U6_1=pars[1313];
    long double k_1_2U3U4U5U6=pars[1314];
    long double kr_1_2U3U4U5U6=pars[1315];
    long double b6_2U3U4U5U6_2=pars[1316];
    long double k_2_2U3U4U5U6=pars[1317];
    long double b6_2U3U4U5U6_3=pars[1318];
    long double k_3_2U3U4U5U6=pars[1319];
    long double a1_2U3U4U6_1=pars[1320];
    long double a1_2U3U4U6_2=pars[1321];
    long double a1_2U3U4U6_3=pars[1322];
    long double b1_1U2U3U4U6_1=pars[1323];
    long double b1_1U2U3U4U6_2=pars[1324];
    long double b1_1U2U3U4U6_3=pars[1325];
    long double a5_2U3U4U6_1=pars[1326];
    long double a5_2U3U4U6_2=pars[1327];
    long double a5_2U3U4U6_3=pars[1328];
    long double b5_2U3U4U5U6_1=pars[1329];
    long double b5_2U3U4U5U6_2=pars[1330];
    long double b5_2U3U4U5U6_3=pars[1331];
    long double a1_2U3U5U6_1=pars[1332];
    long double a1_2U3U5U6_2=pars[1333];
    long double a1_2U3U5U6_3=pars[1334];
    long double b1_1U2U3U5U6_1=pars[1335];
    long double b1_1U2U3U5U6_2=pars[1336];
    long double b1_1U2U3U5U6_3=pars[1337];
    long double a4_2U3U5U6_1=pars[1338];
    long double a4_2U3U5U6_2=pars[1339];
    long double a4_2U3U5U6_3=pars[1340];
    long double b4_2U3U4U5U6_1=pars[1341];
    long double b4_2U3U4U5U6_2=pars[1342];
    long double b4_2U3U4U5U6_3=pars[1343];
    long double a1_2U4U5U6_1=pars[1344];
    long double a1_2U4U5U6_2=pars[1345];
    long double a1_2U4U5U6_3=pars[1346];
    long double b1_1U2U4U5U6_1=pars[1347];
    long double b1_1U2U4U5U6_2=pars[1348];
    long double b1_1U2U4U5U6_3=pars[1349];
    long double a3_2U4U5U6_1=pars[1350];
    long double a3_2U4U5U6_2=pars[1351];
    long double a3_2U4U5U6_3=pars[1352];
    long double b3_2U3U4U5U6_1=pars[1353];
    long double b3_2U3U4U5U6_2=pars[1354];
    long double b3_2U3U4U5U6_3=pars[1355];
    long double a1_3U4U5U6_1=pars[1356];
    long double a1_3U4U5U6_2=pars[1357];
    long double a1_3U4U5U6_3=pars[1358];
    long double b1_1U3U4U5U6_1=pars[1359];
    long double b1_1U3U4U5U6_2=pars[1360];
    long double b1_1U3U4U5U6_3=pars[1361];
    long double a2_3U4U5U6_1=pars[1362];
    long double a2_3U4U5U6_2=pars[1363];
    long double a2_3U4U5U6_3=pars[1364];
    long double b2_2U3U4U5U6_1=pars[1365];
    long double b2_2U3U4U5U6_2=pars[1366];
    long double b2_2U3U4U5U6_3=pars[1367];
    long double a6_1U2U3U4U5_1=pars[1368];
    long double a6_1U2U3U4U5_2=pars[1369];
    long double a6_1U2U3U4U5_3=pars[1370];
    long double b6_1U2U3U4U5U6_1=pars[1371];
    long double k_1_1U2U3U4U5U6=pars[1372];
    long double kr_1_1U2U3U4U5U6=pars[1373];
    long double b6_1U2U3U4U5U6_2=pars[1374];
    long double k_2_1U2U3U4U5U6=pars[1375];
    long double b6_1U2U3U4U5U6_3=pars[1376];
    long double k_3_1U2U3U4U5U6=pars[1377];
    long double a5_1U2U3U4U6_1=pars[1378];
    long double a5_1U2U3U4U6_2=pars[1379];
    long double a5_1U2U3U4U6_3=pars[1380];
    long double b5_1U2U3U4U5U6_1=pars[1381];
    long double b5_1U2U3U4U5U6_2=pars[1382];
    long double b5_1U2U3U4U5U6_3=pars[1383];
    long double a4_1U2U3U5U6_1=pars[1384];
    long double a4_1U2U3U5U6_2=pars[1385];
    long double a4_1U2U3U5U6_3=pars[1386];
    long double b4_1U2U3U4U5U6_1=pars[1387];
    long double b4_1U2U3U4U5U6_2=pars[1388];
    long double b4_1U2U3U4U5U6_3=pars[1389];
    long double a3_1U2U4U5U6_1=pars[1390];
    long double a3_1U2U4U5U6_2=pars[1391];
    long double a3_1U2U4U5U6_3=pars[1392];
    long double b3_1U2U3U4U5U6_1=pars[1393];
    long double b3_1U2U3U4U5U6_2=pars[1394];
    long double b3_1U2U3U4U5U6_3=pars[1395];
    long double a2_1U3U4U5U6_1=pars[1396];
    long double a2_1U3U4U5U6_2=pars[1397];
    long double a2_1U3U4U5U6_3=pars[1398];
    long double b2_1U2U3U4U5U6_1=pars[1399];
    long double b2_1U2U3U4U5U6_2=pars[1400];
    long double b2_1U2U3U4U5U6_3=pars[1401];
    long double a1_2U3U4U5U6_1=pars[1402];
    long double a1_2U3U4U5U6_2=pars[1403];
    long double a1_2U3U4U5U6_3=pars[1404];
    long double b1_1U2U3U4U5U6_1=pars[1405];
    long double b1_1U2U3U4U5U6_2=pars[1406];
    long double b1_1U2U3U4U5U6_3=pars[1407];
    Matrix<InternalType, Dynamic, Dynamic> L = Matrix<InternalType, Dynamic, Dynamic>::Zero(192, 192);
    L<<0,b1_1_1,b2_2_1,b3_3_1,b4_4_1,b5_5_1,b6_6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a1_0_1*x,0,0,0,0,0,0,b2_1U2_1,b3_1U3_1,b4_1U4_1,b5_1U5_1,b6_1U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a2_0_1*x,0,0,0,0,0,0,b1_1U2_1,0,0,0,0,b3_2U3_1,b4_2U4_1,b5_2U5_1,b6_2U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a3_0_1*x,0,0,0,0,0,0,0,b1_1U3_1,0,0,0,b2_2U3_1,0,0,0,b4_3U4_1,b5_3U5_1,b6_3U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a4_0_1*x,0,0,0,0,0,0,0,0,b1_1U4_1,0,0,0,b2_2U4_1,0,0,b3_3U4_1,0,0,b5_4U5_1,b6_4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a5_0_1*x,0,0,0,0,0,0,0,0,0,b1_1U5_1,0,0,0,b2_2U5_1,0,0,b3_3U5_1,0,b4_4U5_1,0,b6_5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
a6_0_1*x,0,0,0,0,0,0,0,0,0,0,b1_1U6_1,0,0,0,b2_2U6_1,0,0,b3_3U6_1,0,b4_4U6_1,b5_5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,a2_1_1*x,a1_2_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3_1,b4_1U2U4_1,b5_1U2U5_1,b6_1U2U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,a3_1_1*x,0,a1_3_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3_1,0,0,0,b4_1U3U4_1,b5_1U3U5_1,b6_1U3U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,a4_1_1*x,0,0,a1_4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4_1,0,0,b3_1U3U4_1,0,0,b5_1U4U5_1,b6_1U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,a5_1_1*x,0,0,0,a1_5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5_1,0,0,b3_1U3U5_1,0,b4_1U4U5_1,0,b6_1U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,a6_1_1*x,0,0,0,0,a1_6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U6_1,0,0,b3_1U3U6_1,0,b4_1U4U6_1,b5_1U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,a3_2_1*x,a2_3_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3_1,0,0,0,0,0,0,0,0,0,b4_2U3U4_1,b5_2U3U5_1,b6_2U3U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,a4_2_1*x,0,a2_4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4_1,0,0,0,0,0,0,0,0,b3_2U3U4_1,0,0,b5_2U4U5_1,b6_2U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,a5_2_1*x,0,0,a2_5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5_1,0,0,0,0,0,0,0,0,b3_2U3U5_1,0,b4_2U4U5_1,0,b6_2U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,a6_2_1*x,0,0,0,a2_6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U6_1,0,0,0,0,0,0,0,0,b3_2U3U6_1,0,b4_2U4U6_1,b5_2U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,a4_3_1*x,a3_4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4_1,0,0,0,0,0,b2_2U3U4_1,0,0,0,0,0,b5_3U4U5_1,b6_3U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,a5_3_1*x,0,a3_5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5_1,0,0,0,0,0,b2_2U3U5_1,0,0,0,0,b4_3U4U5_1,0,b6_3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,a6_3_1*x,0,0,a3_6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U6_1,0,0,0,0,0,b2_2U3U6_1,0,0,0,0,b4_3U4U6_1,b5_3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,a5_4_1*x,a4_5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5_1,0,0,0,0,0,b2_2U4U5_1,0,0,b3_3U4U5_1,0,0,b6_4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,a6_4_1*x,0,a4_6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U6_1,0,0,0,0,0,b2_2U4U6_1,0,0,b3_3U4U6_1,0,b5_4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,a6_5_1*x,a5_6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U5U6_1,0,0,0,0,0,b2_2U5U6_1,0,0,b3_3U5U6_1,b4_4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,a3_1U2_1*x,a2_1U3_1*x,0,0,0,a1_2U3_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4_1,b5_1U2U3U5_1,b6_1U2U3U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,a4_1U2_1*x,0,a2_1U4_1*x,0,0,0,a1_2U4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4_1,0,0,b5_1U2U4U5_1,b6_1U2U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,a5_1U2_1*x,0,0,a2_1U5_1*x,0,0,0,a1_2U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5_1,0,b4_1U2U4U5_1,0,b6_1U2U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,a6_1U2_1*x,0,0,0,a2_1U6_1*x,0,0,0,a1_2U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U6_1,0,b4_1U2U4U6_1,b5_1U2U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,a4_1U3_1*x,a3_1U4_1*x,0,0,0,0,0,0,a1_3U4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4_1,0,0,0,0,0,b5_1U3U4U5_1,b6_1U3U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,a5_1U3_1*x,0,a3_1U5_1*x,0,0,0,0,0,0,a1_3U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5_1,0,0,0,0,b4_1U3U4U5_1,0,b6_1U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,a6_1U3_1*x,0,0,a3_1U6_1*x,0,0,0,0,0,0,a1_3U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U6_1,0,0,0,0,b4_1U3U4U6_1,b5_1U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,a5_1U4_1*x,a4_1U5_1*x,0,0,0,0,0,0,0,0,a1_4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5_1,0,0,b3_1U3U4U5_1,0,0,b6_1U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,a6_1U4_1*x,0,a4_1U6_1*x,0,0,0,0,0,0,0,0,a1_4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U6_1,0,0,b3_1U3U4U6_1,0,b5_1U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,a6_1U5_1*x,a5_1U6_1*x,0,0,0,0,0,0,0,0,0,a1_5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5U6_1,0,0,b3_1U3U5U6_1,b4_1U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,a4_2U3_1*x,a3_2U4_1*x,0,0,a2_3U4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4_1,0,0,0,0,0,0,0,0,0,b5_2U3U4U5_1,b6_2U3U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3_1*x,0,a3_2U5_1*x,0,0,a2_3U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5_1,0,0,0,0,0,0,0,0,b4_2U3U4U5_1,0,b6_2U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3_1*x,0,0,a3_2U6_1*x,0,0,a2_3U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U6_1,0,0,0,0,0,0,0,0,b4_2U3U4U6_1,b5_2U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U4_1*x,a4_2U5_1*x,0,0,0,0,a2_4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5_1,0,0,0,0,0,0,b3_2U3U4U5_1,0,0,b6_2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4_1*x,0,a4_2U6_1*x,0,0,0,0,a2_4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U6_1,0,0,0,0,0,0,b3_2U3U4U6_1,0,b5_2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U5_1*x,a5_2U6_1*x,0,0,0,0,0,a2_5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5U6_1,0,0,0,0,0,0,b3_2U3U5U6_1,b4_2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_3U4_1*x,a4_3U5_1*x,0,a3_4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5_1,0,0,0,b2_2U3U4U5_1,0,0,0,b6_3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4_1*x,0,a4_3U6_1*x,0,a3_4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U6_1,0,0,0,b2_2U3U4U6_1,0,0,b5_3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U5_1*x,a5_3U6_1*x,0,0,a3_5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5U6_1,0,0,0,b2_2U3U5U6_1,0,b4_3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_4U5_1*x,a5_4U6_1*x,a4_5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5U6_1,0,0,0,b2_2U4U5U6_1,b3_3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_1*x,a3_1U2U4_1*x,0,0,a2_1U3U4_1*x,0,0,0,0,0,a1_2U3U4_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5_1,b6_1U2U3U4U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3_1*x,0,a3_1U2U5_1*x,0,0,a2_1U3U5_1*x,0,0,0,0,0,a1_2U3U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U5_1,0,b6_1U2U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3_1*x,0,0,a3_1U2U6_1*x,0,0,a2_1U3U6_1*x,0,0,0,0,0,a1_2U3U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U6_1,b5_1U2U3U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U4_1*x,a4_1U2U5_1*x,0,0,0,0,a2_1U4U5_1*x,0,0,0,0,0,a1_2U4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U5_1,0,0,b6_1U2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4_1*x,0,a4_1U2U6_1*x,0,0,0,0,a2_1U4U6_1*x,0,0,0,0,0,a1_2U4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U6_1,0,b5_1U2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U5_1*x,a5_1U2U6_1*x,0,0,0,0,0,a2_1U5U6_1*x,0,0,0,0,0,a1_2U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5U6_1,b4_1U2U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U3U4_1*x,a4_1U3U5_1*x,0,a3_1U4U5_1*x,0,0,0,0,0,0,0,0,a1_3U4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U5_1,0,0,0,b6_1U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4_1*x,0,a4_1U3U6_1*x,0,a3_1U4U6_1*x,0,0,0,0,0,0,0,0,a1_3U4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U6_1,0,0,b5_1U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U5_1*x,a5_1U3U6_1*x,0,0,a3_1U5U6_1*x,0,0,0,0,0,0,0,0,a1_3U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5U6_1,0,b4_1U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U4U5_1*x,a5_1U4U6_1*x,a4_1U5U6_1*x,0,0,0,0,0,0,0,0,0,a1_4U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5U6_1,b3_1U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3U4_1*x,a4_2U3U5_1*x,0,a3_2U4U5_1*x,0,0,a2_3U4U5_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U5_1,0,0,0,0,b6_2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4_1*x,0,a4_2U3U6_1*x,0,a3_2U4U6_1*x,0,0,a2_3U4U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U6_1,0,0,0,b5_2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U4U6,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U5_1*x,a5_2U3U6_1*x,0,0,a3_2U5U6_1*x,0,0,a2_3U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5U6_1,0,0,b4_2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U5U6,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4U5_1*x,a5_2U4U6_1*x,a4_2U5U6_1*x,0,0,0,a2_4U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5U6_1,0,b3_2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U4U5U6,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4U5_1*x,a5_3U4U6_1*x,a4_3U5U6_1*x,a3_4U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5U6_1,b2_2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_3U4U5U6,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3U4_1*x,a4_1U2U3U5_1*x,0,a3_1U2U4U5_1*x,0,0,a2_1U3U4U5_1*x,0,0,0,a1_2U3U4U5_1*x,0,0,0,0,0,0,0,0,0,0,b6_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U4U5,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4_1*x,0,a4_1U2U3U6_1*x,0,a3_1U2U4U6_1*x,0,0,a2_1U3U4U6_1*x,0,0,0,a1_2U3U4U6_1*x,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U4U6,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U5_1*x,a5_1U2U3U6_1*x,0,0,a3_1U2U5U6_1*x,0,0,a2_1U3U5U6_1*x,0,0,0,a1_2U3U5U6_1*x,0,0,0,0,0,0,0,0,b4_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U5U6,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4U5_1*x,a5_1U2U4U6_1*x,a4_1U2U5U6_1*x,0,0,0,a2_1U4U5U6_1*x,0,0,0,a1_2U4U5U6_1*x,0,0,0,0,0,0,0,b3_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U4U5U6,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4U5_1*x,a5_1U3U4U6_1*x,a4_1U3U5U6_1*x,a3_1U4U5U6_1*x,0,0,0,0,a1_3U4U5U6_1*x,0,0,0,0,0,0,b2_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U3U4U5U6,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4U5_1*x,a5_2U3U4U6_1*x,a4_2U3U5U6_1*x,a3_2U4U5U6_1*x,a2_3U4U5U6_1*x,0,0,0,0,0,0,b1_1U2U3U4U5U6_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_2U3U4U5U6,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4U5_1*x,a5_1U2U3U4U6_1*x,a4_1U2U3U5U6_1*x,a3_1U2U4U5U6_1*x,a2_1U3U4U5U6_1*x,a1_2U3U4U5U6_1*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,kr_1_1U2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_3_1U2U3U4U5U6,
k_1_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1_2,b2_2_2,b3_3_2,b4_4_2,b5_5_2,b6_6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,k_1_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a1_0_2*x,0,0,0,0,0,0,b2_1U2_2,b3_1U3_2,b4_1U4_2,b5_1U5_2,b6_1U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,k_1_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_0_2*x,0,0,0,0,0,0,b1_1U2_2,0,0,0,0,b3_2U3_2,b4_2U4_2,b5_2U5_2,b6_2U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,k_1_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_0_2*x,0,0,0,0,0,0,0,b1_1U3_2,0,0,0,b2_2U3_2,0,0,0,b4_3U4_2,b5_3U5_2,b6_3U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,k_1_4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_0_2*x,0,0,0,0,0,0,0,0,b1_1U4_2,0,0,0,b2_2U4_2,0,0,b3_3U4_2,0,0,b5_4U5_2,b6_4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,k_1_5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_0_2*x,0,0,0,0,0,0,0,0,0,b1_1U5_2,0,0,0,b2_2U5_2,0,0,b3_3U5_2,0,b4_4U5_2,0,b6_5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,k_1_6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_0_2*x,0,0,0,0,0,0,0,0,0,0,b1_1U6_2,0,0,0,b2_2U6_2,0,0,b3_3U6_2,0,b4_4U6_2,b5_5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,k_1_1U2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_1_2*x,a1_2_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3_2,b4_1U2U4_2,b5_1U2U5_2,b6_1U2U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,k_1_1U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_1_2*x,0,a1_3_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3_2,0,0,0,b4_1U3U4_2,b5_1U3U5_2,b6_1U3U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,k_1_1U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1_2*x,0,0,a1_4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4_2,0,0,b3_1U3U4_2,0,0,b5_1U4U5_2,b6_1U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,k_1_1U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1_2*x,0,0,0,a1_5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5_2,0,0,b3_1U3U5_2,0,b4_1U4U5_2,0,b6_1U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,k_1_1U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1_2*x,0,0,0,0,a1_6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U6_2,0,0,b3_1U3U6_2,0,b4_1U4U6_2,b5_1U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_2_2*x,a2_3_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3_2,0,0,0,0,0,0,0,0,0,b4_2U3U4_2,b5_2U3U5_2,b6_2U3U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_2_2*x,0,a2_4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4_2,0,0,0,0,0,0,0,0,b3_2U3U4_2,0,0,b5_2U4U5_2,b6_2U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2_2*x,0,0,a2_5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5_2,0,0,0,0,0,0,0,0,b3_2U3U5_2,0,b4_2U4U5_2,0,b6_2U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2_2*x,0,0,0,a2_6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U6_2,0,0,0,0,0,0,0,0,b3_2U3U6_2,0,b4_2U4U6_2,b5_2U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_3_2*x,a3_4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4_2,0,0,0,0,0,b2_2U3U4_2,0,0,0,0,0,b5_3U4U5_2,b6_3U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_3_2*x,0,a3_5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5_2,0,0,0,0,0,b2_2U3U5_2,0,0,0,0,b4_3U4U5_2,0,b6_3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3_2*x,0,0,a3_6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U6_2,0,0,0,0,0,b2_2U3U6_2,0,0,0,0,b4_3U4U6_2,b5_3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_4_2*x,a4_5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5_2,0,0,0,0,0,b2_2U4U5_2,0,0,b3_3U4U5_2,0,0,b6_4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_4_2*x,0,a4_6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U6_2,0,0,0,0,0,b2_2U4U6_2,0,0,b3_3U4U6_2,0,b5_4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_5_2*x,a5_6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U5U6_2,0,0,0,0,0,b2_2U5U6_2,0,0,b3_3U5U6_2,b4_4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_1U2_2*x,a2_1U3_2*x,0,0,0,a1_2U3_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4_2,b5_1U2U3U5_2,b6_1U2U3U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U2_2*x,0,a2_1U4_2*x,0,0,0,a1_2U4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4_2,0,0,b5_1U2U4U5_2,b6_1U2U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2_2*x,0,0,a2_1U5_2*x,0,0,0,a1_2U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5_2,0,b4_1U2U4U5_2,0,b6_1U2U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2_2*x,0,0,0,a2_1U6_2*x,0,0,0,a1_2U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U6_2,0,b4_1U2U4U6_2,b5_1U2U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U3_2*x,a3_1U4_2*x,0,0,0,0,0,0,a1_3U4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4_2,0,0,0,0,0,b5_1U3U4U5_2,b6_1U3U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U3_2*x,0,a3_1U5_2*x,0,0,0,0,0,0,a1_3U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5_2,0,0,0,0,b4_1U3U4U5_2,0,b6_1U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3_2*x,0,0,a3_1U6_2*x,0,0,0,0,0,0,a1_3U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U6_2,0,0,0,0,b4_1U3U4U6_2,b5_1U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U4_2*x,a4_1U5_2*x,0,0,0,0,0,0,0,0,a1_4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5_2,0,0,b3_1U3U4U5_2,0,0,b6_1U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U4_2*x,0,a4_1U6_2*x,0,0,0,0,0,0,0,0,a1_4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U6_2,0,0,b3_1U3U4U6_2,0,b5_1U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U5_2*x,a5_1U6_2*x,0,0,0,0,0,0,0,0,0,a1_5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5U6_2,0,0,b3_1U3U5U6_2,b4_1U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_2U3_2*x,a3_2U4_2*x,0,0,a2_3U4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4_2,0,0,0,0,0,0,0,0,0,b5_2U3U4U5_2,b6_2U3U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3_2*x,0,a3_2U5_2*x,0,0,a2_3U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5_2,0,0,0,0,0,0,0,0,b4_2U3U4U5_2,0,b6_2U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3_2*x,0,0,a3_2U6_2*x,0,0,a2_3U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U6_2,0,0,0,0,0,0,0,0,b4_2U3U4U6_2,b5_2U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U4_2*x,a4_2U5_2*x,0,0,0,0,a2_4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5_2,0,0,0,0,0,0,b3_2U3U4U5_2,0,0,b6_2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4_2*x,0,a4_2U6_2*x,0,0,0,0,a2_4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U6_2,0,0,0,0,0,0,b3_2U3U4U6_2,0,b5_2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U5_2*x,a5_2U6_2*x,0,0,0,0,0,a2_5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5U6_2,0,0,0,0,0,0,b3_2U3U5U6_2,b4_2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_3U4_2*x,a4_3U5_2*x,0,a3_4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5_2,0,0,0,b2_2U3U4U5_2,0,0,0,b6_3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4_2*x,0,a4_3U6_2*x,0,a3_4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U6_2,0,0,0,b2_2U3U4U6_2,0,0,b5_3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U5_2*x,a5_3U6_2*x,0,0,a3_5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5U6_2,0,0,0,b2_2U3U5U6_2,0,b4_3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_4U5_2*x,a5_4U6_2*x,a4_5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5U6_2,0,0,0,b2_2U4U5U6_2,b3_3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_2*x,a3_1U2U4_2*x,0,0,a2_1U3U4_2*x,0,0,0,0,0,a1_2U3U4_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5_2,b6_1U2U3U4U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3_2*x,0,a3_1U2U5_2*x,0,0,a2_1U3U5_2*x,0,0,0,0,0,a1_2U3U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U5_2,0,b6_1U2U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3_2*x,0,0,a3_1U2U6_2*x,0,0,a2_1U3U6_2*x,0,0,0,0,0,a1_2U3U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U6_2,b5_1U2U3U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U4_2*x,a4_1U2U5_2*x,0,0,0,0,a2_1U4U5_2*x,0,0,0,0,0,a1_2U4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U5_2,0,0,b6_1U2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4_2*x,0,a4_1U2U6_2*x,0,0,0,0,a2_1U4U6_2*x,0,0,0,0,0,a1_2U4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U6_2,0,b5_1U2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U5_2*x,a5_1U2U6_2*x,0,0,0,0,0,a2_1U5U6_2*x,0,0,0,0,0,a1_2U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5U6_2,b4_1U2U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U3U4_2*x,a4_1U3U5_2*x,0,a3_1U4U5_2*x,0,0,0,0,0,0,0,0,a1_3U4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U5_2,0,0,0,b6_1U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4_2*x,0,a4_1U3U6_2*x,0,a3_1U4U6_2*x,0,0,0,0,0,0,0,0,a1_3U4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U6_2,0,0,b5_1U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U5_2*x,a5_1U3U6_2*x,0,0,a3_1U5U6_2*x,0,0,0,0,0,0,0,0,a1_3U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5U6_2,0,b4_1U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U4U5_2*x,a5_1U4U6_2*x,a4_1U5U6_2*x,0,0,0,0,0,0,0,0,0,a1_4U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5U6_2,b3_1U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3U4_2*x,a4_2U3U5_2*x,0,a3_2U4U5_2*x,0,0,a2_3U4U5_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U5_2,0,0,0,0,b6_2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4_2*x,0,a4_2U3U6_2*x,0,a3_2U4U6_2*x,0,0,a2_3U4U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U6_2,0,0,0,b5_2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U5_2*x,a5_2U3U6_2*x,0,0,a3_2U5U6_2*x,0,0,a2_3U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5U6_2,0,0,b4_2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4U5_2*x,a5_2U4U6_2*x,a4_2U5U6_2*x,0,0,0,a2_4U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5U6_2,0,b3_2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4U5_2*x,a5_3U4U6_2*x,a4_3U5U6_2*x,a3_4U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5U6_2,b2_2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3U4_2*x,a4_1U2U3U5_2*x,0,a3_1U2U4U5_2*x,0,0,a2_1U3U4U5_2*x,0,0,0,a1_2U3U4U5_2*x,0,0,0,0,0,0,0,0,0,0,b6_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4_2*x,0,a4_1U2U3U6_2*x,0,a3_1U2U4U6_2*x,0,0,a2_1U3U4U6_2*x,0,0,0,a1_2U3U4U6_2*x,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U5_2*x,a5_1U2U3U6_2*x,0,0,a3_1U2U5U6_2*x,0,0,a2_1U3U5U6_2*x,0,0,0,a1_2U3U5U6_2*x,0,0,0,0,0,0,0,0,b4_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4U5_2*x,a5_1U2U4U6_2*x,a4_1U2U5U6_2*x,0,0,0,a2_1U4U5U6_2*x,0,0,0,a1_2U4U5U6_2*x,0,0,0,0,0,0,0,b3_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4U5_2*x,a5_1U3U4U6_2*x,a4_1U3U5U6_2*x,a3_1U4U5U6_2*x,0,0,0,0,a1_3U4U5U6_2*x,0,0,0,0,0,0,b2_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4U5_2*x,a5_2U3U4U6_2*x,a4_2U3U5U6_2*x,a3_2U4U5U6_2*x,a2_3U4U5U6_2*x,0,0,0,0,0,0,b1_1U2U3U4U5U6_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_1_1U2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4U5_2*x,a5_1U2U3U4U6_2*x,a4_1U2U3U5U6_2*x,a3_1U2U4U5U6_2*x,a2_1U3U4U5U6_2*x,a1_2U3U4U5U6_2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1_3,b2_2_3,b3_3_3,b4_4_3,b5_5_3,b6_6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a1_0_3*x,0,0,0,0,0,0,b2_1U2_3,b3_1U3_3,b4_1U4_3,b5_1U5_3,b6_1U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_0_3*x,0,0,0,0,0,0,b1_1U2_3,0,0,0,0,b3_2U3_3,b4_2U4_3,b5_2U5_3,b6_2U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_0_3*x,0,0,0,0,0,0,0,b1_1U3_3,0,0,0,b2_2U3_3,0,0,0,b4_3U4_3,b5_3U5_3,b6_3U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_0_3*x,0,0,0,0,0,0,0,0,b1_1U4_3,0,0,0,b2_2U4_3,0,0,b3_3U4_3,0,0,b5_4U5_3,b6_4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_0_3*x,0,0,0,0,0,0,0,0,0,b1_1U5_3,0,0,0,b2_2U5_3,0,0,b3_3U5_3,0,b4_4U5_3,0,b6_5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_0_3*x,0,0,0,0,0,0,0,0,0,0,b1_1U6_3,0,0,0,b2_2U6_3,0,0,b3_3U6_3,0,b4_4U6_3,b5_5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a2_1_3*x,a1_2_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3_3,b4_1U2U4_3,b5_1U2U5_3,b6_1U2U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_1_3*x,0,a1_3_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3_3,0,0,0,b4_1U3U4_3,b5_1U3U5_3,b6_1U3U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1_3*x,0,0,a1_4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4_3,0,0,b3_1U3U4_3,0,0,b5_1U4U5_3,b6_1U4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1_3*x,0,0,0,a1_5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5_3,0,0,b3_1U3U5_3,0,b4_1U4U5_3,0,b6_1U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1_3*x,0,0,0,0,a1_6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U6_3,0,0,b3_1U3U6_3,0,b4_1U4U6_3,b5_1U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_2_3*x,a2_3_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3_3,0,0,0,0,0,0,0,0,0,b4_2U3U4_3,b5_2U3U5_3,b6_2U3U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_2_3*x,0,a2_4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4_3,0,0,0,0,0,0,0,0,b3_2U3U4_3,0,0,b5_2U4U5_3,b6_2U4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2_3*x,0,0,a2_5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5_3,0,0,0,0,0,0,0,0,b3_2U3U5_3,0,b4_2U4U5_3,0,b6_2U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2_3*x,0,0,0,a2_6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U6_3,0,0,0,0,0,0,0,0,b3_2U3U6_3,0,b4_2U4U6_3,b5_2U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_3_3*x,a3_4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4_3,0,0,0,0,0,b2_2U3U4_3,0,0,0,0,0,b5_3U4U5_3,b6_3U4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_3_3*x,0,a3_5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5_3,0,0,0,0,0,b2_2U3U5_3,0,0,0,0,b4_3U4U5_3,0,b6_3U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3_3*x,0,0,a3_6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U6_3,0,0,0,0,0,b2_2U3U6_3,0,0,0,0,b4_3U4U6_3,b5_3U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_4_3*x,a4_5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5_3,0,0,0,0,0,b2_2U4U5_3,0,0,b3_3U4U5_3,0,0,b6_4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_4_3*x,0,a4_6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U6_3,0,0,0,0,0,b2_2U4U6_3,0,0,b3_3U4U6_3,0,b5_4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_5_3*x,a5_6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U5U6_3,0,0,0,0,0,b2_2U5U6_3,0,0,b3_3U5U6_3,b4_4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a3_1U2_3*x,a2_1U3_3*x,0,0,0,a1_2U3_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4_3,b5_1U2U3U5_3,b6_1U2U3U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U2_3*x,0,a2_1U4_3*x,0,0,0,a1_2U4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4_3,0,0,b5_1U2U4U5_3,b6_1U2U4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2_3*x,0,0,a2_1U5_3*x,0,0,0,a1_2U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5_3,0,b4_1U2U4U5_3,0,b6_1U2U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2_3*x,0,0,0,a2_1U6_3*x,0,0,0,a1_2U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U6_3,0,b4_1U2U4U6_3,b5_1U2U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U3_3*x,a3_1U4_3*x,0,0,0,0,0,0,a1_3U4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4_3,0,0,0,0,0,b5_1U3U4U5_3,b6_1U3U4U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U3_3*x,0,a3_1U5_3*x,0,0,0,0,0,0,a1_3U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5_3,0,0,0,0,b4_1U3U4U5_3,0,b6_1U3U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3_3*x,0,0,a3_1U6_3*x,0,0,0,0,0,0,a1_3U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U6_3,0,0,0,0,b4_1U3U4U6_3,b5_1U3U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U4_3*x,a4_1U5_3*x,0,0,0,0,0,0,0,0,a1_4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5_3,0,0,b3_1U3U4U5_3,0,0,b6_1U4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U4_3*x,0,a4_1U6_3*x,0,0,0,0,0,0,0,0,a1_4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U6_3,0,0,b3_1U3U4U6_3,0,b5_1U4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U5_3*x,a5_1U6_3*x,0,0,0,0,0,0,0,0,0,a1_5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U5U6_3,0,0,b3_1U3U5U6_3,b4_1U4U5U6_3,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_2U3_3*x,a3_2U4_3*x,0,0,a2_3U4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4_3,0,0,0,0,0,0,0,0,0,b5_2U3U4U5_3,b6_2U3U4U6_3,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3_3*x,0,a3_2U5_3*x,0,0,a2_3U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5_3,0,0,0,0,0,0,0,0,b4_2U3U4U5_3,0,b6_2U3U5U6_3,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3_3*x,0,0,a3_2U6_3*x,0,0,a2_3U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U6_3,0,0,0,0,0,0,0,0,b4_2U3U4U6_3,b5_2U3U5U6_3,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U4_3*x,a4_2U5_3*x,0,0,0,0,a2_4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5_3,0,0,0,0,0,0,b3_2U3U4U5_3,0,0,b6_2U4U5U6_3,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4_3*x,0,a4_2U6_3*x,0,0,0,0,a2_4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U6_3,0,0,0,0,0,0,b3_2U3U4U6_3,0,b5_2U4U5U6_3,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U5_3*x,a5_2U6_3*x,0,0,0,0,0,a2_5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U5U6_3,0,0,0,0,0,0,b3_2U3U5U6_3,b4_2U4U5U6_3,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_3U4_3*x,a4_3U5_3*x,0,a3_4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5_3,0,0,0,b2_2U3U4U5_3,0,0,0,b6_3U4U5U6_3,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4_3*x,0,a4_3U6_3*x,0,a3_4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U6_3,0,0,0,b2_2U3U4U6_3,0,0,b5_3U4U5U6_3,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U5_3*x,a5_3U6_3*x,0,0,a3_5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U5U6_3,0,0,0,b2_2U3U5U6_3,0,b4_3U4U5U6_3,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_4U5_3*x,a5_4U6_3*x,a4_5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U4U5U6_3,0,0,0,b2_2U4U5U6_3,b3_3U4U5U6_3,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a4_1U2U3_3*x,a3_1U2U4_3*x,0,0,a2_1U3U4_3*x,0,0,0,0,0,a1_2U3U4_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5_3,b6_1U2U3U4U6_3,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3_3*x,0,a3_1U2U5_3*x,0,0,a2_1U3U5_3*x,0,0,0,0,0,a1_2U3U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U5_3,0,b6_1U2U3U5U6_3,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3_3*x,0,0,a3_1U2U6_3*x,0,0,a2_1U3U6_3*x,0,0,0,0,0,a1_2U3U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b4_1U2U3U4U6_3,b5_1U2U3U5U6_3,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U4_3*x,a4_1U2U5_3*x,0,0,0,0,a2_1U4U5_3*x,0,0,0,0,0,a1_2U4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U5_3,0,0,b6_1U2U4U5U6_3,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4_3*x,0,a4_1U2U6_3*x,0,0,0,0,a2_1U4U6_3*x,0,0,0,0,0,a1_2U4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U4U6_3,0,b5_1U2U4U5U6_3,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U5_3*x,a5_1U2U6_3*x,0,0,0,0,0,a2_1U5U6_3*x,0,0,0,0,0,a1_2U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b3_1U2U3U5U6_3,b4_1U2U4U5U6_3,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U3U4_3*x,a4_1U3U5_3*x,0,a3_1U4U5_3*x,0,0,0,0,0,0,0,0,a1_3U4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U5_3,0,0,0,b6_1U3U4U5U6_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4_3*x,0,a4_1U3U6_3*x,0,a3_1U4U6_3*x,0,0,0,0,0,0,0,0,a1_3U4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U4U6_3,0,0,b5_1U3U4U5U6_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U5_3*x,a5_1U3U6_3*x,0,0,a3_1U5U6_3*x,0,0,0,0,0,0,0,0,a1_3U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U3U5U6_3,0,b4_1U3U4U5U6_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U4U5_3*x,a5_1U4U6_3*x,a4_1U5U6_3*x,0,0,0,0,0,0,0,0,0,a1_4U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b2_1U2U4U5U6_3,b3_1U3U4U5U6_3,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_2U3U4_3*x,a4_2U3U5_3*x,0,a3_2U4U5_3*x,0,0,a2_3U4U5_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U5_3,0,0,0,0,b6_2U3U4U5U6_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4_3*x,0,a4_2U3U6_3*x,0,a3_2U4U6_3*x,0,0,a2_3U4U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U4U6_3,0,0,0,b5_2U3U4U5U6_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U5_3*x,a5_2U3U6_3*x,0,0,a3_2U5U6_3*x,0,0,a2_3U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U3U5U6_3,0,0,b4_2U3U4U5U6_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U4U5_3*x,a5_2U4U6_3*x,a4_2U5U6_3*x,0,0,0,a2_4U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U2U4U5U6_3,0,b3_2U3U4U5U6_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_3U4U5_3*x,a5_3U4U6_3*x,a4_3U5U6_3*x,a3_4U5U6_3*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,b1_1U3U4U5U6_3,b2_2U3U4U5U6_3,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U4U5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a5_1U2U3U4_3*x,a4_1U2U3U5_3*x,0,a3_1U2U4U5_3*x,0,0,a2_1U3U4U5_3*x,0,0,0,a1_2U3U4U5_3*x,0,0,0,0,0,0,0,0,0,0,b6_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U4U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4_3*x,0,a4_1U2U3U6_3*x,0,a3_1U2U4U6_3*x,0,0,a2_1U3U4U6_3*x,0,0,0,a1_2U3U4U6_3*x,0,0,0,0,0,0,0,0,0,b5_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U5_3*x,a5_1U2U3U6_3*x,0,0,a3_1U2U5U6_3*x,0,0,a2_1U3U5U6_3*x,0,0,0,a1_2U3U5U6_3*x,0,0,0,0,0,0,0,0,b4_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U4U5_3*x,a5_1U2U4U6_3*x,a4_1U2U5U6_3*x,0,0,0,a2_1U4U5U6_3*x,0,0,0,a1_2U4U5U6_3*x,0,0,0,0,0,0,0,b3_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U3U4U5_3*x,a5_1U3U4U6_3*x,a4_1U3U5U6_3*x,a3_1U4U5U6_3*x,0,0,0,0,a1_3U4U5U6_3*x,0,0,0,0,0,0,b2_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_2U3U4U5_3*x,a5_2U3U4U6_3*x,a4_2U3U5U6_3*x,a3_2U4U5U6_3*x,a2_3U4U5U6_3*x,0,0,0,0,0,0,b1_1U2U3U4U5U6_3,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,k_2_1U2U3U4U5U6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,a6_1U2U3U4U5_3*x,a5_1U2U3U4U6_3*x,a4_1U2U3U5U6_3*x,a3_1U2U4U5U6_3*x,a2_1U3U4U5U6_3*x,a1_2U3U4U5U6_3*x,0;

    for (unsigned i = 0; i < 192; ++i){

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
    py::array_t<double> resultpy = py::array_t<double>(192);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    
    for (int i=0;i<192;i++){
        ptrresultpy[i]=rhos[i].template convert_to<double>();
    }
    return resultpy;}
double interfacess(py::array_t<double> parsar, py::array_t<double> varvals){
    std::vector<int> indicesC={128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191};
    std::vector<int> coeffsC={6,13,23,33,43,53,63,73,83,93,103,113,129,139,149,159,181,191,201,229,239,273,313,323,333,343,359,369,379,401,411,439,479,489,499,521,531,559,605,615,643,695,753,763,773,789,799,821,855,865,887,927,979,989,1011,1051,1109,1173,1183,1199,1227,1267,1319,1377};

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
    PYBIND11_MODULE(base_6sites_highprec,m){

    m.def("interfacerhos", &interfacerhos, "A function which returns the normalised ss rhos.",
    py::arg("parsar"), py::arg("varvals"));
    m.def("interfacess", &interfacess, "A function which returns the ss output, where appropriate rhos are multiplied by rates.",
    py::arg("parsar"), py::arg("varvals"));
}

    