#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <boost/multiprecision/mpfr.hpp>

using namespace std;
using namespace Eigen;
namespace py=pybind11;

typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>> InternalType;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>> T; 
typedef Eigen::SparseMatrix<InternalType> SM;

void insert_L_Lx_atval(SM& L, const std::vector<Eigen::Triplet<InternalType>>& Lx, double xval){
    Eigen::Triplet<InternalType> trd;
    for (int j=0;j<Lx.size();j++){
        trd=Lx[j];
        L.insert(trd.row(),trd.col())=trd.value()*xval;
    }
}

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

void pre_laplacianpars(py::array_t<double> parsar, SM& L , std::vector<Eigen::Triplet<T>>& Lxx ){
//laplacian matrix with only parameter values, no variables. Also no values for those labels that depend upon variables.

    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    int n=48;

    T a1_0_1=pars[0];
    T k_1_0=pars[1];
    T kr_1_0=pars[2];
    T a1_0_2=pars[3];
    T k_2_0=pars[4];
    T a1_0_3=pars[5];
    T k_3_0=pars[6];
    T b1_1_1=pars[7];
    T k_1_1=pars[8];
    T kr_1_1=pars[9];
    T b1_1_2=pars[10];
    T k_2_1=pars[11];
    T b1_1_3=pars[12];
    T k_3_1=pars[13];
    T a2_0_1=pars[14];
    T a2_0_2=pars[15];
    T a2_0_3=pars[16];
    T b2_2_1=pars[17];
    T k_1_2=pars[18];
    T kr_1_2=pars[19];
    T b2_2_2=pars[20];
    T k_2_2=pars[21];
    T b2_2_3=pars[22];
    T k_3_2=pars[23];
    T a3_0_1=pars[24];
    T a3_0_2=pars[25];
    T a3_0_3=pars[26];
    T b3_3_1=pars[27];
    T k_1_3=pars[28];
    T kr_1_3=pars[29];
    T b3_3_2=pars[30];
    T k_2_3=pars[31];
    T b3_3_3=pars[32];
    T k_3_3=pars[33];
    T a4_0_1=pars[34];
    T a4_0_2=pars[35];
    T a4_0_3=pars[36];
    T b4_4_1=pars[37];
    T k_1_4=pars[38];
    T kr_1_4=pars[39];
    T b4_4_2=pars[40];
    T k_2_4=pars[41];
    T b4_4_3=pars[42];
    T k_3_4=pars[43];
    T a2_1_1=pars[44];
    T a2_1_2=pars[45];
    T a2_1_3=pars[46];
    T b2_1U2_1=pars[47];
    T k_1_1U2=pars[48];
    T kr_1_1U2=pars[49];
    T b2_1U2_2=pars[50];
    T k_2_1U2=pars[51];
    T b2_1U2_3=pars[52];
    T k_3_1U2=pars[53];
    T a3_1_1=pars[54];
    T a3_1_2=pars[55];
    T a3_1_3=pars[56];
    T b3_1U3_1=pars[57];
    T k_1_1U3=pars[58];
    T kr_1_1U3=pars[59];
    T b3_1U3_2=pars[60];
    T k_2_1U3=pars[61];
    T b3_1U3_3=pars[62];
    T k_3_1U3=pars[63];
    T a4_1_1=pars[64];
    T a4_1_2=pars[65];
    T a4_1_3=pars[66];
    T b4_1U4_1=pars[67];
    T k_1_1U4=pars[68];
    T kr_1_1U4=pars[69];
    T b4_1U4_2=pars[70];
    T k_2_1U4=pars[71];
    T b4_1U4_3=pars[72];
    T k_3_1U4=pars[73];
    T a1_2_1=pars[74];
    T a1_2_2=pars[75];
    T a1_2_3=pars[76];
    T b1_1U2_1=pars[77];
    T b1_1U2_2=pars[78];
    T b1_1U2_3=pars[79];
    T a3_2_1=pars[80];
    T a3_2_2=pars[81];
    T a3_2_3=pars[82];
    T b3_2U3_1=pars[83];
    T k_1_2U3=pars[84];
    T kr_1_2U3=pars[85];
    T b3_2U3_2=pars[86];
    T k_2_2U3=pars[87];
    T b3_2U3_3=pars[88];
    T k_3_2U3=pars[89];
    T a4_2_1=pars[90];
    T a4_2_2=pars[91];
    T a4_2_3=pars[92];
    T b4_2U4_1=pars[93];
    T k_1_2U4=pars[94];
    T kr_1_2U4=pars[95];
    T b4_2U4_2=pars[96];
    T k_2_2U4=pars[97];
    T b4_2U4_3=pars[98];
    T k_3_2U4=pars[99];
    T a1_3_1=pars[100];
    T a1_3_2=pars[101];
    T a1_3_3=pars[102];
    T b1_1U3_1=pars[103];
    T b1_1U3_2=pars[104];
    T b1_1U3_3=pars[105];
    T a2_3_1=pars[106];
    T a2_3_2=pars[107];
    T a2_3_3=pars[108];
    T b2_2U3_1=pars[109];
    T b2_2U3_2=pars[110];
    T b2_2U3_3=pars[111];
    T a4_3_1=pars[112];
    T a4_3_2=pars[113];
    T a4_3_3=pars[114];
    T b4_3U4_1=pars[115];
    T k_1_3U4=pars[116];
    T kr_1_3U4=pars[117];
    T b4_3U4_2=pars[118];
    T k_2_3U4=pars[119];
    T b4_3U4_3=pars[120];
    T k_3_3U4=pars[121];
    T a1_4_1=pars[122];
    T a1_4_2=pars[123];
    T a1_4_3=pars[124];
    T b1_1U4_1=pars[125];
    T b1_1U4_2=pars[126];
    T b1_1U4_3=pars[127];
    T a2_4_1=pars[128];
    T a2_4_2=pars[129];
    T a2_4_3=pars[130];
    T b2_2U4_1=pars[131];
    T b2_2U4_2=pars[132];
    T b2_2U4_3=pars[133];
    T a3_4_1=pars[134];
    T a3_4_2=pars[135];
    T a3_4_3=pars[136];
    T b3_3U4_1=pars[137];
    T b3_3U4_2=pars[138];
    T b3_3U4_3=pars[139];
    T a3_1U2_1=pars[140];
    T a3_1U2_2=pars[141];
    T a3_1U2_3=pars[142];
    T b3_1U2U3_1=pars[143];
    T k_1_1U2U3=pars[144];
    T kr_1_1U2U3=pars[145];
    T b3_1U2U3_2=pars[146];
    T k_2_1U2U3=pars[147];
    T b3_1U2U3_3=pars[148];
    T k_3_1U2U3=pars[149];
    T a4_1U2_1=pars[150];
    T a4_1U2_2=pars[151];
    T a4_1U2_3=pars[152];
    T b4_1U2U4_1=pars[153];
    T k_1_1U2U4=pars[154];
    T kr_1_1U2U4=pars[155];
    T b4_1U2U4_2=pars[156];
    T k_2_1U2U4=pars[157];
    T b4_1U2U4_3=pars[158];
    T k_3_1U2U4=pars[159];
    T a2_1U3_1=pars[160];
    T a2_1U3_2=pars[161];
    T a2_1U3_3=pars[162];
    T b2_1U2U3_1=pars[163];
    T b2_1U2U3_2=pars[164];
    T b2_1U2U3_3=pars[165];
    T a4_1U3_1=pars[166];
    T a4_1U3_2=pars[167];
    T a4_1U3_3=pars[168];
    T b4_1U3U4_1=pars[169];
    T k_1_1U3U4=pars[170];
    T kr_1_1U3U4=pars[171];
    T b4_1U3U4_2=pars[172];
    T k_2_1U3U4=pars[173];
    T b4_1U3U4_3=pars[174];
    T k_3_1U3U4=pars[175];
    T a2_1U4_1=pars[176];
    T a2_1U4_2=pars[177];
    T a2_1U4_3=pars[178];
    T b2_1U2U4_1=pars[179];
    T b2_1U2U4_2=pars[180];
    T b2_1U2U4_3=pars[181];
    T a3_1U4_1=pars[182];
    T a3_1U4_2=pars[183];
    T a3_1U4_3=pars[184];
    T b3_1U3U4_1=pars[185];
    T b3_1U3U4_2=pars[186];
    T b3_1U3U4_3=pars[187];
    T a1_2U3_1=pars[188];
    T a1_2U3_2=pars[189];
    T a1_2U3_3=pars[190];
    T b1_1U2U3_1=pars[191];
    T b1_1U2U3_2=pars[192];
    T b1_1U2U3_3=pars[193];
    T a4_2U3_1=pars[194];
    T a4_2U3_2=pars[195];
    T a4_2U3_3=pars[196];
    T b4_2U3U4_1=pars[197];
    T k_1_2U3U4=pars[198];
    T kr_1_2U3U4=pars[199];
    T b4_2U3U4_2=pars[200];
    T k_2_2U3U4=pars[201];
    T b4_2U3U4_3=pars[202];
    T k_3_2U3U4=pars[203];
    T a1_2U4_1=pars[204];
    T a1_2U4_2=pars[205];
    T a1_2U4_3=pars[206];
    T b1_1U2U4_1=pars[207];
    T b1_1U2U4_2=pars[208];
    T b1_1U2U4_3=pars[209];
    T a3_2U4_1=pars[210];
    T a3_2U4_2=pars[211];
    T a3_2U4_3=pars[212];
    T b3_2U3U4_1=pars[213];
    T b3_2U3U4_2=pars[214];
    T b3_2U3U4_3=pars[215];
    T a1_3U4_1=pars[216];
    T a1_3U4_2=pars[217];
    T a1_3U4_3=pars[218];
    T b1_1U3U4_1=pars[219];
    T b1_1U3U4_2=pars[220];
    T b1_1U3U4_3=pars[221];
    T a2_3U4_1=pars[222];
    T a2_3U4_2=pars[223];
    T a2_3U4_3=pars[224];
    T b2_2U3U4_1=pars[225];
    T b2_2U3U4_2=pars[226];
    T b2_2U3U4_3=pars[227];
    T a4_1U2U3_1=pars[228];
    T a4_1U2U3_2=pars[229];
    T a4_1U2U3_3=pars[230];
    T b4_1U2U3U4_1=pars[231];
    T k_1_1U2U3U4=pars[232];
    T kr_1_1U2U3U4=pars[233];
    T b4_1U2U3U4_2=pars[234];
    T k_2_1U2U3U4=pars[235];
    T b4_1U2U3U4_3=pars[236];
    T k_3_1U2U3U4=pars[237];
    T a3_1U2U4_1=pars[238];
    T a3_1U2U4_2=pars[239];
    T a3_1U2U4_3=pars[240];
    T b3_1U2U3U4_1=pars[241];
    T b3_1U2U3U4_2=pars[242];
    T b3_1U2U3U4_3=pars[243];
    T a2_1U3U4_1=pars[244];
    T a2_1U3U4_2=pars[245];
    T a2_1U3U4_3=pars[246];
    T b2_1U2U3U4_1=pars[247];
    T b2_1U2U3U4_2=pars[248];
    T b2_1U2U3U4_3=pars[249];
    T a1_2U3U4_1=pars[250];
    T a1_2U3U4_2=pars[251];
    T a1_2U3U4_3=pars[252];
    T b1_1U2U3U4_1=pars[253];
    T b1_1U2U3U4_2=pars[254];
    T b1_1U2U3U4_3=pars[255];
    std::vector<Eigen::Triplet<T>> clist;
    clist.push_back(Eigen::Triplet<T>(0,1,b1_1_1));
    clist.push_back(Eigen::Triplet<T>(0,2,b2_2_1));
    clist.push_back(Eigen::Triplet<T>(0,3,b3_3_1));
    clist.push_back(Eigen::Triplet<T>(0,4,b4_4_1));
    clist.push_back(Eigen::Triplet<T>(0,16,kr_1_0));
    clist.push_back(Eigen::Triplet<T>(0,32,k_3_0));
    clist.push_back(Eigen::Triplet<T>(1,5,b2_1U2_1));
    clist.push_back(Eigen::Triplet<T>(1,6,b3_1U3_1));
    clist.push_back(Eigen::Triplet<T>(1,7,b4_1U4_1));
    clist.push_back(Eigen::Triplet<T>(1,17,kr_1_1));
    clist.push_back(Eigen::Triplet<T>(1,33,k_3_1));
    clist.push_back(Eigen::Triplet<T>(2,5,b1_1U2_1));
    clist.push_back(Eigen::Triplet<T>(2,8,b3_2U3_1));
    clist.push_back(Eigen::Triplet<T>(2,9,b4_2U4_1));
    clist.push_back(Eigen::Triplet<T>(2,18,kr_1_2));
    clist.push_back(Eigen::Triplet<T>(2,34,k_3_2));
    clist.push_back(Eigen::Triplet<T>(3,6,b1_1U3_1));
    clist.push_back(Eigen::Triplet<T>(3,8,b2_2U3_1));
    clist.push_back(Eigen::Triplet<T>(3,10,b4_3U4_1));
    clist.push_back(Eigen::Triplet<T>(3,19,kr_1_3));
    clist.push_back(Eigen::Triplet<T>(3,35,k_3_3));
    clist.push_back(Eigen::Triplet<T>(4,7,b1_1U4_1));
    clist.push_back(Eigen::Triplet<T>(4,9,b2_2U4_1));
    clist.push_back(Eigen::Triplet<T>(4,10,b3_3U4_1));
    clist.push_back(Eigen::Triplet<T>(4,20,kr_1_4));
    clist.push_back(Eigen::Triplet<T>(4,36,k_3_4));
    clist.push_back(Eigen::Triplet<T>(5,11,b3_1U2U3_1));
    clist.push_back(Eigen::Triplet<T>(5,12,b4_1U2U4_1));
    clist.push_back(Eigen::Triplet<T>(5,21,kr_1_1U2));
    clist.push_back(Eigen::Triplet<T>(5,37,k_3_1U2));
    clist.push_back(Eigen::Triplet<T>(6,11,b2_1U2U3_1));
    clist.push_back(Eigen::Triplet<T>(6,13,b4_1U3U4_1));
    clist.push_back(Eigen::Triplet<T>(6,22,kr_1_1U3));
    clist.push_back(Eigen::Triplet<T>(6,38,k_3_1U3));
    clist.push_back(Eigen::Triplet<T>(7,12,b2_1U2U4_1));
    clist.push_back(Eigen::Triplet<T>(7,13,b3_1U3U4_1));
    clist.push_back(Eigen::Triplet<T>(7,23,kr_1_1U4));
    clist.push_back(Eigen::Triplet<T>(7,39,k_3_1U4));
    clist.push_back(Eigen::Triplet<T>(8,11,b1_1U2U3_1));
    clist.push_back(Eigen::Triplet<T>(8,14,b4_2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(8,24,kr_1_2U3));
    clist.push_back(Eigen::Triplet<T>(8,40,k_3_2U3));
    clist.push_back(Eigen::Triplet<T>(9,12,b1_1U2U4_1));
    clist.push_back(Eigen::Triplet<T>(9,14,b3_2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(9,25,kr_1_2U4));
    clist.push_back(Eigen::Triplet<T>(9,41,k_3_2U4));
    clist.push_back(Eigen::Triplet<T>(10,13,b1_1U3U4_1));
    clist.push_back(Eigen::Triplet<T>(10,14,b2_2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(10,26,kr_1_3U4));
    clist.push_back(Eigen::Triplet<T>(10,42,k_3_3U4));
    clist.push_back(Eigen::Triplet<T>(11,15,b4_1U2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(11,27,kr_1_1U2U3));
    clist.push_back(Eigen::Triplet<T>(11,43,k_3_1U2U3));
    clist.push_back(Eigen::Triplet<T>(12,15,b3_1U2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(12,28,kr_1_1U2U4));
    clist.push_back(Eigen::Triplet<T>(12,44,k_3_1U2U4));
    clist.push_back(Eigen::Triplet<T>(13,15,b2_1U2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(13,29,kr_1_1U3U4));
    clist.push_back(Eigen::Triplet<T>(13,45,k_3_1U3U4));
    clist.push_back(Eigen::Triplet<T>(14,15,b1_1U2U3U4_1));
    clist.push_back(Eigen::Triplet<T>(14,30,kr_1_2U3U4));
    clist.push_back(Eigen::Triplet<T>(14,46,k_3_2U3U4));
    clist.push_back(Eigen::Triplet<T>(15,31,kr_1_1U2U3U4));
    clist.push_back(Eigen::Triplet<T>(15,47,k_3_1U2U3U4));
    clist.push_back(Eigen::Triplet<T>(16,0,k_1_0));
    clist.push_back(Eigen::Triplet<T>(16,17,b1_1_2));
    clist.push_back(Eigen::Triplet<T>(16,18,b2_2_2));
    clist.push_back(Eigen::Triplet<T>(16,19,b3_3_2));
    clist.push_back(Eigen::Triplet<T>(16,20,b4_4_2));
    clist.push_back(Eigen::Triplet<T>(17,1,k_1_1));
    clist.push_back(Eigen::Triplet<T>(17,21,b2_1U2_2));
    clist.push_back(Eigen::Triplet<T>(17,22,b3_1U3_2));
    clist.push_back(Eigen::Triplet<T>(17,23,b4_1U4_2));
    clist.push_back(Eigen::Triplet<T>(18,2,k_1_2));
    clist.push_back(Eigen::Triplet<T>(18,21,b1_1U2_2));
    clist.push_back(Eigen::Triplet<T>(18,24,b3_2U3_2));
    clist.push_back(Eigen::Triplet<T>(18,25,b4_2U4_2));
    clist.push_back(Eigen::Triplet<T>(19,3,k_1_3));
    clist.push_back(Eigen::Triplet<T>(19,22,b1_1U3_2));
    clist.push_back(Eigen::Triplet<T>(19,24,b2_2U3_2));
    clist.push_back(Eigen::Triplet<T>(19,26,b4_3U4_2));
    clist.push_back(Eigen::Triplet<T>(20,4,k_1_4));
    clist.push_back(Eigen::Triplet<T>(20,23,b1_1U4_2));
    clist.push_back(Eigen::Triplet<T>(20,25,b2_2U4_2));
    clist.push_back(Eigen::Triplet<T>(20,26,b3_3U4_2));
    clist.push_back(Eigen::Triplet<T>(21,5,k_1_1U2));
    clist.push_back(Eigen::Triplet<T>(21,27,b3_1U2U3_2));
    clist.push_back(Eigen::Triplet<T>(21,28,b4_1U2U4_2));
    clist.push_back(Eigen::Triplet<T>(22,6,k_1_1U3));
    clist.push_back(Eigen::Triplet<T>(22,27,b2_1U2U3_2));
    clist.push_back(Eigen::Triplet<T>(22,29,b4_1U3U4_2));
    clist.push_back(Eigen::Triplet<T>(23,7,k_1_1U4));
    clist.push_back(Eigen::Triplet<T>(23,28,b2_1U2U4_2));
    clist.push_back(Eigen::Triplet<T>(23,29,b3_1U3U4_2));
    clist.push_back(Eigen::Triplet<T>(24,8,k_1_2U3));
    clist.push_back(Eigen::Triplet<T>(24,27,b1_1U2U3_2));
    clist.push_back(Eigen::Triplet<T>(24,30,b4_2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(25,9,k_1_2U4));
    clist.push_back(Eigen::Triplet<T>(25,28,b1_1U2U4_2));
    clist.push_back(Eigen::Triplet<T>(25,30,b3_2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(26,10,k_1_3U4));
    clist.push_back(Eigen::Triplet<T>(26,29,b1_1U3U4_2));
    clist.push_back(Eigen::Triplet<T>(26,30,b2_2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(27,11,k_1_1U2U3));
    clist.push_back(Eigen::Triplet<T>(27,31,b4_1U2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(28,12,k_1_1U2U4));
    clist.push_back(Eigen::Triplet<T>(28,31,b3_1U2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(29,13,k_1_1U3U4));
    clist.push_back(Eigen::Triplet<T>(29,31,b2_1U2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(30,14,k_1_2U3U4));
    clist.push_back(Eigen::Triplet<T>(30,31,b1_1U2U3U4_2));
    clist.push_back(Eigen::Triplet<T>(31,15,k_1_1U2U3U4));
    clist.push_back(Eigen::Triplet<T>(32,16,k_2_0));
    clist.push_back(Eigen::Triplet<T>(32,33,b1_1_3));
    clist.push_back(Eigen::Triplet<T>(32,34,b2_2_3));
    clist.push_back(Eigen::Triplet<T>(32,35,b3_3_3));
    clist.push_back(Eigen::Triplet<T>(32,36,b4_4_3));
    clist.push_back(Eigen::Triplet<T>(33,17,k_2_1));
    clist.push_back(Eigen::Triplet<T>(33,37,b2_1U2_3));
    clist.push_back(Eigen::Triplet<T>(33,38,b3_1U3_3));
    clist.push_back(Eigen::Triplet<T>(33,39,b4_1U4_3));
    clist.push_back(Eigen::Triplet<T>(34,18,k_2_2));
    clist.push_back(Eigen::Triplet<T>(34,37,b1_1U2_3));
    clist.push_back(Eigen::Triplet<T>(34,40,b3_2U3_3));
    clist.push_back(Eigen::Triplet<T>(34,41,b4_2U4_3));
    clist.push_back(Eigen::Triplet<T>(35,19,k_2_3));
    clist.push_back(Eigen::Triplet<T>(35,38,b1_1U3_3));
    clist.push_back(Eigen::Triplet<T>(35,40,b2_2U3_3));
    clist.push_back(Eigen::Triplet<T>(35,42,b4_3U4_3));
    clist.push_back(Eigen::Triplet<T>(36,20,k_2_4));
    clist.push_back(Eigen::Triplet<T>(36,39,b1_1U4_3));
    clist.push_back(Eigen::Triplet<T>(36,41,b2_2U4_3));
    clist.push_back(Eigen::Triplet<T>(36,42,b3_3U4_3));
    clist.push_back(Eigen::Triplet<T>(37,21,k_2_1U2));
    clist.push_back(Eigen::Triplet<T>(37,43,b3_1U2U3_3));
    clist.push_back(Eigen::Triplet<T>(37,44,b4_1U2U4_3));
    clist.push_back(Eigen::Triplet<T>(38,22,k_2_1U3));
    clist.push_back(Eigen::Triplet<T>(38,43,b2_1U2U3_3));
    clist.push_back(Eigen::Triplet<T>(38,45,b4_1U3U4_3));
    clist.push_back(Eigen::Triplet<T>(39,23,k_2_1U4));
    clist.push_back(Eigen::Triplet<T>(39,44,b2_1U2U4_3));
    clist.push_back(Eigen::Triplet<T>(39,45,b3_1U3U4_3));
    clist.push_back(Eigen::Triplet<T>(40,24,k_2_2U3));
    clist.push_back(Eigen::Triplet<T>(40,43,b1_1U2U3_3));
    clist.push_back(Eigen::Triplet<T>(40,46,b4_2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(41,25,k_2_2U4));
    clist.push_back(Eigen::Triplet<T>(41,44,b1_1U2U4_3));
    clist.push_back(Eigen::Triplet<T>(41,46,b3_2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(42,26,k_2_3U4));
    clist.push_back(Eigen::Triplet<T>(42,45,b1_1U3U4_3));
    clist.push_back(Eigen::Triplet<T>(42,46,b2_2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(43,27,k_2_1U2U3));
    clist.push_back(Eigen::Triplet<T>(43,47,b4_1U2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(44,28,k_2_1U2U4));
    clist.push_back(Eigen::Triplet<T>(44,47,b3_1U2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(45,29,k_2_1U3U4));
    clist.push_back(Eigen::Triplet<T>(45,47,b2_1U2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(46,30,k_2_2U3U4));
    clist.push_back(Eigen::Triplet<T>(46,47,b1_1U2U3U4_3));
    clist.push_back(Eigen::Triplet<T>(47,31,k_2_1U2U3U4));

    Eigen::Triplet<T> trd;
    for (int j=0;j<clist.size();j++){
        trd=clist[j];
        L.insert(trd.row(),trd.col())=trd.value();
    }
    L.makeCompressed();
    Lxx.push_back(Eigen::Triplet<T>(1,0,a1_0_1));
    Lxx.push_back(Eigen::Triplet<T>(2,0,a2_0_1));
    Lxx.push_back(Eigen::Triplet<T>(3,0,a3_0_1));
    Lxx.push_back(Eigen::Triplet<T>(4,0,a4_0_1));
    Lxx.push_back(Eigen::Triplet<T>(5,1,a2_1_1));
    Lxx.push_back(Eigen::Triplet<T>(5,2,a1_2_1));
    Lxx.push_back(Eigen::Triplet<T>(6,1,a3_1_1));
    Lxx.push_back(Eigen::Triplet<T>(6,3,a1_3_1));
    Lxx.push_back(Eigen::Triplet<T>(7,1,a4_1_1));
    Lxx.push_back(Eigen::Triplet<T>(7,4,a1_4_1));
    Lxx.push_back(Eigen::Triplet<T>(8,2,a3_2_1));
    Lxx.push_back(Eigen::Triplet<T>(8,3,a2_3_1));
    Lxx.push_back(Eigen::Triplet<T>(9,2,a4_2_1));
    Lxx.push_back(Eigen::Triplet<T>(9,4,a2_4_1));
    Lxx.push_back(Eigen::Triplet<T>(10,3,a4_3_1));
    Lxx.push_back(Eigen::Triplet<T>(10,4,a3_4_1));
    Lxx.push_back(Eigen::Triplet<T>(11,5,a3_1U2_1));
    Lxx.push_back(Eigen::Triplet<T>(11,6,a2_1U3_1));
    Lxx.push_back(Eigen::Triplet<T>(11,8,a1_2U3_1));
    Lxx.push_back(Eigen::Triplet<T>(12,5,a4_1U2_1));
    Lxx.push_back(Eigen::Triplet<T>(12,7,a2_1U4_1));
    Lxx.push_back(Eigen::Triplet<T>(12,9,a1_2U4_1));
    Lxx.push_back(Eigen::Triplet<T>(13,6,a4_1U3_1));
    Lxx.push_back(Eigen::Triplet<T>(13,7,a3_1U4_1));
    Lxx.push_back(Eigen::Triplet<T>(13,10,a1_3U4_1));
    Lxx.push_back(Eigen::Triplet<T>(14,8,a4_2U3_1));
    Lxx.push_back(Eigen::Triplet<T>(14,9,a3_2U4_1));
    Lxx.push_back(Eigen::Triplet<T>(14,10,a2_3U4_1));
    Lxx.push_back(Eigen::Triplet<T>(15,11,a4_1U2U3_1));
    Lxx.push_back(Eigen::Triplet<T>(15,12,a3_1U2U4_1));
    Lxx.push_back(Eigen::Triplet<T>(15,13,a2_1U3U4_1));
    Lxx.push_back(Eigen::Triplet<T>(15,14,a1_2U3U4_1));
    Lxx.push_back(Eigen::Triplet<T>(17,16,a1_0_2));
    Lxx.push_back(Eigen::Triplet<T>(18,16,a2_0_2));
    Lxx.push_back(Eigen::Triplet<T>(19,16,a3_0_2));
    Lxx.push_back(Eigen::Triplet<T>(20,16,a4_0_2));
    Lxx.push_back(Eigen::Triplet<T>(21,17,a2_1_2));
    Lxx.push_back(Eigen::Triplet<T>(21,18,a1_2_2));
    Lxx.push_back(Eigen::Triplet<T>(22,17,a3_1_2));
    Lxx.push_back(Eigen::Triplet<T>(22,19,a1_3_2));
    Lxx.push_back(Eigen::Triplet<T>(23,17,a4_1_2));
    Lxx.push_back(Eigen::Triplet<T>(23,20,a1_4_2));
    Lxx.push_back(Eigen::Triplet<T>(24,18,a3_2_2));
    Lxx.push_back(Eigen::Triplet<T>(24,19,a2_3_2));
    Lxx.push_back(Eigen::Triplet<T>(25,18,a4_2_2));
    Lxx.push_back(Eigen::Triplet<T>(25,20,a2_4_2));
    Lxx.push_back(Eigen::Triplet<T>(26,19,a4_3_2));
    Lxx.push_back(Eigen::Triplet<T>(26,20,a3_4_2));
    Lxx.push_back(Eigen::Triplet<T>(27,21,a3_1U2_2));
    Lxx.push_back(Eigen::Triplet<T>(27,22,a2_1U3_2));
    Lxx.push_back(Eigen::Triplet<T>(27,24,a1_2U3_2));
    Lxx.push_back(Eigen::Triplet<T>(28,21,a4_1U2_2));
    Lxx.push_back(Eigen::Triplet<T>(28,23,a2_1U4_2));
    Lxx.push_back(Eigen::Triplet<T>(28,25,a1_2U4_2));
    Lxx.push_back(Eigen::Triplet<T>(29,22,a4_1U3_2));
    Lxx.push_back(Eigen::Triplet<T>(29,23,a3_1U4_2));
    Lxx.push_back(Eigen::Triplet<T>(29,26,a1_3U4_2));
    Lxx.push_back(Eigen::Triplet<T>(30,24,a4_2U3_2));
    Lxx.push_back(Eigen::Triplet<T>(30,25,a3_2U4_2));
    Lxx.push_back(Eigen::Triplet<T>(30,26,a2_3U4_2));
    Lxx.push_back(Eigen::Triplet<T>(31,27,a4_1U2U3_2));
    Lxx.push_back(Eigen::Triplet<T>(31,28,a3_1U2U4_2));
    Lxx.push_back(Eigen::Triplet<T>(31,29,a2_1U3U4_2));
    Lxx.push_back(Eigen::Triplet<T>(31,30,a1_2U3U4_2));
    Lxx.push_back(Eigen::Triplet<T>(33,32,a1_0_3));
    Lxx.push_back(Eigen::Triplet<T>(34,32,a2_0_3));
    Lxx.push_back(Eigen::Triplet<T>(35,32,a3_0_3));
    Lxx.push_back(Eigen::Triplet<T>(36,32,a4_0_3));
    Lxx.push_back(Eigen::Triplet<T>(37,33,a2_1_3));
    Lxx.push_back(Eigen::Triplet<T>(37,34,a1_2_3));
    Lxx.push_back(Eigen::Triplet<T>(38,33,a3_1_3));
    Lxx.push_back(Eigen::Triplet<T>(38,35,a1_3_3));
    Lxx.push_back(Eigen::Triplet<T>(39,33,a4_1_3));
    Lxx.push_back(Eigen::Triplet<T>(39,36,a1_4_3));
    Lxx.push_back(Eigen::Triplet<T>(40,34,a3_2_3));
    Lxx.push_back(Eigen::Triplet<T>(40,35,a2_3_3));
    Lxx.push_back(Eigen::Triplet<T>(41,34,a4_2_3));
    Lxx.push_back(Eigen::Triplet<T>(41,36,a2_4_3));
    Lxx.push_back(Eigen::Triplet<T>(42,35,a4_3_3));
    Lxx.push_back(Eigen::Triplet<T>(42,36,a3_4_3));
    Lxx.push_back(Eigen::Triplet<T>(43,37,a3_1U2_3));
    Lxx.push_back(Eigen::Triplet<T>(43,38,a2_1U3_3));
    Lxx.push_back(Eigen::Triplet<T>(43,40,a1_2U3_3));
    Lxx.push_back(Eigen::Triplet<T>(44,37,a4_1U2_3));
    Lxx.push_back(Eigen::Triplet<T>(44,39,a2_1U4_3));
    Lxx.push_back(Eigen::Triplet<T>(44,41,a1_2U4_3));
    Lxx.push_back(Eigen::Triplet<T>(45,38,a4_1U3_3));
    Lxx.push_back(Eigen::Triplet<T>(45,39,a3_1U4_3));
    Lxx.push_back(Eigen::Triplet<T>(45,42,a1_3U4_3));
    Lxx.push_back(Eigen::Triplet<T>(46,40,a4_2U3_3));
    Lxx.push_back(Eigen::Triplet<T>(46,41,a3_2U4_3));
    Lxx.push_back(Eigen::Triplet<T>(46,42,a2_3U4_3));
    Lxx.push_back(Eigen::Triplet<T>(47,43,a4_1U2U3_3));
    Lxx.push_back(Eigen::Triplet<T>(47,44,a3_1U2U4_3));
    Lxx.push_back(Eigen::Triplet<T>(47,45,a2_1U3U4_3));
    Lxx.push_back(Eigen::Triplet<T>(47,46,a1_2U3U4_3));
    return;
}


double interfacess(py::array_t<double> parsar, py::array_t<double> varvals, double tolerance){
    int i;
    const int n=48;
    SM L(n,n);
    L.reserve(VectorXi::Constant(n,5));    auto varsbuf=varvals.request();
    double *vars=(double *) varsbuf.ptr;
    double xval=vars[0];
    std::vector<Eigen::Triplet<T>> Lxx;
    pre_laplacianpars(parsar,L , Lxx);
    insert_L_Lx_atval(L, Lxx, xval);

    T cs;
    for (int k=0; k<L.outerSize(); ++k) {
        cs=0;
        for(typename Eigen::SparseMatrix<T>::InnerIterator it (L,k); it; ++it){
            cs+=it.value();
        }
        L.insert(k,k)=-cs;
        }
    double ssval=0.0;
    vector<int>indicesC={32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    vector<double>coeffsC={pars[6],pars[13],pars[23],pars[33],pars[43],pars[53],pars[63],pars[73],pars[89],pars[99],pars[121],pars[149],pars[159],pars[175],pars[203],pars[237]};
    Eigen::Matrix< InternalType, Eigen::Dynamic, Eigen::Dynamic > Ld=Eigen::Matrix< InternalType, Eigen::Dynamic, Eigen::Dynamic >(L);
    Matrix<InternalType, Dynamic, 1> steady_state;
    try
    {
    steady_state = getOneDimNullspaceFromSVD<InternalType>(Ld);
    //now doublecheck quality of nullspace
        Eigen::Matrix< InternalType, Eigen::Dynamic, Eigen::Dynamic > out(n,1);
        out=Ld*steady_state;
        i=0;
        bool goodnullspace=true;

        for (i=0;i<n;i++){
            if (abs(out(i,0))>tolerance){
                cout << "inaccurate nullspace to tolerance" << tolerance << "\n";
                cout << out;
                goodnullspace=false;
                break;
            }
        }
        cout << "Goodnullspace? " << goodnullspace << "\n";
    }
    catch (const std::runtime_error& e)
    {
        throw;
    }
    InternalType norm = steady_state.sum();
    
    for (i=0; i<steady_state.size();i++){
        steady_state[i]=steady_state[i] / norm;
    }
    for (i=0;i<indicesC.size();i++){
        ssval+=(steady_state[indicesC[i]]*coeffsC[i]).template convert_to<double>();
    }

    return  ssval;
    }



PYBIND11_MODULE(final_4sites_highprec,m){

       m.def("interfacess", &interfacess, "A function which returns ss.",
            py::arg("parsar"), py::arg("varvals"),py::arg("tolerance")=1e-10);

}


