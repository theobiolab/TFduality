
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
        long double A=vars[0];
    long double B=vars[1];
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
    long double a2_1_1=pars[24];
    long double a2_1_2=pars[25];
    long double a2_1_3=pars[26];
    long double b2_1U2_1=pars[27];
    long double k_1_1U2=pars[28];
    long double kr_1_1U2=pars[29];
    long double b2_1U2_2=pars[30];
    long double k_2_1U2=pars[31];
    long double b2_1U2_3=pars[32];
    long double k_3_1U2=pars[33];
    long double a1_2_1=pars[34];
    long double a1_2_2=pars[35];
    long double a1_2_3=pars[36];
    long double b1_1U2_1=pars[37];
    long double b1_1U2_2=pars[38];
    long double b1_1U2_3=pars[39];
    Matrix<InternalType, Dynamic, Dynamic> L = Matrix<InternalType, Dynamic, Dynamic>::Zero(12, 12);
    L<<0,b1_1_1,b2_2_1,0,kr_1_0,0,0,0,k_3_0,0,0,0,
a1_0_1*A,0,0,b2_1U2_1,0,kr_1_1,0,0,0,k_3_1,0,0,
a2_0_1*B,0,0,b1_1U2_1,0,0,kr_1_2,0,0,0,k_3_2,0,
0,a2_1_1*B,a1_2_1*A,0,0,0,0,kr_1_1U2,0,0,0,k_3_1U2,
k_1_0,0,0,0,0,b1_1_2,b2_2_2,0,0,0,0,0,
0,k_1_1,0,0,a1_0_2*A,0,0,b2_1U2_2,0,0,0,0,
0,0,k_1_2,0,a2_0_2*B,0,0,b1_1U2_2,0,0,0,0,
0,0,0,k_1_1U2,0,a2_1_2*B,a1_2_2*A,0,0,0,0,0,
0,0,0,0,k_2_0,0,0,0,0,b1_1_3,b2_2_3,0,
0,0,0,0,0,k_2_1,0,0,a1_0_3*A,0,0,b2_1U2_3,
0,0,0,0,0,0,k_2_2,0,a2_0_3*B,0,0,b1_1U2_3,
0,0,0,0,0,0,0,k_2_1U2,0,a2_1_3*B,a1_2_3*A,0;

    for (unsigned i = 0; i < 12; ++i){

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
    py::array_t<double> resultpy = py::array_t<double>(12);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    
    for (int i=0;i<12;i++){
        ptrresultpy[i]=rhos[i].template convert_to<double>();
    }
    return resultpy;}
double interfacess(py::array_t<double> parsar, py::array_t<double> varvals){
    std::vector<int> indicesC={8,9,10,11};
    std::vector<int> coeffsC={6,13,23,33};

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
    PYBIND11_MODULE(bPcycle_2bs_3states_rev1_A_B_prec50,m){

    m.def("interfacerhos", &interfacerhos, "A function which returns the normalised ss rhos.",
    py::arg("parsar"), py::arg("varvals"));
    m.def("interfacess", &interfacess, "A function which returns the ss output, where appropriate rhos are multiplied by rates.",
    py::arg("parsar"), py::arg("varvals"));
}

    