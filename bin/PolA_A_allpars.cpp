#include "commonincludes.hpp"
template <typename T>
void GRF_PolA_A(py::array_t<double> parsar, vector<T> &num, vector<T> &den, py::array_t<double>othervars){

    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    T ktia0=pars[0];
    T ktan0=pars[1];
    T ktin0=pars[2];
    T ktni0=pars[3];
    T ktiaA=pars[4];
    T ktanA=pars[5];
    T ktinA=pars[6];
    T ktniA=pars[7];
    T kbAa=pars[8];
    T kuAa=pars[9];
    T kbAi=pars[10];
    T kuAi=pars[11];
    T kbAn=pars[12];
    T kuAn=pars[13];

    auto varsarbuf=othervars.request();
    double *varsar=(double *) varsarbuf.ptr;
    vector<T> coeffs_1={{ktanA*ktia0*ktiaA*ktni0*kuAn + ktanA*ktia0*ktinA*ktni0*kuAn + ktanA*ktia0*ktni0*ktniA*kuAi + ktanA*ktia0*ktni0*kuAi*kuAn + ktia0*ktiaA*ktni0*ktniA*kuAa + ktia0*ktiaA*ktni0*kuAa*kuAn + ktia0*ktinA*ktni0*kuAa*kuAn + ktia0*ktni0*ktniA*kuAa*kuAi + ktia0*ktni0*kuAa*kuAi*kuAn, kbAi*ktiaA*ktni0*ktniA*kuAa + kbAi*ktiaA*ktni0*kuAa*kuAn + kbAn*ktanA*ktia0*ktniA*kuAi + kbAn*ktia0*ktiaA*ktniA*kuAa + kbAn*ktia0*ktniA*kuAa*kuAi + kbAn*ktiaA*ktin0*ktniA*kuAa, kbAi*kbAn*ktiaA*ktniA*kuAa}};
    vector<T> coeffs_2={{ktan0*ktanA*ktiaA*ktni0*kuAn + ktan0*ktanA*ktinA*ktni0*kuAn + ktan0*ktanA*ktni0*ktniA*kuAi + ktan0*ktanA*ktni0*kuAi*kuAn + ktan0*ktiaA*ktni0*ktniA*kuAa + ktan0*ktiaA*ktni0*kuAa*kuAn + ktan0*ktinA*ktni0*kuAa*kuAn + ktan0*ktni0*ktniA*kuAa*kuAi + ktan0*ktni0*kuAa*kuAi*kuAn, kbAa*ktanA*ktiaA*ktni0*kuAn + kbAa*ktanA*ktinA*ktni0*kuAn + kbAa*ktanA*ktni0*ktniA*kuAi + kbAa*ktanA*ktni0*kuAi*kuAn + kbAn*ktan0*ktanA*ktniA*kuAi + kbAn*ktan0*ktniA*kuAa*kuAi, kbAa*kbAn*ktanA*ktniA*kuAi}};
    vector<T> coeffs_3={{ktan0*ktanA*ktia0*ktiaA*kuAn + ktan0*ktanA*ktia0*ktinA*kuAn + ktan0*ktanA*ktia0*ktniA*kuAi + ktan0*ktanA*ktia0*kuAi*kuAn + ktan0*ktanA*ktiaA*ktin0*kuAn + ktan0*ktanA*ktin0*ktinA*kuAn + ktan0*ktanA*ktin0*ktniA*kuAi + ktan0*ktanA*ktin0*kuAi*kuAn + ktan0*ktia0*ktiaA*ktniA*kuAa + ktan0*ktia0*ktiaA*kuAa*kuAn + ktan0*ktia0*ktinA*kuAa*kuAn + ktan0*ktia0*ktniA*kuAa*kuAi + ktan0*ktia0*kuAa*kuAi*kuAn + ktan0*ktiaA*ktin0*ktniA*kuAa + ktan0*ktiaA*ktin0*kuAa*kuAn + ktan0*ktin0*ktinA*kuAa*kuAn + ktan0*ktin0*ktniA*kuAa*kuAi + ktan0*ktin0*kuAa*kuAi*kuAn, kbAa*ktanA*ktia0*ktiaA*kuAn + kbAa*ktanA*ktia0*ktinA*kuAn + kbAa*ktanA*ktia0*kuAi*kuAn + kbAa*ktanA*ktiaA*ktin0*kuAn + kbAa*ktanA*ktin0*ktinA*kuAn + kbAa*ktanA*ktin0*ktniA*kuAi + kbAa*ktanA*ktin0*kuAi*kuAn + kbAi*ktan0*ktanA*ktiaA*kuAn + kbAi*ktan0*ktanA*ktinA*kuAn + kbAi*ktan0*ktiaA*ktniA*kuAa + kbAi*ktan0*ktiaA*kuAa*kuAn + kbAi*ktan0*ktinA*kuAa*kuAn, kbAa*kbAi*ktanA*ktiaA*kuAn + kbAa*kbAi*ktanA*ktinA*kuAn}};
    vector<T> coeffs_4={{0, kbAa*ktia0*ktiaA*ktni0*ktniA + kbAa*ktia0*ktiaA*ktni0*kuAn + kbAa*ktia0*ktinA*ktni0*kuAn + kbAa*ktia0*ktni0*ktniA*kuAi + kbAa*ktia0*ktni0*kuAi*kuAn + kbAi*ktan0*ktiaA*ktni0*ktniA + kbAi*ktan0*ktiaA*ktni0*kuAn + kbAn*ktan0*ktia0*ktiaA*ktniA + kbAn*ktan0*ktiaA*ktin0*ktniA, kbAa*kbAi*ktiaA*ktni0*ktniA + kbAa*kbAi*ktiaA*ktni0*kuAn + kbAa*kbAn*ktia0*ktiaA*ktniA + kbAa*kbAn*ktia0*ktniA*kuAi + kbAa*kbAn*ktiaA*ktin0*ktniA + kbAi*kbAn*ktan0*ktiaA*ktniA, kbAa*kbAi*kbAn*ktiaA*ktniA}};
    vector<T> coeffs_5={{0, kbAa*ktanA*ktia0*ktni0*ktniA + kbAi*ktan0*ktanA*ktni0*ktniA + kbAi*ktan0*ktanA*ktni0*kuAn + kbAi*ktan0*ktni0*ktniA*kuAa + kbAi*ktan0*ktni0*kuAa*kuAn + kbAn*ktan0*ktanA*ktia0*ktniA + kbAn*ktan0*ktanA*ktin0*ktniA + kbAn*ktan0*ktia0*ktniA*kuAa + kbAn*ktan0*ktin0*ktniA*kuAa, kbAa*kbAi*ktanA*ktni0*ktniA + kbAa*kbAi*ktanA*ktni0*kuAn + kbAa*kbAn*ktanA*ktia0*ktniA + kbAa*kbAn*ktanA*ktin0*ktniA + kbAi*kbAn*ktan0*ktanA*ktniA + kbAi*kbAn*ktan0*ktniA*kuAa, kbAa*kbAi*kbAn*ktanA*ktniA}};
    vector<T> coeffs_6={{0, kbAa*ktanA*ktia0*ktiaA*ktni0 + kbAa*ktanA*ktia0*ktinA*ktni0 + kbAa*ktanA*ktia0*ktni0*kuAi + kbAi*ktan0*ktanA*ktiaA*ktni0 + kbAi*ktan0*ktanA*ktinA*ktni0 + kbAi*ktan0*ktinA*ktni0*kuAa + kbAn*ktan0*ktanA*ktia0*ktiaA + kbAn*ktan0*ktanA*ktia0*ktinA + kbAn*ktan0*ktanA*ktia0*kuAi + kbAn*ktan0*ktanA*ktiaA*ktin0 + kbAn*ktan0*ktanA*ktin0*ktinA + kbAn*ktan0*ktanA*ktin0*kuAi + kbAn*ktan0*ktia0*ktiaA*kuAa + kbAn*ktan0*ktia0*ktinA*kuAa + kbAn*ktan0*ktia0*kuAa*kuAi + kbAn*ktan0*ktiaA*ktin0*kuAa + kbAn*ktan0*ktin0*ktinA*kuAa + kbAn*ktan0*ktin0*kuAa*kuAi, kbAa*kbAi*ktanA*ktiaA*ktni0 + kbAa*kbAi*ktanA*ktinA*ktni0 + kbAa*kbAn*ktanA*ktia0*ktiaA + kbAa*kbAn*ktanA*ktia0*ktinA + kbAa*kbAn*ktanA*ktia0*kuAi + kbAa*kbAn*ktanA*ktiaA*ktin0 + kbAa*kbAn*ktanA*ktin0*ktinA + kbAa*kbAn*ktanA*ktin0*kuAi + kbAi*kbAn*ktan0*ktanA*ktiaA + kbAi*kbAn*ktan0*ktanA*ktinA + kbAi*kbAn*ktan0*ktiaA*kuAa + kbAi*kbAn*ktan0*ktinA*kuAa, kbAa*kbAi*kbAn*ktanA*ktiaA + kbAa*kbAi*kbAn*ktanA*ktinA}};
    T numdeg0=ktan0*coeffs_1[0];
    T numdeg1=ktanA*coeffs_4[1]+ktan0*coeffs_1[1];
    T numdeg2=ktanA*coeffs_4[2]+ktan0*coeffs_1[2];
    T numdeg3=ktanA*coeffs_4[3];
    T dendeg0=coeffs_1[0]+coeffs_2[0]+coeffs_3[0];
    T dendeg1=coeffs_1[1]+coeffs_2[1]+coeffs_3[1]+coeffs_4[1]+coeffs_5[1]+coeffs_6[1];
    T dendeg2=coeffs_1[2]+coeffs_2[2]+coeffs_3[2]+coeffs_4[2]+coeffs_5[2]+coeffs_6[2];
    T dendeg3=coeffs_4[3]+coeffs_5[3]+coeffs_6[3];
    num={numdeg0,numdeg1,numdeg2,numdeg3};
    den={dendeg0,dendeg1,dendeg2,dendeg3};
}

template <typename T>
void rhos_GRF_PolA_A(py::array_t<double> parsar, vector<T> &rhos, py::array_t<double>othervars, double valGRF){

    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    T ktia0=pars[0];
    T ktan0=pars[1];
    T ktin0=pars[2];
    T ktni0=pars[3];
    T ktiaA=pars[4];
    T ktanA=pars[5];
    T ktinA=pars[6];
    T ktniA=pars[7];
    T kbAa=pars[8];
    T kuAa=pars[9];
    T kbAi=pars[10];
    T kuAi=pars[11];
    T kbAn=pars[12];
    T kuAn=pars[13];
    T rho_1=(ktanA*ktia0*ktiaA*ktni0*kuAn + ktanA*ktia0*ktinA*ktni0*kuAn + ktanA*ktia0*ktni0*ktniA*kuAi + ktanA*ktia0*ktni0*kuAi*kuAn + ktia0*ktiaA*ktni0*ktniA*kuAa + ktia0*ktiaA*ktni0*kuAa*kuAn + ktia0*ktinA*ktni0*kuAa*kuAn + ktia0*ktni0*ktniA*kuAa*kuAi + ktia0*ktni0*kuAa*kuAi*kuAn)*pow(valGRF,0)+(kbAi*ktiaA*ktni0*ktniA*kuAa + kbAi*ktiaA*ktni0*kuAa*kuAn + kbAn*ktanA*ktia0*ktniA*kuAi + kbAn*ktia0*ktiaA*ktniA*kuAa + kbAn*ktia0*ktniA*kuAa*kuAi + kbAn*ktiaA*ktin0*ktniA*kuAa)*pow(valGRF,1)+(kbAi*kbAn*ktiaA*ktniA*kuAa)*pow(valGRF,2);
    T rho_2=(ktan0*ktanA*ktiaA*ktni0*kuAn + ktan0*ktanA*ktinA*ktni0*kuAn + ktan0*ktanA*ktni0*ktniA*kuAi + ktan0*ktanA*ktni0*kuAi*kuAn + ktan0*ktiaA*ktni0*ktniA*kuAa + ktan0*ktiaA*ktni0*kuAa*kuAn + ktan0*ktinA*ktni0*kuAa*kuAn + ktan0*ktni0*ktniA*kuAa*kuAi + ktan0*ktni0*kuAa*kuAi*kuAn)*pow(valGRF,0)+(kbAa*ktanA*ktiaA*ktni0*kuAn + kbAa*ktanA*ktinA*ktni0*kuAn + kbAa*ktanA*ktni0*ktniA*kuAi + kbAa*ktanA*ktni0*kuAi*kuAn + kbAn*ktan0*ktanA*ktniA*kuAi + kbAn*ktan0*ktniA*kuAa*kuAi)*pow(valGRF,1)+(kbAa*kbAn*ktanA*ktniA*kuAi)*pow(valGRF,2);
    T rho_3=(ktan0*ktanA*ktia0*ktiaA*kuAn + ktan0*ktanA*ktia0*ktinA*kuAn + ktan0*ktanA*ktia0*ktniA*kuAi + ktan0*ktanA*ktia0*kuAi*kuAn + ktan0*ktanA*ktiaA*ktin0*kuAn + ktan0*ktanA*ktin0*ktinA*kuAn + ktan0*ktanA*ktin0*ktniA*kuAi + ktan0*ktanA*ktin0*kuAi*kuAn + ktan0*ktia0*ktiaA*ktniA*kuAa + ktan0*ktia0*ktiaA*kuAa*kuAn + ktan0*ktia0*ktinA*kuAa*kuAn + ktan0*ktia0*ktniA*kuAa*kuAi + ktan0*ktia0*kuAa*kuAi*kuAn + ktan0*ktiaA*ktin0*ktniA*kuAa + ktan0*ktiaA*ktin0*kuAa*kuAn + ktan0*ktin0*ktinA*kuAa*kuAn + ktan0*ktin0*ktniA*kuAa*kuAi + ktan0*ktin0*kuAa*kuAi*kuAn)*pow(valGRF,0)+(kbAa*ktanA*ktia0*ktiaA*kuAn + kbAa*ktanA*ktia0*ktinA*kuAn + kbAa*ktanA*ktia0*kuAi*kuAn + kbAa*ktanA*ktiaA*ktin0*kuAn + kbAa*ktanA*ktin0*ktinA*kuAn + kbAa*ktanA*ktin0*ktniA*kuAi + kbAa*ktanA*ktin0*kuAi*kuAn + kbAi*ktan0*ktanA*ktiaA*kuAn + kbAi*ktan0*ktanA*ktinA*kuAn + kbAi*ktan0*ktiaA*ktniA*kuAa + kbAi*ktan0*ktiaA*kuAa*kuAn + kbAi*ktan0*ktinA*kuAa*kuAn)*pow(valGRF,1)+(kbAa*kbAi*ktanA*ktiaA*kuAn + kbAa*kbAi*ktanA*ktinA*kuAn)*pow(valGRF,2);
    T rho_4=(kbAa*ktia0*ktiaA*ktni0*ktniA + kbAa*ktia0*ktiaA*ktni0*kuAn + kbAa*ktia0*ktinA*ktni0*kuAn + kbAa*ktia0*ktni0*ktniA*kuAi + kbAa*ktia0*ktni0*kuAi*kuAn + kbAi*ktan0*ktiaA*ktni0*ktniA + kbAi*ktan0*ktiaA*ktni0*kuAn + kbAn*ktan0*ktia0*ktiaA*ktniA + kbAn*ktan0*ktiaA*ktin0*ktniA)*pow(valGRF,1)+(kbAa*kbAi*ktiaA*ktni0*ktniA + kbAa*kbAi*ktiaA*ktni0*kuAn + kbAa*kbAn*ktia0*ktiaA*ktniA + kbAa*kbAn*ktia0*ktniA*kuAi + kbAa*kbAn*ktiaA*ktin0*ktniA + kbAi*kbAn*ktan0*ktiaA*ktniA)*pow(valGRF,2)+(kbAa*kbAi*kbAn*ktiaA*ktniA)*pow(valGRF,3);
    T rho_5=(kbAa*ktanA*ktia0*ktni0*ktniA + kbAi*ktan0*ktanA*ktni0*ktniA + kbAi*ktan0*ktanA*ktni0*kuAn + kbAi*ktan0*ktni0*ktniA*kuAa + kbAi*ktan0*ktni0*kuAa*kuAn + kbAn*ktan0*ktanA*ktia0*ktniA + kbAn*ktan0*ktanA*ktin0*ktniA + kbAn*ktan0*ktia0*ktniA*kuAa + kbAn*ktan0*ktin0*ktniA*kuAa)*pow(valGRF,1)+(kbAa*kbAi*ktanA*ktni0*ktniA + kbAa*kbAi*ktanA*ktni0*kuAn + kbAa*kbAn*ktanA*ktia0*ktniA + kbAa*kbAn*ktanA*ktin0*ktniA + kbAi*kbAn*ktan0*ktanA*ktniA + kbAi*kbAn*ktan0*ktniA*kuAa)*pow(valGRF,2)+(kbAa*kbAi*kbAn*ktanA*ktniA)*pow(valGRF,3);
    T rho_6=(kbAa*ktanA*ktia0*ktiaA*ktni0 + kbAa*ktanA*ktia0*ktinA*ktni0 + kbAa*ktanA*ktia0*ktni0*kuAi + kbAi*ktan0*ktanA*ktiaA*ktni0 + kbAi*ktan0*ktanA*ktinA*ktni0 + kbAi*ktan0*ktinA*ktni0*kuAa + kbAn*ktan0*ktanA*ktia0*ktiaA + kbAn*ktan0*ktanA*ktia0*ktinA + kbAn*ktan0*ktanA*ktia0*kuAi + kbAn*ktan0*ktanA*ktiaA*ktin0 + kbAn*ktan0*ktanA*ktin0*ktinA + kbAn*ktan0*ktanA*ktin0*kuAi + kbAn*ktan0*ktia0*ktiaA*kuAa + kbAn*ktan0*ktia0*ktinA*kuAa + kbAn*ktan0*ktia0*kuAa*kuAi + kbAn*ktan0*ktiaA*ktin0*kuAa + kbAn*ktan0*ktin0*ktinA*kuAa + kbAn*ktan0*ktin0*kuAa*kuAi)*pow(valGRF,1)+(kbAa*kbAi*ktanA*ktiaA*ktni0 + kbAa*kbAi*ktanA*ktinA*ktni0 + kbAa*kbAn*ktanA*ktia0*ktiaA + kbAa*kbAn*ktanA*ktia0*ktinA + kbAa*kbAn*ktanA*ktia0*kuAi + kbAa*kbAn*ktanA*ktiaA*ktin0 + kbAa*kbAn*ktanA*ktin0*ktinA + kbAa*kbAn*ktanA*ktin0*kuAi + kbAi*kbAn*ktan0*ktanA*ktiaA + kbAi*kbAn*ktan0*ktanA*ktinA + kbAi*kbAn*ktan0*ktiaA*kuAa + kbAi*kbAn*ktan0*ktinA*kuAa)*pow(valGRF,2)+(kbAa*kbAi*kbAn*ktanA*ktiaA + kbAa*kbAi*kbAn*ktanA*ktinA)*pow(valGRF,3);
    rhos={rho_1,rho_2,rho_3,rho_4,rho_5,rho_6};
}

template <typename T, typename Tpolyc, typename polytype, typename thresholdtype>
class GRF_PolA_GRFCalculations: public GRFCalculations <T, Tpolyc, polytype, thresholdtype>{

    public:

    void fill_num_den(py::array_t<double> parsar, py::array_t<double>othervars){
            
        GRF_PolA_A<T>(parsar, this->num,this->den,othervars);
    }

    void fill_rhos(py::array_t<double> parsar, py::array_t<double>othervars, double xval){

        rhos_GRF_PolA_A<T>(parsar, this->rhos,othervars, xval);
    }
    
};
typedef number<mpc_complex_backend<50> > mpc_50;
typedef Polynomial<50> polytype_50;
typedef number<mpfr_float_backend<15> > mpfr_15;
typedef number<mpfr_float_backend<50> > mpfr_50;
typedef number<mpfr_float_backend<100> > mpfr_100;
typedef number<mpc_complex_backend<100> > mpc_100;
typedef Polynomial<100> polytype_100;


PYBIND11_MODULE(PolA_A_allpars, m) {
    py::class_<GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>>(m, "GRFCalculations_ld_50_15", py::module_local())
    .def(py::init())
    .def("fill_num_den", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::fill_num_den)
    .def("fill_rhos", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::fill_rhos)
    .def("interfaceGRF", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::interfaceGRF)
    .def("getrhos", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::getrhos)
    .def("interfaceps", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::interfaceps,py::arg("verbose")=false, py::arg("thresholdimag")=1e-15,py::arg("minx0")=false,py::arg("maxx1")=false,py::arg("writex05coeffs")=false, py::arg("absder")=false, py::arg("normalisefirst")=true, py::arg("fnamecoeffs")="coefficients.txt")
    .def("interface_return_num_den", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::interface_return_num_den)
    .def("interfacemonotonic", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::interfacemonotonic, py::arg("thresholdimag")=1e-15)
    .def("print_num_den", &GRF_PolA_GRFCalculations<long double,mpc_50,polytype_50,mpfr_15>::print_num_den);

    py::class_<GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>>(m, "GRFCalculations_50_50_15", py::module_local())
    .def(py::init())
    .def("fill_num_den", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::fill_num_den)
    .def("fill_rhos", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::fill_rhos)
    .def("interfaceGRF", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::interfaceGRF)
    .def("getrhos", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::getrhos)
    .def("interfaceps", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::interfaceps,py::arg("verbose")=false, py::arg("thresholdimag")=1e-15,py::arg("minx0")=false,py::arg("maxx1")=false,py::arg("writex05coeffs")=false, py::arg("absder")=false, py::arg("normalisefirst")=true, py::arg("fnamecoeffs")="coefficients.txt")
    .def("interface_return_num_den", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::interface_return_num_den)
    .def("interfacemonotonic", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::interfacemonotonic, py::arg("thresholdimag")=1e-15)
    .def("print_num_den", &GRF_PolA_GRFCalculations<mpfr_50,mpc_50,polytype_50,mpfr_15>::print_num_den);

    py::class_<GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>>(m, "GRFCalculations_100_100_15", py::module_local())
    .def(py::init())
    .def("fill_num_den", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::fill_num_den)
    .def("fill_rhos", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::fill_rhos)
    .def("interfaceGRF", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::interfaceGRF)
    .def("getrhos", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::getrhos)
    .def("interfaceps", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::interfaceps,py::arg("verbose")=false, py::arg("thresholdimag")=1e-15,py::arg("minx0")=false,py::arg("maxx1")=false,py::arg("writex05coeffs")=false, py::arg("absder")=false, py::arg("normalisefirst")=true, py::arg("fnamecoeffs")="coefficients.txt")
    .def("interface_return_num_den", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::interface_return_num_den)
    .def("interfacemonotonic", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::interfacemonotonic, py::arg("thresholdimag")=1e-15)
    .def("print_num_den", &GRF_PolA_GRFCalculations<mpfr_100,mpc_100,polytype_100,mpfr_15>::print_num_den);

}
