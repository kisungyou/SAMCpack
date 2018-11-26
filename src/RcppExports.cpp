// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// RSAarma
Rcpp::List RSAarma(arma::vec pData, int pDataCol, int pDataNum, int pSampleNum, int pStepscale, int pTotal_Iteration, int pWarm);
RcppExport SEXP _SAMCpack_RSAarma(SEXP pDataSEXP, SEXP pDataColSEXP, SEXP pDataNumSEXP, SEXP pSampleNumSEXP, SEXP pStepscaleSEXP, SEXP pTotal_IterationSEXP, SEXP pWarmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pData(pDataSEXP);
    Rcpp::traits::input_parameter< int >::type pDataCol(pDataColSEXP);
    Rcpp::traits::input_parameter< int >::type pDataNum(pDataNumSEXP);
    Rcpp::traits::input_parameter< int >::type pSampleNum(pSampleNumSEXP);
    Rcpp::traits::input_parameter< int >::type pStepscale(pStepscaleSEXP);
    Rcpp::traits::input_parameter< int >::type pTotal_Iteration(pTotal_IterationSEXP);
    Rcpp::traits::input_parameter< int >::type pWarm(pWarmSEXP);
    rcpp_result_gen = Rcpp::wrap(RSAarma(pData, pDataCol, pDataNum, pSampleNum, pStepscale, pTotal_Iteration, pWarm));
    return rcpp_result_gen;
END_RCPP
}
// exec_SAMC
Rcpp::List exec_SAMC(Function func, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec& init);
RcppExport SEXP _SAMCpack_exec_SAMC(SEXP funcSEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type func(funcSEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_SAMC(func, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init));
    return rcpp_result_gen;
END_RCPP
}
// exec_samcfast_type0
Rcpp::List exec_samcfast_type0(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec init);
RcppExport SEXP _SAMCpack_exec_samcfast_type0(SEXP func_SEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_samcfast_type0(func_, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init));
    return rcpp_result_gen;
END_RCPP
}
// exec_samcfast_type1
Rcpp::List exec_samcfast_type1(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec init, arma::vec data);
RcppExport SEXP _SAMCpack_exec_samcfast_type1(SEXP func_SEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_samcfast_type1(func_, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init, data));
    return rcpp_result_gen;
END_RCPP
}
// exec_samcfast_type2
Rcpp::List exec_samcfast_type2(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec init, arma::mat data);
RcppExport SEXP _SAMCpack_exec_samcfast_type2(SEXP func_SEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_samcfast_type2(func_, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init, data));
    return rcpp_result_gen;
END_RCPP
}
// exec_samcfast_type3
Rcpp::List exec_samcfast_type3(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec init, Rcpp::List data);
RcppExport SEXP _SAMCpack_exec_samcfast_type3(SEXP func_SEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_samcfast_type3(func_, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init, data));
    return rcpp_result_gen;
END_RCPP
}
// exec_samcfast_sexpdata
Rcpp::List exec_samcfast_sexpdata(SEXP func_, const int nv, arma::vec& energy, arma::mat& domain, const double tau, const int niter, arma::vec& vecpi, const double t0, const double xi, arma::vec stepsize, arma::mat& trange, arma::vec init, SEXP data);
RcppExport SEXP _SAMCpack_exec_samcfast_sexpdata(SEXP func_SEXP, SEXP nvSEXP, SEXP energySEXP, SEXP domainSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP vecpiSEXP, SEXP t0SEXP, SEXP xiSEXP, SEXP stepsizeSEXP, SEXP trangeSEXP, SEXP initSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type func_(func_SEXP);
    Rcpp::traits::input_parameter< const int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type energy(energySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type domain(domainSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type vecpi(vecpiSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type trange(trangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(exec_samcfast_sexpdata(func_, nv, energy, domain, tau, niter, vecpi, t0, xi, stepsize, trange, init, data));
    return rcpp_result_gen;
END_RCPP
}
// rescale_vert2
arma::mat rescale_vert2(arma::mat A);
RcppExport SEXP _SAMCpack_rescale_vert2(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(rescale_vert2(A));
    return rcpp_result_gen;
END_RCPP
}
// rescale_hori2
arma::mat rescale_hori2(arma::mat A);
RcppExport SEXP _SAMCpack_rescale_hori2(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(rescale_hori2(A));
    return rcpp_result_gen;
END_RCPP
}
