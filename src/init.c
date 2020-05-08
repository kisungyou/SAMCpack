/* RUN AT THE STAGE OF IMPORTING RSAGEO 
 * FROM : tools::package_native_routine_registration_skeleton(".")
 *        tools::package_native_routine_registration_skeleton(".", character_only = FALSE) for repeated use
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SAMCpack_cpp_ttest_statistic(SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_SAMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_samcfast_sexpdata(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_samcfast_type0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_samcfast_type1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_samcfast_type2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_samcfast_type3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_exec_SAMCoptim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_find_mingeq(SEXP, SEXP);
extern SEXP _SAMCpack_ksSAMCrsa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_rescale_hori2(SEXP);
extern SEXP _SAMCpack_rescale_vert2(SEXP);
extern SEXP _SAMCpack_samc_cpp_ttest(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SAMCpack_sample_int(SEXP, SEXP);
extern SEXP _SAMCpack_two_perm_vec(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SAMCpack_cpp_ttest_statistic",    (DL_FUNC) &_SAMCpack_cpp_ttest_statistic,     3},
    {"_SAMCpack_exec_SAMC",              (DL_FUNC) &_SAMCpack_exec_SAMC,              12},
    {"_SAMCpack_exec_samcfast_sexpdata", (DL_FUNC) &_SAMCpack_exec_samcfast_sexpdata, 13},
    {"_SAMCpack_exec_samcfast_type0",    (DL_FUNC) &_SAMCpack_exec_samcfast_type0,    12},
    {"_SAMCpack_exec_samcfast_type1",    (DL_FUNC) &_SAMCpack_exec_samcfast_type1,    13},
    {"_SAMCpack_exec_samcfast_type2",    (DL_FUNC) &_SAMCpack_exec_samcfast_type2,    13},
    {"_SAMCpack_exec_samcfast_type3",    (DL_FUNC) &_SAMCpack_exec_samcfast_type3,    13},
    {"_SAMCpack_exec_SAMCoptim",         (DL_FUNC) &_SAMCpack_exec_SAMCoptim,         15},
    {"_SAMCpack_find_mingeq",            (DL_FUNC) &_SAMCpack_find_mingeq,             2},
    {"_SAMCpack_ksSAMCrsa",              (DL_FUNC) &_SAMCpack_ksSAMCrsa,               9},
    {"_SAMCpack_rescale_hori2",          (DL_FUNC) &_SAMCpack_rescale_hori2,           1},
    {"_SAMCpack_rescale_vert2",          (DL_FUNC) &_SAMCpack_rescale_vert2,           1},
    {"_SAMCpack_samc_cpp_ttest",         (DL_FUNC) &_SAMCpack_samc_cpp_ttest,         11},
    {"_SAMCpack_sample_int",             (DL_FUNC) &_SAMCpack_sample_int,              2},
    {"_SAMCpack_two_perm_vec",           (DL_FUNC) &_SAMCpack_two_perm_vec,            2},
    {NULL, NULL, 0}
};

void R_init_SAMCpack(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
