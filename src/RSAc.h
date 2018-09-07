#ifndef _SAMCpack_RSAC_H
#define _SAMCpack_RSAC_H

#ifdef __cplusplus
extern "C" {
#endif
  
  void RSAc(double *pData,int* pDataCol,int* pDataNum, int* pSampleNum,int* pStepscale, int* pTotal_Iteration, int* pWarm,
            double *pbeta, double* pPhi,double* pSigmasq, double*  pTausq);
  

  
#ifdef __cplusplus
}
#endif


#endif
