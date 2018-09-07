#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>

//#include "nrutil.h"//.c"
#define pi 3.14159265
#include "lib.h"//.c"
static double  rho=0.01, eta=1.0, Rho=100, Eta=0.55;

void RSAc(double *pData,int* pDataCol,int* pDataNum, int* pSampleNum,int* pStepscale, int* pTotal_Iteration, int* pWarm,
	double *pbeta, double* pPhi,double* pSigmasq, double*  pTausq)
{

double **ppData;

//FILE *ins;
int stime;
long ltime;
int i, j, *sam, iter, truncation, truncount;
double logphi, logsigmasq, logtausq;//logkappa, 
double phi, kappa, sigmasq, tausq;
double Dphi, Dsigmasq, Dtausq;//Dkappa,
double *avebeta, avephi, avesigmasq, avetausq;
double *Dbeta, *beta;
double **Dist,**R,**V,**IV,**DRphi, **DRkappa, *mu,**C,*zmu,*b;
double delta, sum, a; //max;
double low[6],up[6],sze[6],lower[6],upper[6], add[*pDataCol+1], delctrl;

double rnormval;
GetRNGstate();
rnormval = norm_rand();
PutRNGstate();



int sampleNum = *pSampleNum;
int dataCol = *pDataCol;
int dataNum = *pDataNum;
int stepscale=*pStepscale;
int total_iteration=*pTotal_Iteration;
int warm = *pWarm;

//Rprintf("initialize\n");

ppData=dmatrix(1,*pDataNum,1,*pDataCol);
for(i = 1;i<=*pDataNum;i++)
for(j = 1;j<=*pDataCol;j++)//{
	ppData[i][j] = pData[(j-1)*(*pDataNum)+(i-1)];
//Rprintf("ppData %d %d : %lf\n",i,j,ppData[i][j]);}


// initialize the random number generator 
/*
 * ltime=time(NULL);
 stime=(unsigned int)ltime/2;
 srand(stime);
 */
  
  sam=ivector(1,sampleNum); 
  R=dmatrix(1,sampleNum,1,sampleNum);
  V=dmatrix(1,sampleNum,1,sampleNum);
  IV=dmatrix(1,sampleNum,1,sampleNum);
  Dist=dmatrix(1,sampleNum,1,sampleNum);
  DRphi=dmatrix(1,sampleNum,1,sampleNum);
  DRkappa=dmatrix(1,sampleNum,1,sampleNum);
  mu=dvector(1,sampleNum);
//  C=dvector(1,sampleNum);
  C=dmatrix(1,sampleNum,1,dataCol);
  zmu=dvector(1,sampleNum);
  b=dvector(1,sampleNum);

  Dbeta = dvector(1,dataCol-2);
  beta = dvector(1,dataCol-2);
  avebeta = dvector(1,dataCol-2);
 // Initialization of the parameters 

	low[1]=-2.0; up[1]=2.0; low[2]=-2.0; up[2]=2.0;
  low[3]=0.0;  up[3]=4.0; low[4]=-2.0; up[4]=2.0; low[5]=-2.0; up[5]=2.0;
  sze[1]=1.0; sze[2]=1.0; sze[3]=1.0; sze[4]=1.0; sze[5]=1.0;
  
  avephi=avesigmasq=avetausq=0.0;
//  for(j=1;j<=dataCol-1;j++) avebeta[j]=0;
  for(j=1;j<=dataCol-2;j++) avebeta[j]=0;
  truncount=0; truncation=1;

  // Algorithm Iteration  
  for(iter=1; iter<=total_iteration; iter++){
ABC:
    if(truncation==1){ 
      
      //        beta0=rand()*1.0/RAND_MAX*(up[1]-low[1])+low[1];
      //        beta1=rand()*1.0/RAND_MAX*(up[2]-low[2])+low[2];
      for(j=1;j<=dataCol-2;j++){
        GetRNGstate();
        rnormval = norm_rand();
        PutRNGstate();
        beta[j]=rnormval*1.0/RAND_MAX*(up[1]-low[1])+low[1];
      }
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logphi=rnormval*1.0/RAND_MAX*(up[3]-low[3])+low[3];
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logsigmasq=rnormval*1.0/RAND_MAX*(up[4]-low[4])+low[4]; 
      GetRNGstate();
      rnormval = norm_rand();
      PutRNGstate();
      logtausq=rnormval*1.0/RAND_MAX*(up[5]-low[5])+low[5];
      lower[1]=low[1]-sze[1]*truncount; upper[1]=up[1]+sze[1]*truncount;
      lower[2]=low[2]-sze[2]*truncount; upper[2]=up[2]+sze[2]*truncount;
      lower[3]=low[3]-sze[3]*truncount; upper[3]=up[3]+sze[3]*truncount;
      lower[4]=low[4]-sze[4]*truncount; upper[4]=up[4]+sze[4]*truncount;
      lower[5]=low[5]-sze[5]*truncount; upper[5]=up[5]+sze[5]*truncount;
    }
    

    if(iter<=stepscale){ delta=rho; delctrl=Rho; }
          else{
            delta=rho*exp(eta*log(1.0*stepscale/iter));
            delctrl=Rho*exp(Eta*log(1.0*stepscale/iter));
           }

    subset_sample(sam, sampleNum, dataNum);
 //  Rprintf("subset_sample\n");
    phi=exp(logphi); kappa=1.0; sigmasq=exp(logsigmasq); tausq=exp(logtausq);

    // V=sigmasq*R+ tausq*I 
    for(i=1; i<=sampleNum; i++){
     R[i][i]=1.0; Dist[i][i]=0.0; V[i][i]=sigmasq+tausq; DRphi[i][i]=0.0; DRkappa[i][i]=0.0;
     for(j=1; j<i; j++){
         Dist[i][j]=Dist[j][i]=sqrt(pow(ppData[sam[i]][1]-ppData[sam[j]][1],2)+pow(ppData[sam[i]][2]-ppData[sam[j]][2],2));
         a=pow(Dist[i][j]/phi,kappa);
         R[i][j]=R[j][i]=exp(-a); 
         V[i][j]=V[j][i]=R[i][j]*sigmasq; 
         DRphi[i][j]=DRphi[j][i]=a*R[i][j]*kappa/phi;
         DRkappa[i][j]=DRkappa[j][i]=-a*R[i][j]*log(Dist[i][j]/phi); 
        }
     }
    matrix_inverse(V,IV,sampleNum);
    
    for(i=1; i<=sampleNum; i++){ for(j=1;j<=dataCol-2;j++) C[i][j]=ppData[sam[i]][j+2];
                                 C[i][1]=1;
                                 mu[i]=0;
                                 for(j=1;j<=dataCol-2;j++) mu[i]+=beta[j]*C[i][j];
                                 zmu[i]=ppData[sam[i]][3]-mu[i];
                               }

    for(i=1; i<=sampleNum; i++){
        b[i]=0.0;
        for(j=1; j<=sampleNum; j++) b[i]+=IV[i][j]*zmu[j];
       }
      

    // calculate derivatives    
    for(j=1;j<=dataCol-2;j++) for(Dbeta[j]=0.0, i=1; i<=sampleNum; i++) Dbeta[j]+=b[i]*C[i][j];
//    Rprintf("Dbeta[1]=%lf\n",Dbeta[1]);
//    for(Dbeta1=0.0, i=1; i<=sampleNum; i++) Dbeta1+=C[i]*b[i];
   
    Dphi=0.0;
    for(i=1; i<=sampleNum; i++) 
       for(j=1; j<=sampleNum; j++) Dphi+=b[i]*DRphi[i][j]*b[j];
    for(sum=0.0, i=1; i<=sampleNum; i++)
        for(j=1; j<=sampleNum; j++) sum+=IV[i][j]*DRphi[j][i];
    Dphi-=sum;
    Dphi*=0.5*sigmasq;
    
    Dsigmasq=0.0;
    for(i=1; i<=sampleNum; i++)
       for(j=1; j<=sampleNum; j++) Dsigmasq+=b[i]*R[i][j]*b[j];
    for(sum=0.0,i=1; i<=sampleNum; i++)
        for(j=1; j<=sampleNum; j++) sum+=IV[i][j]*R[j][i];
    Dsigmasq-=sum;
    Dsigmasq*=0.5; 
     
    Dtausq=0.0;
    for(i=1; i<=sampleNum; i++) Dtausq+=b[i]*b[i];
    for(sum=0.0,i=1; i<=sampleNum; i++) sum+=IV[i][i];
    Dtausq-=sum;
    Dtausq*=0.5;


    for(j=1;j<=dataCol-2;j++) add[j]=delta*Dbeta[j]; 
    add[dataCol-1]=delta*Dphi*phi; add[dataCol]=delta*Dsigmasq*sigmasq; add[dataCol+1]=delta*Dtausq*tausq;

    //beta0+=add[1];
    //beta1+=add[2];
    for(j=1;j<=dataCol-2;j++) beta[j]+=add[j];
//    Rprintf("add[1]=%lf,beta[1]=%lf\n",add[1],beta[1]);
    logphi+=add[dataCol-1];
    logsigmasq+=add[dataCol]; 
    logtausq+=add[dataCol+1];  

    truncation=0;
    //for(sum=0.0, i=1; i<=5; i++) sum+=add[i]*add[i];
    for(sum=0.0, j=1; j<=dataCol+1; j++) sum+=add[j]*add[j];
    if(sqrt(sum)>delctrl) truncation=1;
       else{
         for(j=1;j<=dataCol-2;j++) if(beta[j]<lower[1] || beta[j]>upper[1]) truncation=1;
 //          else if(beta1<lower[2] || beta1>upper[2]) truncation=1;
              else if(logphi<lower[3] || logphi>upper[3]) truncation=1;
                else if(logsigmasq<lower[4] || logsigmasq>upper[4]) truncation=1;
                 else if(logtausq<lower[5] || logtausq>upper[5]) truncation=1;
        }
    if(truncation==1){ truncount++; 
                      goto ABC; 
                     }

    if(iter>warm){ avephi+=logphi/(total_iteration-warm); avesigmasq+=logsigmasq/(total_iteration-warm); avetausq+=logtausq/(total_iteration-warm); 
			for(j=1;j<=dataCol-2;j++) avebeta[j]+=beta[j]/(total_iteration-warm);
		}
//      else if(iter==warm){
//              ins=fopen("bb100.sum","a");
//              fprintf(ins, " %g %g %g %g %g\n",  beta0,beta1,exp(logphi),exp(logsigmasq),exp(logtausq));
//              fclose(ins);
//           }
    }  // end iteration 
  
//    avebeta0/=(total_iteration-warm); avebeta1/=(total_iteration-warm);
//    avephi/=(total_iteration-warm); avesigmasq/=(total_iteration-warm);
//    avetausq/=(total_iteration-warm);
    
//    ins=fopen("aa100.sum","a");
//    fprintf(ins, "%g %g %g %g %g %g %g %g %g %g\n",  beta0,beta1,exp(logphi),exp(logsigmasq),exp(logtausq), 
//            avebeta0,avebeta1,exp(avephi),exp(avesigmasq),exp(avetausq));
//    Rprintf("%g %g %g %g %g %g %g %g %g %g\n",  beta0,beta1,exp(logphi),exp(logsigmasq),exp(logtausq), 
//           avebeta0,avebeta1,exp(avephi),exp(avesigmasq),exp(avetausq));
//    fclose(ins);

    for(j=1;j<=dataCol-2;j++) pbeta[j-1]=beta[j];
//    *pbeta0=beta0; 
//    *pbeta1=beta1; 
    *pPhi=exp(logphi),
    *pSigmasq=exp(logsigmasq); 
    *pTausq=exp(logtausq);

  //  Rprintf("%g %g %g %g %g \n",  *pbeta0,*pbeta1,*pPhi,*pSigmasq,*pTausq);


  free_ivector(sam,1,sampleNum);
  free_dvector(Dbeta,1,dataCol-2);
  free_dvector(beta,1,dataCol-2);
  free_dvector(avebeta,1,dataCol);
  free_dmatrix(R,1,sampleNum,1,sampleNum);
  free_dmatrix(V,1,sampleNum,1,sampleNum);
  free_dmatrix(IV,1,sampleNum,1,sampleNum);
  free_dmatrix(Dist,1,sampleNum,1,sampleNum);
  free_dmatrix(DRphi,1,sampleNum,1,sampleNum);
  free_dmatrix(DRkappa,1,sampleNum,1,sampleNum);
  free_dvector(mu,1,sampleNum);
  // free_dvector(C,1,sampleNum);
  free_dmatrix(C,1,sampleNum,1,dataCol-2);
  free_dvector(zmu,1,sampleNum);
  free_dvector(b,1,sampleNum);
  free_dmatrix(ppData,1,*pDataNum,1,*pDataCol);
}
