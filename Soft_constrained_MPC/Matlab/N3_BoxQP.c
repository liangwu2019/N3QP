#include <math.h>
#include <string.h>
#include "mex.h"
#include "lapack.h"
#include "blas.h"

const ptrdiff_t ione = 1;
const double fone = 1.0;
const double mone = -1.0;
const double fzero = 0.0;

void N3_BoxQP(double *H, double *h, double epsilon, double alpha, double delta, ptrdiff_t n, double *z);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	ptrdiff_t n, i;
	double *z;
	
	double *H = mxGetPr(prhs[0]);
	double *h = mxGetPr(prhs[1]);
	double epsilon = mxGetScalar(prhs[2]);
    double alpha = mxGetScalar(prhs[3]);
    double delta = mxGetScalar(prhs[4]);
	n = mxGetM(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	z = mxGetPr(plhs[0]);
	memset(z,0,sizeof(double)*n);
	N3_BoxQP(H,h,epsilon,alpha,delta,n,z);
	return;
}

void N3_BoxQP(double *H, double *h, double epsilon, double alpha, double delta, ptrdiff_t n, double *z)
{
	ptrdiff_t i, j, k, iter, info;
	double h_max_norm = 0.0;
	
	double ub = 1.0+delta;
	double lb = 1.0/(ub);
	double sigma = sqrt(2)*delta*(1+delta)*(1+delta)*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)*(1+delta)*alpha*alpha/(1-alpha);        
	double beta = (alpha-sigma)/(1.0+alpha/sqrt(2*n));
	double decrease_rate = 1.0-beta/sqrt(2*n);
    double tau = 1.0;
	double lambda = alpha/sqrt(2*n);
    ptrdiff_t max_iter = ceil(-log(2*n*(1+alpha)/epsilon)/log(decrease_rate));
	double *gamma, *theta, *phi, *psi, *temp, *temp1, *temp2, *temp_r1, *temp_r2, *delta_z_temp, *delta_z;
	double *gamma_hat, *theta_hat, *phi_hat, *psi_hat;
	double *M;
	int *Idx;
	double ratio_gamma, ratio_theta, ratio_phi, ratio_psi;
	int flag_gamma, flag_theta, flag_phi, flag_psi;
	int count;
	double *d;
	double coeff;
	gamma = malloc(sizeof(double)*n);
	theta = malloc(sizeof(double)*n);
	phi = malloc(sizeof(double)*n);
	psi = malloc(sizeof(double)*n);
    temp = malloc(sizeof(double)*n);
	temp1 = malloc(sizeof(double)*n);
	temp2 = malloc(sizeof(double)*n);
	temp_r1 = malloc(sizeof(double)*n);
	temp_r2 = malloc(sizeof(double)*n);
	delta_z_temp = malloc(sizeof(double)*n);
	delta_z = malloc(sizeof(double)*n);
	gamma_hat = malloc(sizeof(double)*n);
	theta_hat = malloc(sizeof(double)*n);
	phi_hat = malloc(sizeof(double)*n);
	psi_hat = malloc(sizeof(double)*n);
	M = malloc(sizeof(double)*n*n);
	Idx = malloc(sizeof(int)*n);
	d = malloc(sizeof(double)*n);

	for(i=0;i<n;i++)
        h_max_norm = h_max_norm > fabs(h[i]) ? h_max_norm : fabs(h[i]);
    if(h_max_norm==0.0)
    	return;
    else
    {
    	for(i=0;i<n;i++)
        {
            gamma[i] = 1.0 - lambda/h_max_norm * h[i];
            theta[i] = 1.0 + lambda/h_max_norm * h[i];
            phi[i] = 1.0;
            psi[i] = 1.0;
        }
        memcpy(gamma_hat,gamma,sizeof(double)*n);
        memcpy(theta_hat,theta,sizeof(double)*n);
        memcpy(phi_hat,phi,sizeof(double)*n);
        memcpy(psi_hat,psi,sizeof(double)*n);
    	for(i=0;i<n;i++)
    		for(j=0;j<n;j++)
    			M[i*n+j] = 2.0*lambda/h_max_norm*H[i*n+j];
    	for(i=0;i<n;i++)
        {
    		M[i*n+i] += 2.0;
            d[i] = 2.0;
        }
    	dpotrf("l",&n,M,&n,&info);
    	dpotri("l",&n,M,&n,&info);
        memset(Idx,0,sizeof(int)*n);
        for(iter=0;iter<max_iter;iter++)
        {     
        	count = -1;
        	for(i=0;i<n;i++)
        	{
        		ratio_gamma = gamma_hat[i]/gamma[i]; 
        		ratio_theta = theta_hat[i]/theta[i];
        		ratio_phi = phi_hat[i]/phi[i];
        		ratio_psi = psi_hat[i]/psi[i];
        		flag_gamma = (ratio_gamma<lb || ratio_gamma>ub);
        		flag_theta = (ratio_theta<lb || ratio_theta>ub);
        		flag_phi = (ratio_phi<lb || ratio_phi>ub);
        		flag_psi = (ratio_psi<lb || ratio_psi>ub);
        		if(flag_gamma)
        			gamma_hat[i] = gamma[i];
        		if(flag_theta)
        			theta_hat[i] = theta[i];
        		if(flag_phi)
        			phi_hat[i] = phi[i];
        		if(flag_psi)
        			psi_hat[i] = psi[i];
        		if(flag_gamma || flag_theta || flag_phi || flag_psi)
        		{
        			count += 1;
        			Idx[count] = i;
        		}          		  	
        	}
            for(i=0;i<n;i++)
            {
                temp1[i] = gamma_hat[i]/phi_hat[i];
                temp2[i] = theta_hat[i]/psi_hat[i];
                temp_r1[i] = (tau-gamma[i]*phi[i])/phi_hat[i];
                temp_r2[i] = (tau-theta[i]*psi[i])/psi_hat[i];             
                delta_z_temp[i] = temp_r2[i] - temp_r1[i];  
            }

        	for(j=0;j<=count;j++)
        	{
        		i = Idx[j];
        		coeff = -(temp1[i]+temp2[i]-d[i])/(1+(temp1[i]+temp2[i]-d[i])*M[i*n+i]);
                 for(k=0;k<n;k++) 
                    temp[k] = (k < i) ? M[k*n + i] : M[k + i*n];
        		dsyr("l",&n,&coeff,temp,&ione,M,&n);
                d[i] = temp1[i]+temp2[i];
        	}
        	// calculate delta_z
        	dsymv("l",&n,&fone,M,&n,delta_z_temp,&ione,&fzero,delta_z,&ione);

        	daxpy(&n,&fone,delta_z,&ione,z,&ione);    
            for(i=0;i<n;i++)
            {
                gamma[i] += (temp_r1[i] + temp1[i]*delta_z[i]);
                theta[i] += (temp_r2[i] - temp2[i]*delta_z[i]);
            }
            daxpy(&n,&mone,delta_z,&ione,phi,&ione);
			daxpy(&n,&fone,delta_z,&ione,psi,&ione);
            tau = decrease_rate*tau;
        }
    }
    return;
}