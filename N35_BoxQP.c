#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>

const int ione = 1;
const double fone = 1.0;
const double mone = -1.0;
const double fzero = 0.0;

void N35_BoxQP(double *H, double *h, double epsilon, double alpha, int n, double *z)
{
    int i, j, iter, info;
    double h_max_norm = 0.0;
    double sigma = 0.5*alpha*alpha/(1-alpha);
    double beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
    double decrease_rate = 1.0-beta/sqrt(2*n);
    double tau = 1.0;
    double lambda = alpha/sqrt(2*n);
    int max_iter = ceil(-log(2*n*(1+alpha)/epsilon)/log(decrease_rate));
    double *gamma, *theta, *phi, *psi, *temp1, *temp2, *temp_r1, *temp_r2, *delta_z, *H_hat_old, *H_hat;
    gamma = malloc(sizeof(double)*n);
    theta = malloc(sizeof(double)*n);
    phi = malloc(sizeof(double)*n);
    psi = malloc(sizeof(double)*n);
    temp1 = malloc(sizeof(double)*n);
    temp2 = malloc(sizeof(double)*n);
    temp_r1 = malloc(sizeof(double)*n);
    temp_r2 = malloc(sizeof(double)*n);
    delta_z = malloc(sizeof(double)*n);
    H_hat_old = malloc(sizeof(double)*n*n);
    H_hat = malloc(sizeof(double)*n*n);
    memset(z,0,sizeof(double)*n);
    for(i=0;i<n;i++)
        h_max_norm = h_max_norm > fabs(h[i]) ? h_max_norm : fabs(h[i]);
    if(h_max_norm<=epsilon)
        return;
    else
    {
        // initialize H_hat_old
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                H_hat_old[i*n+j] = 2.0*lambda/h_max_norm*H[i*n+j];
        // initialize gamma, theta, phi, psi
        for(i=0;i<n;i++)
        {
            gamma[i] = 1.0 - lambda/h_max_norm * h[i];
            theta[i] = 1.0 + lambda/h_max_norm * h[i];
            phi[i] = 1.0;
            psi[i] = 1.0;
        }
        // start the iterations
        for(iter=0;iter<max_iter;iter++)
        {
            memcpy(H_hat,H_hat_old,sizeof(double)*n*n);
            for(i=0;i<n;i++)
            {
                temp1[i] = gamma[i]/phi[i];
                temp2[i] = theta[i]/psi[i];
                temp_r1[i] = (tau-gamma[i]*phi[i])/phi[i];
                temp_r2[i] = (tau-theta[i]*psi[i])/psi[i];
                H_hat[i*n+i] += temp1[i] + temp2[i];
                delta_z[i] = temp_r2[i] - temp_r1[i];
            }
            info = LAPACKE_dposv(LAPACK_COL_MAJOR,'L',n,ione,H_hat,n,delta_z,n);
            cblas_daxpy(n,fone,delta_z,ione,z,ione);
            for(i=0;i<n;i++)
            {
                gamma[i] += (temp_r1[i]+temp1[i]*delta_z[i]);
                theta[i] += (temp_r2[i]-temp2[i]*delta_z[i]);                
            }
            cblas_daxpy(n,mone,delta_z,ione,phi,ione);
            cblas_daxpy(n,fone,delta_z,ione,psi,ione);
            tau = decrease_rate*tau;
        }
    }
    return;
}