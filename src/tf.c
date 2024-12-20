#include "pwsim.h"

/*////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 * A standard implementation to find the inverse of a complex square matrix
 * using Gauss-Jordan elimination algorithm. 
 * The base of this code was taken from codesansar.com and modified for 
 * using arbitrary complex matrix size.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////*/
double complex** inverter(double complex**a, int n)
{
	double complex ratio;
	int i,j,k;
	double complex **b=(double complex**)malloc(n*sizeof(double complex*));
	
	double complex **ret=(double complex**)malloc(n*sizeof(double complex*));
	
	for(i=0;i<n;i++)
	{
		b[i]  =(double complex*)malloc(2*n*sizeof(double complex));
		ret[i]=(double complex*)malloc(n*sizeof(double complex));
		for(j=0;j<n;j++)
		{
			b[i][j]=a[i][j];
		}
	}
	
	/* Augmenting Identity Matrix of Order n */
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
			{
				b[i][j+n] = 1.;
			}
			else
			{
				b[i][j+n] = 0.;
			}
		}
	}
	/* Applying Gauss Jordan Elimination */
	for(i=0;i<n;i++)
	{
		if(b[i][i] == 0.0)
		{
			printf("Mathematical Error!\n");
			exit(EXIT_FAILURE);
		}
		for(j=0;j<n;j++)
		{
			if(i!=j)
			{
				ratio = b[j][i]/b[i][i];
				for(k=0;k<2*n;k++)
				{
					b[j][k] = b[j][k] - ratio*b[i][k];
				}
			}
		}
	}
	
	/* Row Operation to Make Principal Diagonal to 1 */
	for(i=0;i<n;i++)
	{
		for(j=2*n-1;j>=0;j--)
		{
			b[i][j] = b[i][j]/b[i][i];
		}
	}
	for(i=0;i<n;i++)
		for(j=n;j<2*n;j++)
			ret[i][j-n]=b[i][j];
			
	for(i=0;i<n;i++)
		free(b[i]);
	
	free(b);
	
	return(ret);
}
/*/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////*/


/* Auxiliary function for P and SV transfer functions. It returns a 
 * 2x2 complex matrix.
 * */
double complex**halfspace(double cs,double cp,double ds,double dp,double rho_hs,double k,double w,double fact)
{
	int i,j;
	double complex **matrixhs = (double complex**)malloc(2*sizeof(double complex*));
	for(i=0;i<2;i++)
		matrixhs[i]  = (double complex *)malloc(2*sizeof(double complex));
		
	matrixhs[0][0] = 0.0;
	matrixhs[0][1] = 0.0;
	matrixhs[1][0] = 0.0;
	matrixhs[1][1] = 0.0;
	
	double Cs=cs; 
	double Cp=cp; 
	double dampS=ds; 
	double dampP=dp; 
	double wr=w; 
	double wc; 
	if(wr!=0.0)
	{
		wc = wr;
		double complex CCs = Cs/(1.-I*dampS);
		double complex CCp = Cp/(1.-I*dampP);
		
		if(k!=0)
		{
			double k2 = k*k;
			double complex wcs = (wc/CCs)*(wc/CCs);
			double complex wcp = (wc/CCp)*(wc/CCp);
			double complex ks = csqrt(k2-wcs);
			double complex kp = csqrt(k2-wcp);
			
			if(cimag(ks)<0) ks=-ks;
			if(cimag(kp)<0) kp=-kp;
			
			double complex G = rho_hs*CCs*CCs;
			
			double complex D = G*wcs/(k2-kp*ks);
			
			matrixhs[0][0] = kp*D;
			matrixhs[0][1] = k*(2.*G-D);
			matrixhs[1][0] = matrixhs[0][1];
			matrixhs[1][1] = ks*D;
			
		}
		else
		{
			matrixhs[0][0] = I*rho_hs*CCs*wc;
			matrixhs[1][1] = I*rho_hs*CCp*wc;
		}
		
		
	}
	else if(k!=0.0)
	{
		double alph =  (Cs/Cp)*(Cs/Cp);
		matrixhs[0][0] =  2.*k*rho_hs*Cs*Cs/(1.+alph);
		matrixhs[0][1] =  alph*matrixhs[0][0];
		matrixhs[1][0] =  matrixhs[0][1];
		matrixhs[1][1] =  matrixhs[0][0];
	}
	
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			matrixhs[i][j]  = matrixhs[i][j]*fact;
		
	return matrixhs;
	
}

/* Auxiliary function for P and SV transfer functions. It returns a 
 * 4x4 complex matrix.
 * */
double complex **stiflayer(double k,double wc,double ccs,double ccp,double h,double rho_sti,double ds,double dp)
{
	int i,j;
	double complex **Sa = (double complex**)malloc(4*sizeof(double complex*));
	
	for(i=0;i<4;i++)
	{
		Sa[i]  = (double complex *)malloc(4*sizeof(double complex));
		for(j=0;j<4;j++)
			Sa[i][j]=0.0;
	}
	
	double complex CCp=ccp/csqrt(1.0 - 2.0*I*dp);
	
	//printf("%f%+fi\n", crealf(CCp), cimagf(CCp));
	double complex CCs=ccs/csqrt(1.0 - 2.0*I*ds);
	
	double complex k2  = k*k;
	double complex kp  = csqrt( k2 - (wc/CCp)*(wc/CCp) );    
	double complex ks  = csqrt( k2 - (wc/CCs)*(wc/CCs) );
	
	
	
	if(cimag(ks)<0.0) ks=-ks;
	if(cimag(kp)<0.0) kp=-kp;
	
	double complex p = kp/k;
	double complex s = ks/k;
		
	double complex kph = kp*h;
	double complex ksh = ks*h;
	
	
	double complex ps  = p*s;
	double complex s2  = s*s;
	double complex Es, Ss;
	if (crealf(ksh)>0.0)
	{                            
		Es = cexp(-ksh);                        
		Ss = 0.5*(1.-Es*Es); 
	}
	else
	{
		Es = cexp(ksh);
		Ss = 0.5*(Es*Es-1.);
	}
	
	
	
	double complex Ep, Sp;
	
	if (crealf(kph)>0.0)
	{                            
		Ep = cexp(-kph);                        
		Sp = 0.5*(1.-Ep*Ep); 
	}
	else
	{
		Ep = cexp(kph);
		Sp = 0.5*(Ep*Ep-1.);
	}
	
	
	
	double complex Cs = 0.5*(1.0+Es*Es); 
	double complex Cp = 0.5*(1.0+Ep*Ep);
	
	
	
	double complex D = 2.0*(Ep*Es-Cp*Cs)+(ps+1.0/ps)*Sp*Ss;
	D = 0.5*(1.0 - s2)/D;	
	
	Sa[0][0] = D*(Cp*Ss/s-p*Cs*Sp);
	Sa[0][1] =  D*(Ep*Es-Cp*Cs+ps*Sp*Ss) + 0.5*(1.0+s2);
	Sa[0][2] = D*(p*Sp*Es-Ss*Ep/s);
	Sa[0][3] = D*(Cp*Es-Cs*Ep);
	
	Sa[1][0] = Sa[0][1];
	Sa[1][1] = D*(Cs*Sp/p-s*Cp*Ss);
	Sa[1][2] = -Sa[0][3];
	Sa[1][3] = D*(s*Ss*Ep-Sp*Es/p);
	
	Sa[2][0] = Sa[0][2];
	Sa[2][1] = Sa[1][2];
	Sa[2][2] = Sa[0][0];
	Sa[2][3] = -Sa[0][1];
	
	Sa[3][0] = Sa[0][3];
	Sa[3][1] = Sa[1][3];
	Sa[3][2] = -Sa[1][0];
	Sa[3][3] = Sa[1][1];
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			Sa[i][j]=(2.*k*rho_sti*CCs*CCs)*Sa[i][j];
		}
		
	}
	
	return Sa;
}







/* Auxiliary function for P and SV transfer functions. It returns a 
 * complex matrix which size depends on the local velocity model.
 * */
double complex **globalmatrix(double k, double w, double cs[],double cp[],double h[],int model_size, double rho_satf,double beta1, double beta2)
{
	int i,j,l;
	
	
	double complex **global_matrix = (double complex**)malloc((2*(model_size-1) + 2)*sizeof(double complex*));
	for(int i=0;i<2*(model_size-1) + 2;i++)
		global_matrix[i]  = (double complex *)malloc((2*(model_size-1) + 2)*sizeof(double complex));
	
	for(i=0;i<2*(model_size-1) + 2;i++)
		for(j=0;j<2*(model_size-1) + 2;j++)
			global_matrix[i][j]=0.0;
	
	//khalf is 2x2
	double complex **khalf=halfspace(cs[model_size-1],cp[model_size-1],beta1,beta2,rho_satf,k,w,1.0);
	
	
	double complex **mlayer;
	
	for(j=0;j<model_size-1;j++)
	{
		mlayer=stiflayer(k,w,cs[j],cp[j],h[j],rho_satf,beta1,beta2);
		
		if(j==0)
		{
			for(i=0;i<4;i++)
				for(l=0;l<4;l++)
					global_matrix[i][l]=mlayer[i][l];
		}
		else
		{
			for(i=2*j;i<2*j+4;i++)
				for(l=2*j;l<2*j+4;l++)
					global_matrix[i][l]=global_matrix[i][l]+mlayer[i-2*j][l-2*j];
		}
			
	}
	
	for(i=2*(model_size-1) ; i<2*(model_size-1)+2 ; i++)
		for(l=2*(model_size-1) ; l<2*(model_size-1)+2 ; l++)
			global_matrix[i][l]=global_matrix[i][l]+khalf[i-2*(model_size-1)][l-2*(model_size-1)];
	
	
	for(i=0;i<4;i++)
	{
		if(i<2)
			free(khalf[i]);
		free(mlayer[i]);
	}
	free(khalf);
	free(mlayer);
	
	return global_matrix;
}


/* Calculates the P and SV transfer functions
 * */
double complex ** SATF_P_SV(double *freq,int freq_size,double *vs_sta,double *vp_sta,double *espesor_sta,int model_size,double am_p, double am_sv,double angle,double rho_satf)
{
	int i;
	double w;
	double tp=angle;
	double betaS = am_sv;
	double betaP = am_p;
	double kx,kz,k;
	double complex **K;
	
	double xp=0.;
    double zp=0.;
    double A=1.;
    
    double complex **matrixfs;
    
    double complex vectorfinal[2];
    double complex **K_inv;
    
    
    double complex **ver_hor = (double complex**)malloc(2*sizeof(double complex*));
	
	for(i=0;i<2;i++)
		ver_hor[i]  = (double complex *)malloc(freq_size*sizeof(double complex));
		
	double complex u0,u1,upr;
	
	double complex *ula = (double complex*)malloc(freq_size*sizeof(double complex));
	double complex *hor_aux = (double complex*)malloc(freq_size*sizeof(double complex));
	
	for(i=0;i<freq_size;i++)
	{
		
		w=2.*M_PI*freq[i];
		kx=w/vp_sta[model_size-1]*sin(tp);
		kz=w/vp_sta[model_size-1]*cos(tp);
		k=kx;
		
		K=globalmatrix(k,w,vs_sta,vp_sta,espesor_sta,model_size,rho_satf,betaS,betaP);
		
		//matrixfs is 2x2
		matrixfs=halfspace(vs_sta[model_size-1],vp_sta[model_size-1],betaS,betaP,rho_satf,k,w,2.0);
	
		matrixfs[0][1]=0.0;
		matrixfs[1][0]=0.0;
		
		vectorfinal[0] =  conj(A*matrixfs[0][0]*sin(tp)*cexp(-I*(kx*xp+kz*zp)));
		vectorfinal[1] = conj(-A*matrixfs[1][1]*I*cos(tp)*cexp(-I*(kx*xp+kz*zp)));
		
		K_inv=inverter(K, 2*(model_size-1)+2);
		
		u0  = K_inv[0][2*(model_size-1)]*vectorfinal[0] + K_inv[0][2*(model_size-1)+1]*vectorfinal[1];
		u1  = K_inv[1][2*(model_size-1)]*vectorfinal[0] + K_inv[1][2*(model_size-1)+1]*vectorfinal[1];
		upr = K_inv[2*(model_size-1)][2*(model_size-1)]*vectorfinal[0] + K_inv[2*(model_size-1)][2*(model_size-1)+1]*vectorfinal[1];
		ula[i] = K_inv[2*(model_size-1)+1][2*(model_size-1)]*vectorfinal[0] + K_inv[2*(model_size-1)+1][2*(model_size-1)+1]*vectorfinal[1];
		
		hor_aux[i]=u1;
		ver_hor[1][i]=cabs( u0/upr );
		
	}
	
	double complex cum_sum0=0.0;
	
	for(i=0;i<freq_size;i++)
	{
		ver_hor[0][i]=cabs(smooth(hor_aux,freq_size,i,&cum_sum0,250)/ula[i]);
		//ver_hor[1][i]=smooth(ver_hor[1],freq_size,i,&cum_sum1,250);
	}
	free(hor_aux);
	free(ula);
	
	return ver_hor;
	
}

/* Calculates the SH transfer function
 * */
double complex * SATF_SH(double *freq,int freq_size,double *vs_sta,double *espesor_sta,int model_size,double am_sh, double rho_satf)
{
	int i,j;
	double complex *delta = (double complex *)malloc((model_size-1)*sizeof(double complex) );
	double complex *E = (double complex *)malloc(freq_size*sizeof(double complex) );
	double complex *F = (double complex *)malloc(freq_size*sizeof(double complex) );
	
	
    double complex **ka = (double complex**)malloc((model_size-1)*sizeof(double complex*));	
	
	for(i=0;i<model_size-1;i++)
	{
		//delta(k)=sqrt((rh(k)*G(k)*(1+1i*2*ds(k)))/(rh(k+1)*G(k+1)*(1+1i*2*ds(k+1))));
		delta[i]=csqrt((vs_sta[i]*vs_sta[i]*(1.+I*2.*am_sh))/(vs_sta[i+1]*vs_sta[i+1]*(1.+I*2.*am_sh)));
		
		
	
		ka[i]=(double complex *)malloc(freq_size*sizeof(double complex));		
		for(j=0;j<freq_size;j++)
		{
		
			E[j]=1.0;
			F[j]=1.0;
			ka[i][j]= (2.*M_PI*freq[j])/(vs_sta[i]*csqrt(1.+I*2.*am_sh));

		}
		
	}

	
	double complex a,b,c,d;
	
	for(i=1;i<model_size;i++)
	{
		for(j=0;j<freq_size;j++)
		{
			a=E[j]*(1.+delta[i-1])*cexp( I*ka[i-1][j]*espesor_sta[i-1] ) ;
			b=F[j]*(1.-delta[i-1])*cexp(-I*ka[i-1][j]*espesor_sta[i-1] ) ;
			c=E[j]*(1.-delta[i-1])*cexp( I*ka[i-1][j]*espesor_sta[i-1] ) ;
			d=F[j]*(1.+delta[i-1])*cexp(-I*ka[i-1][j]*espesor_sta[i-1] ) ;
			E[j]=0.5*(a+b);
			F[j]=0.5*(c+d);
			
		}
		
	}
	
	double complex *As1 = (double complex *)malloc(freq_size*sizeof(double complex) );
	
	for(j=0;j<freq_size;j++)
	{
		As1[j]= cabs(2.0/(E[j]+F[j])   )    ;
	}
	
	free(delta);
	free(E);
	free(F);
	for(i=0;i<model_size-1;i++)
		free(ka[i]);
	free(ka);
	
	return As1;
}
