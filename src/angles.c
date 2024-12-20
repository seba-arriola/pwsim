#include "pwsim.h"

/* Function to calculate the angle of incidence and the angle of coordinate's rotation
 * Receives the station coordinates and subfault coordinates, and pointers to
 * variables which stores the calculated angles.
 * */
void angles(double xst, double xso, double yst, double yso, double zst,  double zso, double *theta, double *phi, double *d_total,double *tp, double *ts, double *ie )
{

	int i,j,c;

	double Dhorizontal = sqrt( (xst-xso)*(xst-xso) +   (yst-yso)*(yst-yso));
	//double Dvertical = abs(zst-zso);
	double x =xst-xso;
	double y = yst-yso;
	//double z = zst-zso;

	double phiEW, phi_aux;
	if( y<0 )
		phiEW = 2*M_PI - acos(   x/sqrt(x*x+y*y) );

	else
		phiEW = acos(   x/sqrt(x*x+y*y) );

	if(phiEW<=M_PI/2.0)
		phi_aux=M_PI/2.0-phiEW;
	else
		phi_aux=(2.0*M_PI+M_PI/2.0)-phiEW;

	// added by JO to estimate correctly the angle in the Norcia case
	// correct the azimuth angle from the position of each station to each subfault
	// in this case, is an angle with respect to the North, in the direction X
	phi_aux=phi_aux-(M_PI/2.0);
	if(phi_aux<0)
		*phi=phi_aux+(2*M_PI);
	else
		*phi=phi_aux;
	//

	double *espesor=(double*)malloc(M*sizeof(double));

	espesor[0]=prof[0];

	for(i=1;i<M;i++)
		espesor[i]=prof[i]-prof[i-1];


	int indexs;

	if (-1.0*zso>prof[M-1])
	{
		indexs = M-1;
	}
	else
	{
		for(i=0;i<M;i++)
		{

			if( -1.0*zso<prof[i] )
			{
				indexs = i;
				break;
			}

		}
	}

	double aux_sum=0.0;

	for(i=0;i<indexs;i++)
		aux_sum = aux_sum+espesor[i];


	double d_ang = 0.0001;
	int size_ang = (int)((M_PI/2)/d_ang);
	double *ang = (double*)malloc(size_ang*sizeof(double));

	double *io = (double*)malloc((indexs+1)*sizeof(double));
	double dTOTALu;
	double H;
	double d_aux;
	double theta_aux;
	double dist_min = INFINITY;
	double tpi, tp_aux;
	double tsi, ts_aux;
	double ie_aux;

	for(j=0;j<size_ang;j++)
	{
		dTOTALu = 0.0;
		tpi=0.0;
		tsi=0.0;
		if(j==0)
			ang[j]=0.0;
		else
			ang[j]=ang[j-1]+d_ang;
		for(c=0;c<=indexs;c++)
		{
			if(c==0)
			{
				io[c]=ang[j];
			}
			else
			{
				io[c]=asin(vs[indexs-c]*sin(io[c-1])/vs[indexs+1-c]);
			}
			if(c==0)
			{
				H=fabs(zso)-aux_sum;
			}
			else
			{
				H=espesor[indexs-c];
			}

			dTOTALu=dTOTALu+H*tan(io[c]);
			tpi=tpi+H/(cos(io[c])*vp[indexs-c]);
            tsi=tsi+H/(cos(io[c])*vs[indexs-c]);
		}



		if(fabs(dTOTALu-Dhorizontal )<dist_min )
		{
			d_aux = dTOTALu;
			dist_min=fabs(dTOTALu-Dhorizontal );
			theta_aux = io[indexs];
			ie_aux=ang[j];
			tp_aux = tpi;
			ts_aux = tsi;
		}

	}

	*theta = theta_aux;
	*d_total = d_aux;
	*tp = tp_aux;
	*ts = ts_aux;
	*ie = M_PI-ie_aux; 
	free(espesor);
	free(ang);
	free(io);
	return;
}
