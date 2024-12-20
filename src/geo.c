#include "pwsim.h"

/* Gauss-Kruger's method to convert lat-lon to x-y
 * Receives two pointers x and y to store the x-y calculated coordinates  
 * and the longitude lon_x and latitude lat_y to be converted
 * */
void geoToXY(double *x, double *y, double lon_x, double lat_y)
{

	double A=6378.160, E2=6.6946053*pow(10.0,-3.0), E12=6.7397251*pow(10.0,-3.0), D=57.29578, RD=1.0/57.29578;
	double RLAT = lat_y*RD;
	double SLAT = sin( RLAT );
	double CLAT = cos( RLAT );
	double V2 = 1.0 + E12*CLAT*CLAT;
	double AL = lon_x - ALNGO;
	double PH1 = lat_y + 0.5*V2*AL*AL*SLAT*CLAT*RD;
	double RPH1 = PH1*RD;
	double RPH2 = (PH1 + ALATO)/2.0*RD;
	double SRPH1 = sin( RPH1 );
	double SRPH2 = sin( RPH2 );
	double R = A*(1.0 - E2) /pow(sqrt( 1. - E2*SRPH2*SRPH2 ),3.0);
	double AN = A/sqrt( 1.0 - E2*SRPH1*SRPH1 );
	double C1 = D/R;
	double C2 = D/AN;
	*x = 1000.0*AL*CLAT/C2*( 1.0 + AL*AL*cos(2.0*RLAT)/(6.0*D*D) );
	*y = 1000.0*(PH1 - ALATO)/C1;

	return;
	
}

/* Gauss-Kruger's method to convert x-y to lat-lon
 * Receives two pointers ALNGDG and ALATDG to store the lon-lat 
 * calculated coordinates and the longitude lon_x and latitude lat_y
 * to be converted
 * */
void xyToGEO(double *ALNGDG, double *ALATDG, double x, double y)
{
	double y_aux=y/1000.0;
	double x_aux =x/1000.0;
	double A=6378.160, E2=6.6946053*pow(10.0,-3.0), E12=6.7397251*pow(10.0,-3.0), D=57.29578, RD=1.0/57.29578;
	
	double RLATO = ALATO*RD;
	double SLATO = sin( RLATO );
	double CLATO = cos( RLATO );
	double DEN = sqrt( 1. - E2*SLATO*SLATO );
	double R = A*( 1. - E2 ) /(DEN*DEN*DEN);
	double AN = A / DEN;
	double V2 = 1. + E12*CLATO*CLATO;
	
	double C1 = D / R;
	double C2 = D / AN;
	double PH1 = ALATO + C1*y_aux;
	double RPH1 = PH1*RD;
	double TPH1 = tan(RPH1);
	double CPH1 = cos(RPH1);
	double BL = C2*x_aux;
	*ALATDG = PH1 - 0.5*BL*BL*V2*TPH1*RD;
	*ALNGDG = ALNGO+BL/CPH1*(1.- BL*BL*(1.+2.*TPH1*TPH1)/(6.*D*D));

	return;
	
}
