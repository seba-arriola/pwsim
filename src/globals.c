#include "pwsim.h"

//pthread for mutual exclusion in fftw3 plans and writing P times
pthread_mutex_t mutex;
pthread_mutex_t mutex2;

//File for P waves times
FILE *fppp;

//Number of points in the finite fault
int N;

//number of layers in velocity model
int M;

//number of stations
int V;

//calculate amplification factors due to the free surface
int calcfs;

//seismic moment and moment magnitude
double M0;
double Mw;

//Medium parameters
double alpha;
double beta;
double rho;

//Quality factors
double Qso,Qpo;
double Qexp;
int N_attpar;

//stress drop
double dsigma;

double *R_att;
double *R_att_aux_1;
double *p_att;

//hipocenter coordinates in meters, xhip=yhip=0.0
double xhip, yhip, zhip;

//corner frequency, according to Brune
double f0s, f0p;

//total time in seconds for waveforms
int waveform_time;

//Numbers of cores to use.
int n_cores;

int nfft;
int nnyq;
int Fs;
int T_S;

//Hypocenter geographical coordinates, used as reference
double ALATO;
double ALNGO;

//number of noise waveform to take an average
int N_simul;

//
int only_SH;

//apply transfer functions
int applyTF;
double rho_tf;
double damping_p,damping_sv,damping_sh;

//input filenames
char finite_fault[512];
char vel_model[512];
char stations_file[512];
char envelope_file[512];
char attenuation_file[512];
double envelope_e;
double envelope_n;
double envelope_ft;

//Fault parameters
double *X;
double *Y;
double *Z;
double *slip;
double *strike;
double *dip;
double *rake;
double *length;
double *width;
double *Tr;
double *trup;
double *m0;
double mean_dip;
double mean_rake;
double mean_length;
double mean_width;
double sum_slip;

//corner frequencies for every subfault, P and S waves
double *fcs;
double *fcp;
double *f;

//Radiation pattern to use
int radpat;

// Average radiation patterns according to Onishi Horike
double Rpp_OH;
double Rpsv_OH;
double Rpsh_OH;

//global velocity model parameters
double *prof;
double *vp;
double *vs;

//stations parameters
int sample;
double *ruido;
int *TFstat;
char **sta_names;
char **sta_models;
double *sta_x;
double *sta_y;
double *sta_z;
double *kappa;
double *gamma_sta;

//seed for randoms
int seed;
