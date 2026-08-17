#include "pwsim.h"

/* Initializes subfaults, global velocity model and stations parameters 
 * from files.
 * */
void readInputs()
{
	
	N=countLines(finite_fault);
	
	X=(double*)malloc(N*sizeof(double));
	Y=(double*)malloc(N*sizeof(double));
	Z=(double*)malloc(N*sizeof(double));
	slip=(double*)malloc(N*sizeof(double));
	strike=(double*)malloc(N*sizeof(double));
	dip=(double*)malloc(N*sizeof(double));
	rake=(double*)malloc(N*sizeof(double));
	length=(double*)malloc(N*sizeof(double));
	width=(double*)malloc(N*sizeof(double));
	Tr=(double*)malloc(N*sizeof(double));
	trup=(double*)malloc(N*sizeof(double));

	FILE	*archivo = fopen(finite_fault,"rt"); 
	if(archivo == NULL)
	{
		printf("Finite fault model file %s cannot be opened\n", finite_fault);
		exit(EXIT_FAILURE);
	}
	double A1,A2,A3,A4,A5, A6, A7,A8,A9,A10;
	int C = 0;
	double max_depth=0.0;
	double x_aux,y_aux;
	NR=0;
	sum_slip=0.0;
	mean_dip=0.0;
	mean_rake=0.0;
	mean_length=0.0;
    mean_width=0.0;
    while(fscanf(archivo,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf",&A1,&A2,&A3,&A4,&A5,&A6,&A7,&A8,&A9,&A10)==10)
    {
        geoToXY(&x_aux,&y_aux,A1,A2);
        Y[C]=x_aux;
		X[C]=-y_aux;
		Z[C]=-A3;
		if(A3>max_depth)
			max_depth=A3;
		slip[C]=A4;
		if(A4>0.0)
			NR=NR+1;
		sum_slip+=A4;
		strike[C]=A5*M_PI/180.0;
		dip[C]=A6*M_PI/180.0;
		mean_dip+=A6;
		rake[C]=A7*M_PI/180.0;
		mean_rake+=A7;
		width[C]=A8;
        length[C]=A9;
        mean_length=mean_length+length[C];
        mean_width=mean_width+width[C];
        Tr[C]=A10;
		trup[C]=A10;
        C++;
    }

	fclose(archivo);

	if(C < N)
	{
		printf("Warning: finite fault file %s: read %d of %d subfaults (malformed/truncated lines ignored)\n", finite_fault, C, N);
		N = C;
	}
	if(N == 0)
	{
		printf("Error: no subfaults read from finite fault file %s\n", finite_fault);
		exit(EXIT_FAILURE);
	}
	if(sum_slip == 0.0)
	{
		printf("Error: total slip is zero in finite fault file %s\n", finite_fault);
		exit(EXIT_FAILURE);
	}

	mean_dip=mean_dip/N;
	mean_rake=mean_rake/N;
	mean_length=mean_length/N/1000.; // en km
    mean_width=mean_width/N/1000.;


    M=countLines(vel_model);
    prof=(double*)malloc(M*sizeof(double));
	vp=(double*)malloc(M*sizeof(double));
	vs=(double*)malloc(M*sizeof(double));
    
    archivo = fopen(vel_model,"rt");
    if(archivo == NULL)
    {
		printf("Velocity model file %s cannot be opened\n", vel_model);
		exit(EXIT_FAILURE);
	}
    C=0; 
    
    while(fscanf(archivo,"%lf  %lf  %lf", &A1,&A2,&A3)==3)
    {
        prof[C]=A1;
		vp[C]=A2;
		vs[C]=A3;
		
        C++;
    }
    
    
    fclose(archivo);

	if(C < M)
	{
		printf("Warning: velocity model file %s: read %d of %d layers\n", vel_model, C, M);
		M = C;
	}
	if(M == 0)
	{
		printf("Error: no layers read from velocity model file %s\n", vel_model);
		exit(EXIT_FAILURE);
	}
    
    
    V=countLines(stations_file);
 
    sta_names = (char**)malloc(V*sizeof(char*));
	sta_models = (char**)malloc(V*sizeof(char*));
	sta_x = (double*)malloc(V*sizeof(double));
	sta_y = (double*)malloc(V*sizeof(double));
	sta_z = (double*)malloc(V*sizeof(double));
	kappa = (double*)malloc(V*sizeof(double));
	gamma_sta = (double*)malloc(V*sizeof(double));

	
	TFstat = (int*)malloc(V*sizeof(int));
	for(int i=0;i<V;i++)
	{
		sta_names[i]  = NULL;
		sta_models[i] = NULL;
	}
	
	char A12[256];
	char A11[512];
	int A14;

	archivo = fopen(stations_file,"rt");
	
    if(archivo == NULL)
    {
		printf("Stations file %s cannot be opened\n", stations_file);
		exit(EXIT_FAILURE);
	}
    C=0; 
    
    while(fscanf(archivo,"%255s %lf %lf %lf %lf %lf %511s %d", A12,&A2,&A3,&A4,&A5, &A6,A11,&A14)==8)
    {
        sta_names[C] = (char*)malloc(strlen(A12)+1);
        strcpy(sta_names[C], A12);
		printf("sta_name: %s \n",sta_names[C]);

        sta_models[C] = (char*)malloc(strlen(A11)+1);
        strcpy(sta_models[C],A11);

		geoToXY(&x_aux,&y_aux,A2,A3);
        sta_y[C]=x_aux;
		sta_x[C]=-y_aux;
		sta_z[C]= A4;
		kappa[C]= A5;
		gamma_sta[C]= A6;
		TFstat[C]=A14;
        C++;
    }
	fclose(archivo);

	if(C < V)
	{
		printf("Warning: stations file %s: read %d of %d stations\n", stations_file, C, V);
		V = C;
	}
	if(V == 0)
	{
		printf("Error: no stations read from stations file %s\n", stations_file);
		exit(EXIT_FAILURE);
	}
	
	int noise_size = (int)(sample*waveform_time);
	ruido = (double*)malloc(noise_size*sizeof(double));
	
	for(C=0;C<noise_size;C++)
    {//gaussian white noise
        ruido[C]=0.0;
    }
        
    if(seed==0)
		seed=time(0);
    
    srand(seed);
    printf("seed: %d\n",seed);
    
    
	for(int nnoise=0;nnoise<N_simul;nnoise++)
    {
     
        for(C=0;C<noise_size;C++)
        {//gaussian white noise averaged over N_simul simulations
            ruido[C] += AWGN_generator();
        }
        seed = seed + nnoise*(C+1);
        srand(seed);
        
    }

	if(N_simul>1)
	{
		for(C=0;C<noise_size;C++)
			ruido[C]/=N_simul;
	}

	return;

}


//thread function to calculate waveforms in different cores
void *generateWaveform_thread(void *ptr)
{
	Args *args= (Args*)ptr;
	CalculateWaveform(args->starting, args->ending);
	return NULL;
}



void readParams(int argcasd, char *argvasd)
{

	//checking the input parameters file
	if(argcasd!=2)
	{   
		printf("ERROR: Usage: ./exec_file <input parameters file> \n");
		exit(EXIT_FAILURE);
    }
	else
	{
		//reading input parameters file
		removeEmptyLines(argvasd);
		FILE * fr = fopen(argvasd, "rt");
		if(fr == NULL){
			printf("Parameters file %s not found\n", argvasd);
			exit(EXIT_FAILURE);
		}

		char tmpstr1[256] ;
		char tmpstr2[256] ;

		char tempbuff[1024];
		
		//reading every parameter
		while(fgets(tempbuff,sizeof(tempbuff),fr))
		{
			
			if (tempbuff[0]=='#' || tempbuff[0]=='/' || tempbuff[0]=='\n' || tempbuff[0]=='\r')
				continue;
			
			if (sscanf(tempbuff, "%255s %255s", tmpstr1, tmpstr2) != 2)
				continue;
			
			
			if (strcmp(tmpstr1,"Mw")==0)
			{
				Mw = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"alpha")==0)
			{
				alpha = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"beta")==0)
			{
				beta = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"rho")==0)
			{
				rho = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"dsigma")==0)
			{
				dsigma = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"lathip")==0)
			{
				ALATO = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"lonhip")==0)
			{
				ALNGO = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"zhip")==0)
			{
				zhip = -1000.0*atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"ffm")==0)
			{
				strcpy(finite_fault,tmpstr2);
			}
			else if (strcmp(tmpstr1,"velmodel")==0)
			{
				strcpy(vel_model,tmpstr2);
			}	
			else if (strcmp(tmpstr1,"threads")==0)
			{
				n_cores = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"stations")==0)
			{
				strcpy(stations_file,tmpstr2);
			}
			else if (strcmp(tmpstr1,"ttime")==0)
			{
				waveform_time = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"sps")==0)
			{
				sample = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"applyTF")==0)
			{
				applyTF = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"b_p")==0)
			{
				am_p = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"b_sv")==0)
			{
				am_sv = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"b_sh")==0)
			{
				am_sh = atof(tmpstr2);
			}
			else if (strcmp(tmpstr1,"rho_tf")==0)
			{
				rho_tf = atof(tmpstr2)*1000.0;
			}
			else if (strcmp(tmpstr1,"seed")==0)
			{
				seed = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"radpat")==0)
			{
				radpat = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"calcfs")==0)
			{
				calcfs = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"N_simul")==0)
			{
				N_simul = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"only_SH")==0)
			{
				only_SH = atoi(tmpstr2);
			}
			else if (strcmp(tmpstr1,"envelope_file")==0)
			{
				strcpy(envelope_file,tmpstr2);
			}
			else if (strcmp(tmpstr1,"attenuation_file")==0)
			{
				strcpy(attenuation_file,tmpstr2);
			}
			else
			{
				printf("Unrecognized parameter : \"%s\"\n", tmpstr1);
				exit(EXIT_FAILURE);
			}
		}
		
		fclose(fr);
		
	}
	
}



void readEnvelope()
{
	//reading input parameters file
	removeEmptyLines(envelope_file);
	FILE * fr = fopen(envelope_file, "r");
	if(fr == NULL)
	{
		//Boore (2003) parameters
		printf("Envelope parameters file %s not found, using default values\ne=0.2\nn=0.05\nft=2.0\n", envelope_file);
		envelope_e = 0.2;   //input parameter
		envelope_n = 0.05;  //input parameter
		envelope_ft = 2.;
		return;
	}

	char tmpstr1[256] ;
	char tmpstr2[256] ;

	char tempbuff[1024];
	
	//reading every parameter
	while(fgets(tempbuff,sizeof(tempbuff),fr))
	{
		
		if (tempbuff[0]=='#' || tempbuff[0]=='/' || tempbuff[0]=='\n' || tempbuff[0]=='\r')
			continue;
		
		if (sscanf(tempbuff, "%255s %255s", tmpstr1, tmpstr2) != 2)
			continue;
		
		if (strcmp(tmpstr1,"e")==0)
		{
			envelope_e = atof(tmpstr2);
		}
		else if (strcmp(tmpstr1,"n")==0)
		{
			envelope_n = atof(tmpstr2);
		}
		else if (strcmp(tmpstr1,"ft")==0)
		{
			envelope_ft = atof(tmpstr2);
		}				
		else
		{
			printf("Unrecognized parameter for envelope: \"%s\"\n using default values\ne=0.2\nn=0.05\nft=2.0\n", tmpstr1);
			envelope_e = 0.2;   //input parameter
			envelope_n = 0.05;  //input parameter
			envelope_ft = 2.;
		}
	}
	
	fclose(fr);
	
}


//read attenuation file
void readAttenuation()
{
    FILE *fp = fopen(attenuation_file, "r");
    if (!fp) {
        fprintf(stderr, "Error opening attenuation file '%s'\n", attenuation_file);
        exit(EXIT_FAILURE);
    }

    char line[1024];
    char *p;

    // first line: Qpo, Qso, Qexp
    if (fgets(line, sizeof(line), fp) == NULL ||
        sscanf(line, "%lf %lf %lf", &Qpo, &Qso, &Qexp) != 3) {
        fprintf(stderr, "Error: first line of attenuation file '%s' must contain Qpo Qso Qexp\n", attenuation_file);
        exit(EXIT_FAILURE);
    }

    // second line: number of elements
    if (fgets(line, sizeof(line), fp) == NULL ||
        sscanf(line, "%d", &N_attpar) != 1 || N_attpar <= 0) {
        fprintf(stderr, "Error: second line of attenuation file '%s' must contain a positive integer N\n", attenuation_file);
        exit(EXIT_FAILURE);
    }
    if (N_attpar > 1000) {
        fprintf(stderr, "Error: attenuation file '%s' declares too many parameters (N=%d, max 1000)\n", attenuation_file, N_attpar);
        exit(EXIT_FAILURE);
    }

    R_att=(double *)malloc( N_attpar*sizeof(double ));
    R_att_aux_1=(double *)malloc( (N_attpar+1)*sizeof(double ));
    p_att=(double *)malloc( N_attpar*sizeof(double ));

    int i;

    // third line: R0 R1 ... RN-1
    if (fgets(line, sizeof(line), fp) == NULL) {
        fprintf(stderr, "Error: missing R parameters in attenuation file '%s'\n", attenuation_file);
        exit(EXIT_FAILURE);
    }
    p = line;
    for (i=0; i<N_attpar; i++) {
        char *tok = strtok_r(p, " \t", &p);
        if (tok == NULL) {
            fprintf(stderr, "Error: not enough R parameters in attenuation file '%s'\n", attenuation_file);
            exit(EXIT_FAILURE);
        }
        R_att[i]=atof(tok);
    }

    // fourth line: P1 P2 ... PN
    if (fgets(line, sizeof(line), fp) == NULL) {
        fprintf(stderr, "Error: missing P parameters in attenuation file '%s'\n", attenuation_file);
        exit(EXIT_FAILURE);
    }
    p = line;
    for (i=0; i<N_attpar; i++) {
        char *tok = strtok_r(p, " \t", &p);
        if (tok == NULL) {
            fprintf(stderr, "Error: not enough P parameters in attenuation file '%s'\n", attenuation_file);
            exit(EXIT_FAILURE);
        }
        p_att[i]=atof(tok);
        R_att_aux_1[i]=R_att[i];
    }

    R_att_aux_1[0]=0.0;
    R_att_aux_1[N_attpar]=999999999.9;

    fclose(fp);
    return;
}






int main(int argc, char **argv)
{
	int i;
	
	//time measurement
	time_t tbegin, tend;
	
	tbegin = time(NULL);
	
	if(argc != 2)
	{
		printf("ERROR: Usage: ./pwsim <input parameters file>\n");
		exit(EXIT_FAILURE);
	}
	readParams(argc,argv[1]);
	readEnvelope();
	readAttenuation();

	if(n_cores < 1)
	{
		printf("ERROR: 'threads' must be a positive integer\n");
		exit(EXIT_FAILURE);
	}
	if(sample <= 0 || waveform_time <= 0)
	{
		printf("ERROR: 'sps' and 'ttime' must be positive integers\n");
		exit(EXIT_FAILURE);
	}
	if(N_simul < 1)
		N_simul = 1;
	
	xhip=0.0;
	yhip=0.0;
	
	
	struct stat st = {0};
	
	//create output directory
	if (stat("output", &st) == -1)
		mkdir("output", 0777);
	
	
	readInputs();
	
	calc();	
	
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&mutex2, NULL);
	
	fppp = fopen("output/pwaves_times.dat","w");
	if(fppp == NULL)
	{
		printf("Cannot write p waves time to file pwaves_times.dat\n");
		exit(EXIT_FAILURE);
	}
	
	int C = V/n_cores;
	int D = V%n_cores;
	
	if(C==0)
		n_cores=D;
	
	Args *args_a= (Args*)malloc(n_cores*sizeof(Args));
	
	int aux[n_cores];
	
	for(i=0;i<n_cores;i++)
		aux[i]=C;

	if( D>0 )
		for(i=0;i<D;i++)
			aux[i] = aux[i] + 1;
	
	int ind[n_cores+1], ini[n_cores], end[n_cores];
	ind[0] = 0;

	for(i=1;i<n_cores+1;i++)
	{
		ind[i] = aux[i-1]+ind[i-1];
		ini[i-1] = ind[i-1];
		end[i-1] = ind[i]-1;
	}
	
	Args *args;     
	
	for (i= 0; i<n_cores; i++)
	{
		args= &args_a[i];
		args->starting = ini[i];
		args->ending = end[i];
		pthread_create(&args->pid, NULL, generateWaveform_thread, args);
    }
	
	for (i= 0; i<n_cores; i++)
	{
		args= &args_a[i];
		pthread_join(args->pid, NULL);
	}
	
    fclose(fppp);
    
    //cleaning
    free(args_a);
	
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&mutex2);
	free(X);
	free(Y);
	free(Z);
	free(slip);
	free(strike);
	free(dip);
	free(rake);
	free(length);
	free(width);
	free(Tr);
	free(trup);
	free(prof);
	free(vp);
	free(vs);
	free(m0);
	free(fcs);
	free(fcp);
	free(sta_y);
	free(sta_x);
	free(sta_z);
	free(ruido);
	free(kappa);
	free(gamma_sta);
	free(TFstat);
	free(f);
	free(R_att);
	free(R_att_aux_1);
	free(p_att);

	
	for(i=0; i< V; i++)
	{
		free(sta_names[i]);
		free(sta_models[i]);
	}
	free(sta_names);
	free(sta_models);
	
	
	tend = time(NULL);
	int time_spent = (int)(tend - tbegin);
	
	printf("Total time = %d seconds\n",time_spent);
	
	malloc_trim(0);
	
	return 0;
}
