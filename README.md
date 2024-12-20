## Open Parallel Waveform Simulation: Synthetic accelerogram 

Beta version, implemented and moderately tested. 

# 1.- Introduction

This C implementation was mostly programmed by Sebastian Arriola, any doubt, complaint, claim or bug report, please write to 

sarriola@csn.uchile.cl 

seba.arriola.s@gmail.com


For this project, we used some code obtained from the web, particularly from:

https://www.geeksforgeeks.org/

https://www.embeddedrelated.com/


The detail is mentioned within the code as a comment. 


The Gauss-Krueger projection algorithm was also implemented to convert latitude-longitude coordinates to X-Y and vice versa. 
This algorithm is described in the following link as Fortran code.

http://www.bosai.go.jp/study/application/dc3d/download/DC3Dmanual.pdf


****IMPORTANT****

This program is parallelised using the C pthread library. Each thread is used to calculate the synthetic waveforms of a station. 
The total computation time will be proportional to the number of stations divided by the number of threads.

# 2.- Installation

This version of the code only needs to be compiled to generate an executable. As a prerequisite, you need the FFTW library, which is used to perform Fourier transforms.

http://www.fftw.org/download.html

Once this library is installed, the code can be compiled using the Makefile, simply by using the command

>> make

The Makefile contains minimal compilation options to generate an executable. Other compilation options, e.g. for optimisation or use of another compiler, can be added. 
another compiler, can be added/modified at the user's discretion.

# 3.- Execution

Once compiled, if the compilation was successful, a file called pwsim is generated, which is executed with the command

./pwsim params_file


# 4.- Parameters

params_file corresponds to a file containing a series of parameters. An example for this parameter file is given below:


Mw 7.0 			(float)

ttime 100 		(int)

sps 100 		(int)

alpha 7.0 		(float)

beta 4.0 		(float)

rho 3.0 		(float)

dsigma 100.0 		(float)

lonhip 10.00 		(float)

lathip 10.00 		(float)

zhip 10.0 		(float)

threads 8 		(int)

ffm source.dat 		(file)

velmodel velmodel.dat 	(file)

stations stations.txt 	(file)

applyTF 0 		(int, option)

rho_tf 1.8 		(float)

b_p 0.05 		(float)

b_sv 0.05 		(float)

b_sh 0.05 		(float)

seed 0 			(int)

radpat 0 		(int, option)

calcfs 1 		(int, option)

N_simul 100 	(int)

only_SH 0 	(int, option)

envelope_file params/envelope.dat 	(file)

attenuation_file params/attenuation.dat 	(file)





Important to note that this file should only contain a variable name followed by its value.


Mw => Moment magnitude of the event giving rise to the waveforms.

ttime => Total synthetic trace time for the event. This value does not correspond to the total size of the generated trace. 

sps => Number of samples per second for the simulated trace. It must be an integer

alpha => Overall P-wave velocity in the vicinity of the source, in km/s 

beta => Overall S-wave velocity in the vicinity of the source, in km/s 

rho => Density in the vicinity of the source, in gr/cc

dsigma => Stress drop in bars

lonhip => Hypocentral longitude, used as a reference for conversion

lathip => Hypocentral latitude, used as a reference for conversion

zhip => Hypocentral depth in km

threads => Number of threads used to calculate the waveforms. 

ffm => File containing the finite source model. Each raw represents a sub-fault in which the position corresponds to the center of each sub-fault. This columnar file is ordered as follows:

longitude(°) latitude(°) depth(metres) slip(metres) strike(°) dip(°) rake(°) width(metres) length(metres) Break_time(sec).


velmodel => File containing the layered global velocity model. This file is sorted by columns as follows:

depth(metres) velocity_P(metres/sec) velocity_S(metres/sec) 

stations => File with details of the stations used to generate the waveforms. This file is sorted by columns as follows:


station_name longitude(°) latitude(°) elevation(metres) kappa_station gamma_station TF_model_path apply_function_transfer

- kappa_station corresponds to the empirical high-frequency attenuation factor in the vicinity of the station
- gamma_station corresponds to the spectral decay factor
- TF_model_path is a path to a file containing the information necessary to apply a soil amplification transfer function (SATF). If you want to use a theoretical model, you must specify the model of the layers above the station. This file should contain ONLY the layers above the station. The file should include the following information in columns:

thickness speed_P(metres/sec) speed_S(metres/sec)

Conversely, if you wish to apply pre-defined frequency amplification factors (or generic-rock soil amplification functions or GRSA), the file must contain your information according to the following format:

frequency horizontal_amplification vertical_amplification

- apply_transfer_function [if condition] "0" if you do not want to apply SATF, "1" if you want to apply the SATF, and "2" if you will apply GRSA.

applyTF => [if condition] If "1", SATF is applied to stations marked with 1 in the last value of the station file. "0" otherwise.

rho_tf => global density for SATF.

b_p => P-wave damping for SATF.

b_sv => SV-wave damping for SATF.

b_sh => SH-wave damping for SATF.

seed => Seed for white Gaussian noise generation. If a value of 0 is used, the seed is called with time(0), i.e., a new sequence of numbers is generated. If a numerical value is used (e.g., 11010220) you will generate accelerograms considering that value as the seed. Knowing the value to generate random noise allows you to reproduce results. The program will show you the seed used on the screen at each execution.

radpat => [if condition] Radiation pattern to be used for the calculations. If "0", a value is calculated for each sub-fault, which depends on the geometry of the sub-fault and the take-off and incidence angles of the source-station ray path (Aki & Richards, 2002). If "1", an average value over the entire fault is used (Onishi & Horike, 2004).

calcfs => [if condition] Calculates the free surface (FS) factors. If "0", you will use constant values FS=1. If "1", FS is calculated. If "2", both cases are done, and an acceleration file with six columns is written. You may find three components NS, EW, UD with FS factors calculated and applied (calcfs 1), then the other three components (NS, EW, UD) but without the FS factors applied (calcfs 0).

N_simul => number of simulations of noise' waveforms, to take an average.

only_SH => [if condition] If "1" only calculate SH component for accelaration, the other two columns are 0. If "0" calculate all three cartesian components (NS,EW,Z).

envelope_file => Contains e, n, ft, envelope parameters, each one in a different row.

attenuation_file => Contain attenuation parameters.

Qso Qpo Qexp
N
R0 R1 R2 R3 R4 ... RN-1
P1 P2 P3 P4 P5 ... PN


Qso => S-wave quality factor coefficient, Q(f) = Qso*f**Qexp

Qpo => P-wave quality factor coefficient, Q(f) = Qpo*f**Qexp



If everything is correct, the program runs and displays the total time used in seconds to generate the waveforms.


In the main directory where you execute the program, a directory called "output/" is created with one file per station, called "synt_STATION_NAME.dat" for the accelerations. These files have 3 columns, where the first one corresponds to the NS component, the second one to the EW component, and the third and last one corresponds to the UD component.
