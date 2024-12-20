import sys
import os
import pandas as pd

#THIS CODE SHOULD BE USE LIKE 
# >> python3 check_input_params.py /path_to_params_file/params_file 
#
float_params = ["Mw","alpha" ,"beta","rho" ,"dsigma","lonhip" ,"lathip" ,"zhip" ,"b_p" ,"b_sv" ,"b_sh" ,"rho_tf" ]

int_params = ["applyTF","ttime","sps","threads","seed","radpat","calcfs","N_simul","only_SH"]

string_params = ["ffm","velmodel","stations","envelope_file","attenuation_file"]



#JB_file 0

def check_envelope_file(filename):
	envelope_dict=create_dictionary_from_file(filename)
	envelope_valids = ["e","n","ft"]
	for a_key in envelope_dict:
		if not a_key in envelope_valids:
			print("\nEnvelope file has an incorrect parameter named %s\n" %(a_key))
		else:
			envelope_valids.remove(a_key)
			try:
				asd=float(envelope_dict[a_key])
			except:
				print("\nValue %s associated to parameter %s from envelope file %s can not be converted to float\n" %(envelope_dict[a_key],a_key,filename))
	if len(envelope_valids)>0:
		print("\nEnvelope file %s has no %s parameter\n" %(filename,envelope_valids))
		
	

def check_attenuation_file(filename):
	data = [line_catalog.rstrip('\n') for line_catalog in open(filename)]
	if len(data)!=4:
		print("\nAttenuation file %s must contain four rows\n" %(filename))
	try:
		[Qpo,Qso,Qexp]=data[0].split()
		Qpo=float(Qpo)
		Qso=float(Qso)
		Qexp=float(Qexp)
	except:
		print("\nError at first line in attenuation file %s containing Qpo,Qso,Qexp\n" %(filename))
	
	N_attpar=0
	try:
		N_attpar = int(data[1])
	except:
		print("\nError at second line in attenuation file %s\n" %(filename))
		
	try:
		aux_R=data[2].split()
		if len(aux_R)<N_attpar:
			print("\nError at third line in attenuation file %s: Not enough parameters\n" %(filename))
	except:
		print("\nError at third line in attenuation file %s\n" %(filename))
		
	try:
		aux_p=data[3].split()
		if len(aux_p)<N_attpar:
			print("\nError at fourth line in attenuation file %s: Not enough parameters\n" %(filename))
	except:
		print("\nError at fourth line in attenuation file %s\n" %(filename))
		
	
	
	
	
	
	
def check_stations_file(filename):
	with open(filename) as f2:
		for line in f2:
			try:
				[station_name, longitude, latitude, elevation, kappa_station, gamma_station, TF_model_path, apply_function_transfer] = line.split()
				apply_function_transfer = apply_function_transfer.rstrip('\n')
				longitude=float(longitude)
				latitude=float(latitude)
				elevation=float(elevation)
				kappa_station=float(kappa_station)
				gamma_station=float(gamma_station)
				apply_function_transfer=int(apply_function_transfer)
			
			except:
				print("\nProblem with stations params at line %s\n" %(line))
				return
			
			
			if(apply_function_transfer==1):
				if  os.path.exists(TF_model_path):
					print(TF_model_path)
					df = pd.read_csv(TF_model_path, sep=" ", index_col=False, names=["depth", "velocity_P", "velocity_S"] )
					print("\nSTATISTICS FROM THE LOCAL MODEL FILE %s FOR STATION %s\n" %(TF_model_path,station_name) )
					print(df.describe() )
					print("\n\n")
				else:
					print("\n\nFile model %s does not exist or has incorrect format, you can not set apply transfer function parameter to 1 in station %s\n\n" %(TF_model_path,station_name))
			
			
			if(apply_function_transfer==2):
				if  os.path.exists(TF_model_path):
					print(TF_model_path)
					df = pd.read_csv(TF_model_path, sep=" ", index_col=False, names=['frequency', 'horizontal_amplification', 'vertical_amplification'] )
					print("\nSTATISTICS FROM THE LOCAL MODEL FILE %s FOR STATION %s\n" %(TF_model_path,station_name) )
					print(df.describe() )
					print("\n\n")
				else:
					print("\n\nFile model %s does not exist or has incorrect format, you can not set apply transfer function parameter to 1 in station %s\n\n" %(TF_model_path,station_name))
			
	

def analize_file(parname,filename):
	
	if parname=="ffm":
		try:
			df = pd.read_csv(filename, sep=" ", index_col=False,names=['longitude', 'latitude', 'depth', 'slip' ,'strike', 'dip' ,'rake', 'width', 'length', 'Tr'])
			print("\nSTATISTICS FROM THE FFM FILE %s\n" %(filename) )
			print(df.describe() )
			print("\n\n")
		except:
			print("\nSEEMS TO BE A PROBLEM WITH THE FFM FILE %s\n" %(filename) )
	elif parname=="velmodel":
		try:
			df = pd.read_csv(filename, sep=" ", index_col=False,names=["depth", "velocity_P", "velocity_S"])
			print("\nSTATISTICS FROM THE VELOCITY MODEL FILE %s\n" %(filename) )
			print(df.describe() )
			print("\n\n")
		except:
			print("\nSEEMS TO BE A PROBLEM WITH THE VELMODEL FILE %s\n" %(filename) )
	elif parname=="stations":
		check_stations_file(filename)
	elif parname=="envelope_file":
		check_envelope_file(filename)		
	elif parname=="attenuation_file":
		check_attenuation_file(filename)
	else:
		pass
		


def create_dictionary_from_file(filename):
	params_dictionary = { }
	with open(filename) as f2:
		for line in f2:
			[param_name,param_value] = line.split(" ") 
			param_value = param_value.rstrip('\n')
			params_dictionary[param_name] =  param_value
	return params_dictionary
	
	
def check_parameters_file(dict_params,filename):
	for a_key in dict_params:
		if a_key in float_params:
			try:
				float(dict_params[a_key])
			except:
				print( "\nValue %s from parameter %s couldn't be converted to float\n" %(dict_params[a_key],a_key) )
		elif a_key in int_params:
			try:
				int(dict_params[a_key])
			except:
				print( "\nValue %s from parameter %s couldn't be converted to int\n" %(dict_params[a_key],a_key) )
			
		elif a_key in string_params:
			if os.path.exists(dict_params[a_key]):
				analize_file(a_key,dict_params[a_key])				
			else:
				print( "\nFile %s from parameter %s doesn't exists\n" %(dict_params[a_key],a_key) )
			
		else:
			print("\nUnrecognized parameter %s in file %s\n" %(a_key,filename))
		
	
	
def main():
	params_filename = sys.argv[1]
	a_param_dict = create_dictionary_from_file(params_filename)
	check_parameters_file(a_param_dict,params_filename )
	
	

	#print(a_param_dict)

if __name__ == "__main__":
    main()
