import sys
import os
import time
import statsmodels.stats.multitest as smm
import csv
import copy
import numpy as np

os.environ['R_HOME'] = 'C:/R/R-3.5.1'
os.environ['R_USER'] = 'C:/Python/Python36/Lib/site-packages/rpy2'

from rpy2 import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import r

import rpy2.robjects.numpy2ri
ro.numpy2ri.activate()


def noBootstrap(aracne,expression_file,output_folder,gene_names,samples,message,percentage_of_average_for_threshold):
	start_time = time.time()
	message+='noBootstrap aracne is being executed.\n'

	try:
		result_p1 = os.popen('java -Xmx5G -jar '+aracne+' -e '+expression_file+' -o '+output_folder+' --tfs '+gene_names+' --pvalue 1E-8 --seed 1 --calculateThreshold')
		#print(result_p1)
		result_p1.close()
		threshold_file=open(output_foldername+"miThreshold_p1E-8_samples"+str(samples-1)+".txt",'r')	
		aracne_threshold=float(threshold_file.read())
		threshold_file.close()
		threshold_file=open(output_foldername+"miThreshold_p1E-8_samples"+str(samples-1)+".txt",'w')
		threshold_file.write(str(0.001))
		threshold_file.close()
		result_p2 = os.popen('java -Xmx5G -jar '+aracne+' -e '+expression_file+' -o '+output_folder+' --tfs '+gene_names+' --pvalue 1E-8 --seed 1')
		#print(result_p2)
		result_p2.close()
		dynamic_threshold_network(output_folder+"network.txt",percentage_of_average_for_threshold,aracne_threshold)
	
	except OSError as err:
		message+="os error: {0}".format(err) + "\n"
		return message
		
	except:
		message+="Unexpected error\n"
		return message

	message+='noBootstrap aracne has been successfully executed.\n'
	end_time = time.time()
	message+='Elapsed time was %g seconds.\n' % (end_time-start_time)
	return [message,aracne_threshold]



def bootstrap(aracne,expression_file,output_folder,gene_names,bootstraps,samples,message,percentage_of_average_for_threshold):
	start_time = time.time()	
	message+='Bootstrapped aracne is being executed.\n'
	print(samples)
	try:
		#Delete bootstrap files
		for filename in os.listdir(output_folder):
			if filename.startswith("bootstrap"):
				os.remove(output_folder+filename)
		result_p1 = os.popen('java -Xmx5G -jar '+aracne+' -e '+expression_file+' -o '+output_folder+' --tfs '+gene_names+' --pvalue 1E-8 --seed 1 --calculateThreshold')
		#print(result_p1)
		result_p1.close()
		threshold_file=open(output_foldername+"miThreshold_p1E-8_samples"+str(samples-1)+".txt",'r')	
		aracne_threshold=float(threshold_file.read())
		threshold_file.close()
		threshold_file=open(output_foldername+"miThreshold_p1E-8_samples"+str(samples-1)+".txt",'w')
		threshold_file.write(str(-1))
		threshold_file.close()
		
		for i in range(1,int(bootstraps)+1):
			result_pi = os.popen('java -Xmx5G -jar '+aracne+' -e '+expression_file+' -o '+output_folder+' --tfs '+gene_names+' --pvalue 1E-8 --seed '+str(i))
			result_pi.close()
			#print(result_pi)
		#apply dynamic threshold to all networks
		for filename in os.listdir(output_folder):
			if filename.startswith("bootstrap"):
				dynamic_threshold_network(output_folder+filename,percentage_of_average_for_threshold,aracne_threshold)
			
		
		result_p3 = os.popen('java -Xmx5G -jar '+aracne+' -o '+output_folder+' --consolidate --nobonferroni')
		#print(result_p3)
		result_p3.close()
		
	
	except OSError as err:
		
		message+="os error: {0}".format(err) + "\n"
		print(message)
		return message
		
	except:
		message+="Unexpected error\n"
		print(message)
		return message

	message+='Bootstrapped aracne has been successfully executed.\n'
	end_time = time.time()
	message+='Elapsed time was %g seconds.\n' % (end_time-start_time)
	return [message,aracne_threshold]

def dynamic_threshold_network(filename,dynamic_threshold,aracne_threshold):
	proteins=dict()
	interactions=list()
	num_of_line=0

	
	with  open(filename) as  fname:
		for line in csv.reader(fname, delimiter="\t"):
			if num_of_line>0:
				interactions.append([])
				interactions[num_of_line-1].append(line[0].strip())
				interactions[num_of_line-1].append(line[1].strip())
				interactions[num_of_line-1].append(float(line[2].strip()))
				if line[0].strip() in proteins and line[1].strip() in proteins:
					proteins[line[0].strip()].append(float(line[2].strip()))
					proteins[line[1].strip()].append(float(line[2].strip()))
				elif (line[0].strip() not in proteins) and (line[1].strip() in proteins):
					proteins[line[0].strip()]=list()
					proteins[line[0].strip()].append(float(line[2].strip()))
					proteins[line[1].strip()].append(float(line[2].strip()))
				elif (line[0].strip() in proteins) and (line[1].strip() not in proteins):
					proteins[line[1].strip()]=list()
					proteins[line[0].strip()].append(float(line[2].strip()))
					proteins[line[1].strip()].append(float(line[2].strip()))
				else:
					proteins[line[0].strip()]=list()
					proteins[line[1].strip()]=list()
					proteins[line[0].strip()].append(float(line[2].strip()))
					proteins[line[1].strip()].append(float(line[2].strip()))
				
			num_of_line+=1
	
	new_thresholds=dict()
	for key,vals in proteins.items():
		new_threshold=dynamic_threshold*(sum(vals)/float(len(vals)))
		if new_threshold>aracne_threshold:
			new_thresholds[key]=aracne_threshold
		else:
			new_thresholds[key]=new_threshold

	fname.close()
	os.remove(filename)

	fname=open(filename,"w")
	fname.write(str("Regulator")+"\t"+str("Target")+"\t"+str("MI")+"\n")
	for i in range(len(interactions)):
		if interactions[i][2]>=new_thresholds[interactions[i][0]]:
			fname.write(str(interactions[i][0])+"\t"+str(interactions[i][1])+"\t"+str(interactions[i][2])+"\n")
	fname.close()
	
def correction(pvals,corrected_pvalues,message):
	start_time = time.time()
	message+='Benjamini-Hockberg correction is being used for pvalues correction.\n'

	try:
		output = smm.multipletests(pvals,method='fdr_bh')
		corrected_pvalues=output[1]
		
	except:
		message+="Unexpected error\n"
		return message

	message+='pvalues have been succesfully corrected.\n'
	end_time = time.time()
	message+='Elapsed time was %g seconds.\n' % (end_time-start_time)
	return (message, corrected_pvalues)
		


def correctNetwork(corrected_pvalues,network_file,corrected_network_file,pval_threshold,message):
	start_time = time.time()
	message+='Corrected network is being constructed.\n'
	
	try:
		with open(network_file) as csvfile:
			
			reader = csv.DictReader(csvfile, delimiter='\t')
			
			with open(corrected_network_file, 'w') as newfile:
				fieldnames = ['Regulator', 'Target', 'MI', 'pvalue']
				newfile.write('Regulator'+'\t'+ 'Target'+'\t'+ 'MI'+'\t'+ 'pvalue'+'\n')
				#writer = csv.DictWriter(newfile, fieldnames=fieldnames, delimiter='\t')
				#writer.writeheader()
				for i,line in zip(range(len(corrected_pvalues)),reader):
					#print(str(line['Regulator'])+'\t'+str(line['Target'])+'\t'+str(line['MI'])+'\t'+str(line['pvalue'])+'\n')
					line['pvalue'] = corrected_pvalues[i]
					
					#print(line+"..")
					
						#writer.writerow(line)
					
					
				
					if float(line['pvalue'])<=float(pval_threshold):
						newfile.write(str(line['Regulator'])+'\t'+str(line['Target'])+'\t'+str(line['MI'])+'\t'+str(line['pvalue'])+'\n')
					
	except OSError as err:
		
		message+="os error: {0}".format(err) + "\n"
		
		return message
	except:
		message+='Unexpected error\n'
		
		return message

	message+='Corrected network has been succesfully constructed.\n'
	end_time = time.time()
	message+='Elapsed time was %g seconds.\n' % (end_time-start_time)
	return message




		
		
		 
def parse_expression_data(data_filename):
	num_of_lines=0
	molecules=dict()
	molecules2=dict()
	data=list()
	samples=list()
	
	with open(data_filename) as data_fname:
		for line in csv.reader(data_fname, delimiter="\t"):
			if num_of_lines==0:
				for j in range(len(line)):
					if j>0:
						samples.append(line[j].strip())
			else:
				molecules[line[0].strip()]=num_of_lines
				molecules2[num_of_lines]=line[0].strip()
				data.append([])
				for j in range(len(line)):
					if j>0:
						if line[j]!='':
							data[num_of_lines-1].append(float(line[j]))
						else:
							data[num_of_lines-1].append('')
			num_of_lines+=1
	print('Expression Data were successfully parsed!')
	return [molecules,molecules2,data,samples]


def parse_network_file(network_filename):
	num_of_lines=0
	network_edges=dict()
	
	
	
	with open(network_filename) as data_fname:
		for line in csv.reader(data_fname, delimiter="\t"):
			if num_of_lines==0:
				donothing=0;
			else:
				network_edges[(line[0].strip(),line[1].strip())]=(float(line[2].strip()),float(line[3].strip()))
				
			num_of_lines+=1
	print('Data were successfully parsed!')
	return network_edges

#siren_function(siren_folder_pathway, network_edges, molecules, expression_data )	
def siren_function(siren_folder_pathway,network_edges, molecules,molecules2, expression_data,output_folder,message):
	#create siren expression file
	siren_exp_datafile=open(siren_folder_pathway+'expression_data.txt','w')
	for i in range(len(expression_data)):
		for j in range(len(expression_data[i])):
			siren_exp_datafile.write(str(expression_data[i][j]))
			if j<(len(expression_data[i])-1):
				siren_exp_datafile.write("\t")
		if i<(len(expression_data)-1):
			siren_exp_datafile.write("\n")
	siren_exp_datafile.close()	
	#create sire network file
	siren_network_file=open(siren_folder_pathway+'network.txt','w')
	for key,val in network_edges.items():
		siren_network_file.write(str(molecules[key[0]]) +'\t'+ str(molecules[key[1]]) +'\n')
	siren_network_file.close()
	ro.r("source(\""+siren_folder_pathway+"SIREN.R\")")
	ro.r("exp=as.matrix(read.table(\""+siren_folder_pathway+"expression_data.txt\",sep=\"\t\"))")
	ro.r("net=as.matrix(read.table(\""+siren_folder_pathway+"network.txt\",sep=\"\t\"))")
	ro.r("w=as.matrix(read.table(\""+siren_folder_pathway+"Weighting_Matrix.txt\",sep=\"\t\"))")
	ro.r("Result=SIREN(exp,w,net)")
	ro.r("write.table(Result,file=\""+siren_folder_pathway+"Result.txt\",sep=\"\t\")")
	final_network_file=open(output_folder+"directed_network_file.txt",'w')
	num_of_lines=0
	with open(siren_folder_pathway+"Result.txt") as siren_output_file:
		for line in csv.reader(siren_output_file, delimiter="\t"):
			if num_of_lines==0:
				final_network_file.write("Regulator\tTarget\tMI Score\tPCorrected Pvalue\tSiren Directionality Score\n")
			else:
				final_network_file.write(str(molecules2[int(line[1].strip())])+'\t'+ str(molecules2[int(line[2].strip())])+"\t"+str(network_edges[(str(molecules2[int(line[1].strip())]),str(molecules2[int(line[2].strip())]))][0])+"\t"+str(network_edges[(str(molecules2[int(line[1].strip())]),str(molecules2[int(line[2].strip())]))][1])+"\t"+str(line[3].strip())+"\n")
			num_of_lines+=1
	final_network_file.close()
	message=message+"Siren was executed successfully.\n"
	return message

def filter_multiple_edges(output_folder,network_filename,message):
	interactions=dict()
	num_of_lines=0
	with open(network_filename) as fname:
		for line in csv.reader(fname, delimiter="\t"):
			if num_of_lines>0:
				interactions[(line[0].strip(),line[1].strip())]=list()
				interactions[(line[0].strip(),line[1].strip())].append(float(line[2].strip()))
				interactions[(line[0].strip(),line[1].strip())].append(float(line[3].strip()))
				interactions[(line[0].strip(),line[1].strip())].append(float(line[4].strip()))
			num_of_lines+=1
	new_interactions=copy.deepcopy(interactions)
	for key,vals in interactions.items():
		#remove duplicate edges between the same nodes when they have the same directionality
		if (key[1],key[0]) in interactions:
			if vals[2]*interactions[(key[1],key[0])][2]>0:
				if abs(vals[2])>abs(interactions[(key[1],key[0])][2]):
					if (((key[1],key[0])) in new_interactions) and (key in new_interactions):
						del new_interactions[(key[1],key[0])]
				else:
					if (((key[1],key[0])) in new_interactions) and (key in new_interactions):
						del new_interactions[key]
		#remove self loops
		if vals[0]==vals[1]:
			if (((key[1],key[0])) in new_interactions) and (key in new_interactions):
				del new_interactions[key]
	output_file=open(output_folder+"filtered_directed_network.txt",'w')
	output_file.write("Regulator\tTarget\tMI\tpvalue\tdirectionality\n")
	for key,vals in new_interactions.items():
		output_file.write(str(key[0])+"\t"+str(key[1])+"\t"+str(vals[0])+"\t"+str(vals[1])+"\t"+str(vals[2])+"\n")
	output_file.close()
	message=message+'Network has been successfully filtered.'
	return message
	
	
if __name__ == "__main__":
	jar_path=sys.argv[1]
	jar_path_modified=sys.argv[2]
	expression_filename=sys.argv[3]
	output_foldername=sys.argv[4]
	names_file=sys.argv[5]
	siren_folder_pathway=sys.argv[6]
	percentage_of_average_for_threshold=float(sys.argv[7]) #Percentage of the average correlations of each molecule to use for the dynamic threshold calculation
	pval_threshold=float(sys.argv[8]) 
	print('Pval_threshold='+str(pval_threshold))
	if len(sys.argv)>9:
		bootstraps=sys.argv[9]
	else:
		bootstraps='100'
	

	message=''

	f = open(expression_filename,'r')
	first_line = f.readline()
	samples = len(first_line.split('\t'))
	message+='Step 1: Running Aracne Ap with dynamic correlation thresholds...'
	#Step 1 Running Aracne Ap with dynamic correlation thresholds
	

	if samples>10:	
		[message,aracne_threshold]=bootstrap(jar_path,expression_filename,output_foldername,names_file,bootstraps,samples,message,percentage_of_average_for_threshold)


	else:
		[message,aracne_threshold]=noBootstrap(jar_path_modified,expression_filename,output_foldername,names_file,samples,message,percentage_of_average_for_threshold)
	
	
	#Step 2 Correcting Pvalues with Benjamini Hockberg Method
	message+='Step 2: Correcting Pvalues with Benjamini Hockberg Method...'
	p_vals=[]
	with open(os.path.join(sys.argv[4], 'network.txt')) as csvfile:
		reader = csv.DictReader(csvfile, delimiter='\t')
		data={}
		for row in reader:
			for header, value in row.items():
				try:
					data[header].append(value)
				except KeyError:
					message+='Key error: Accessing invalid key.\n'
					data[header] = [value]
				except:
					message+="Unexpected error\n"
		
		p_vals = data['pvalue'] #list of strings
	#print(p_vals) 
	pvals_str = np.array(p_vals) #string array
	pvals = pvals_str.astype(np.float) #float array
	#print(pvals)
	corrected_pvalues=[]
	message, corrected_pvalues=correction(pvals,corrected_pvalues,message)
	
	network_file = os.path.join(sys.argv[4], 'network.txt')
	corrected_network_file = os.path.join(sys.argv[4], 'corrected_network.txt')
	
	message = correctNetwork(corrected_pvalues,network_file,corrected_network_file,pval_threshold,message)
	#print(message)
	#genes = []                                          
	#expression_data = []                                
	#with open(expression_filename) as csvfile:
	#		reader = csv.reader(csvfile, delimiter='\t')
	#		headers = next(reader)
	#		for row in reader:
	#			genes.append(row[0])
	#			del row[0]
	#			expression_data.append(row[1:])
	#	print(genes)
	#	print(expression_data)
	
	#Step 3: Running Siren to calculate directionality
	message+='Step 3: Running Siren to calculate directionality...'
	#parse expression data
	[molecules,molecules2,expression_data,samples]=parse_expression_data(expression_filename)
	#parse network file
	network_edges=parse_network_file(corrected_network_file)
	message=siren_function(siren_folder_pathway, network_edges, molecules,molecules2, expression_data,output_foldername,message )
	#Step 4: Filtering multiple edges between the same node	
	message+='Step 4: Filtering multiple edges between the same nodes...'
	message=filter_multiple_edges(output_foldername,output_foldername+"directed_network_file.txt",message)
	log_file_sec=open('log_file_sec.txt','a')
	log_file_sec.write(message)
	log_file_sec.close()


