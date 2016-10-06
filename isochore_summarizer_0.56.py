#!/usr/bin/python

##  A script by Alex Hall
##  allopatry@gmail.com
##  August 22, 2016
##  Version 0.56

##  Usage:
##    >>> python isochore_summarizer.py [absolute path to folder containing output of isochore_analyzer.py]
##  Iteratively reads the log files from the isochore_analyzer1.831.py script
##  Script produces output for each of the following categories, separated by tabs. A header row with these variable names is created.
##  Sample	Class	fasta_number	fasta_recovered_percent_1000bp	fasta_ignored_missing_percent_1000bp	fasta_ignored_small_percent_1000bp	fasta_ignored_total_percent_1000bp	mean_1000bp	sd_1000bp	range_low_1000bp	range_high_1000bp	median_gc_1000bp	contigs_analyzed_1000bp	fasta_recovered_percent_3000bp	fasta_ignored_missing_percent_3000bp	fasta_ignored_small_percent_3000bp	fasta_ignored_total_percent_3000bp	mean_3000bp	sd_3000bp	range_low_3000bp	range_high_3000bp	median_gc_3000bp	contigs_analyzed_3000bp	fasta_recovered_percent_5000bp	fasta_ignored_missing_percent_5000bp	fasta_ignored_small_percent_5000bp	fasta_ignored_total_percent_5000bp	mean_5000bp	sd_5000bp	range_low_5000bp	range_high_5000bp	median_gc_5000bp	contigs_analyzed_5000bp	fasta_recovered_percent_20000bp	fasta_ignored_missing_percent_20000bp	fasta_ignored_small_percent_20000bp	fasta_ignored_total_percent_20000bp	mean_20000bp	sd_20000bp	range_low_20000bp	range_high_20000bp	median_gc_20000bp	contigs_analyzed_20000bp	fasta_recovered_percent_80000bp	fasta_ignored_missing_percent_80000bp	fasta_ignored_small_percent_80000bp	fasta_ignored_total_percent_80000bp	mean_80000bp	sd_80000bp	range_low_80000bp	range_high_80000bp	median_gc_80000bp	contigs_analyzed_80000bp	fasta_recovered_percent_320000bp	fasta_ignored_missing_percent_320000bp	fasta_ignored_small_percent_320000bp	fasta_ignored_total_percent_320000bp	mean_320000bp	sd_320000bp	range_low_320000bp	range_high_320000bp	median_gc_320000bp	contigs_analyzed_320000bp

##  Missing values will be recorded as NA; thus, values of 0.0 are real.
##  Output is placed in the current working directory.

##  Example file hierarchy:
##
##  /
##    data/
##      fish/
##        Run1.831_GC_output_20140520_1000_1468015719/
##          gc_content_20140520_1000_1468015719.tsv
##          log_20140520_1000_1468015719.txt
##          unused_contigs_missing_data_20140520_1000_1468015719.tsv
##          unused_short_contigs_content_20140520_1000_1468015719.tsv
##        Run1.831_GC_output_20140520_5000_1467929224/
##          gc_content_20140520_5000_1467929224.tsv
##          log_20140520_5000_1467929224.txt
##          unused_contigs_missing_data_20140520_5000_1467929224.tsv
##          unused_short_contigs_content_20140520_5000_1467929224.tsv
##
##  You would type:
##    python isochore_summarizer.py /data/fish
##  and the script will analyze each file beginning with "log_" in each folder within "/data/fish/"

##  Log file must be in the following format:
##
##    Script started at: 2016-07-07 17:07:04.964499 
##    The file analyzed was: /media/whiptail/Drive3/Alex_Isochores/Fish/GCA_001465895.2_Nfu_20140520_genomic.fna2 
##    The file was a valid fasta file format. 
##    Passed physical memory check. 
##    Window size provided or assumed: 5000 (i.e., 5 kb).
##     5896 fasta headers (i.e., contigs) were detected.
##     As a percent of all contigs recovered, 58.92147983 percent could have their GC content analyzed.
##     As a percent of all contigs recovered, 15.7458371072 percent could not have their GC content analyzed because there was too much ambiguous or missing data. 
##    As a percent of all contigs recovered, 25.3326830628 percent could not have their GC content analyzed because the contig was not larger than or equal to 5000 base pairs.
##    
##    41.07852017 percent of contigs could not be used to calculate GC content.
##     Warning! The amount of missing data exceeded 10 percent as measured by the percentage of contigs that could be used. Be sure this is okay for the interpretation of your results.
##    
##    Summary of results from isochore analysis by non-overlapping 5000 base pair windows:
##     Mean GC content: 43.0004072187 percent.
##     Standard deviation of GC content: 3.64256753304 percent.
##     GC content ranges from 21.88 percent to 60.2884420069 percent.
##     Median GC content: 42.7107591989 percent.
##     Overall, there were 8457 contigs analyzed.
##     Script finished at: 2016-07-07 17:07:34.246061 
##     29.16642 seconds elapsed running this script (i.e., 0.49 minutes).

##  Setup environment:
from  __future__ import print_function
from glob import glob
from collections import OrderedDict
import sys, re, os, time, datetime, errno, math, numpy
inputpath_o, path = sys.argv[1], os.getcwd()
d = {}
inputpath = inputpath_o.rstrip('/')
inputpath_split = inputpath.split('/')
Class = inputpath_split[-1]

##  Record the start time:
starttime = time.clock()
timesinceepoch = int(time.time())
now = datetime.datetime.now()
print('\n======================================\n\nTime at the beginning of the script:', now)

if inputpath.startswith("./"):
	sys.exit('Error: for stability, provide an absolute path to the data and not a relative path.')

##Set up an output file and open it:
outfile_path=path+'/'+'isochore_summary_'+Class+'_'+str(timesinceepoch)+'.tsv'
outfile = open(outfile_path, 'a')
outfile.write("Sample\tClass\tfasta_number\tfasta_recovered_percent_1000bp\tfasta_ignored_missing_percent_1000bp\tfasta_ignored_small_percent_1000bp\tfasta_ignored_total_percent_1000bp\tmean_1000bp\tsd_1000bp\trange_low_1000bp\trange_high_1000bp\tmedian_gc_1000bp\tcontigs_analyzed_1000bp\tfasta_recovered_percent_3000bp\tfasta_ignored_missing_percent_3000bp\tfasta_ignored_small_percent_3000bp\tfasta_ignored_total_percent_3000bp\tmean_3000bp\tsd_3000bp\trange_low_3000bp\trange_high_3000bp\tmedian_gc_3000bp\tcontigs_analyzed_3000bp\tfasta_recovered_percent_5000bp\tfasta_ignored_missing_percent_5000bp\tfasta_ignored_small_percent_5000bp\tfasta_ignored_total_percent_5000bp\tmean_5000bp\tsd_5000bp\trange_low_5000bp\trange_high_5000bp\tmedian_gc_5000bp\tcontigs_analyzed_5000bp\tfasta_recovered_percent_20000bp\tfasta_ignored_missing_percent_20000bp\tfasta_ignored_small_percent_20000bp\tfasta_ignored_total_percent_20000bp\tmean_20000bp\tsd_20000bp\trange_low_20000bp\trange_high_20000bp\tmedian_gc_20000bp\tcontigs_analyzed_20000bp\tfasta_recovered_percent_80000bp\tfasta_ignored_missing_percent_80000bp\tfasta_ignored_small_percent_80000bp\tfasta_ignored_total_percent_80000bp\tmean_80000bp\tsd_80000bp\trange_low_80000bp\trange_high_80000bp\tmedian_gc_80000bp\tcontigs_analyzed_80000bp\tfasta_recovered_percent_320000bp\tfasta_ignored_missing_percent_320000bp\tfasta_ignored_small_percent_320000bp\tfasta_ignored_total_percent_320000bp\tmean_320000bp\tsd_320000bp\trange_low_320000bp\trange_high_320000bp\tmedian_gc_320000bp\tcontigs_analyzed_320000bp\n")
outfile_ggplot_path=path+'/'+'isochore_summary_ggplot_'+Class+'_'+str(timesinceepoch)+'.tsv'
outfile_ggplot = open(outfile_ggplot_path, 'a')
outfile_ggplot.write("Sample\tClass\tfasta_number\twindow_size\tfasta_recovered_percent\tfasta_ignored_missing_percent\ttfasta_ignored_small_percent\tfasta_ignored_total_percent\tmean_gc\tsd_gc\trange_low_gc\trange_high_gc\tmedian_gc\tcontigs_analyzed\n")


#Compile list of samples:
os.chdir(inputpath)
paths = [p.rstrip('/') for p in glob('*/')]  #Suggested here: http://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python#comment63926589_18278257
for p in paths:
	sub_string = str(p)
	if sub_string.startswith('Run'):
		sub_string_sample1 = re.sub('^.*Run1.831_GC_output_', '', sub_string) #Removes beginning of filename
		sub_string_sample2 = re.sub('_\d*_\d*$', '', sub_string_sample1) #Removes timestamp and window size
		if not sub_string_sample2 in d.keys():
			d[sub_string_sample2] = 1
		else:
			d[sub_string_sample2] = d[sub_string_sample2]+1

#Iterate through every sample in the directory you provided. Pulls out relevant info from the log file and prints it one one line to the output file.
for key in d:
	#Initialize these variables with default of "NA" for each value:
	fasta_recovered_percent_1000bp, fasta_ignored_missing_percent_1000bp, fasta_ignored_small_percent_1000bp, fasta_ignored_total_percent_1000bp, mean_1000bp, sd_1000bp, range_low_1000bp, range_high_1000bp, median_gc_1000bp, contigs_analyzed_1000bp, fasta_recovered_percent_3000bp, fasta_ignored_missing_percent_3000bp, fasta_ignored_small_percent_3000bp, fasta_ignored_total_percent_3000bp, mean_3000bp, sd_3000bp, range_low_3000bp, range_high_3000bp, median_gc_3000bp, contigs_analyzed_3000bp, fasta_recovered_percent_5000bp, fasta_ignored_missing_percent_5000bp, fasta_ignored_small_percent_5000bp, fasta_ignored_total_percent_5000bp, mean_5000bp, sd_5000bp, range_low_5000bp, range_high_5000bp, median_gc_5000bp, contigs_analyzed_5000bp, fasta_recovered_percent_20000bp, fasta_ignored_missing_percent_20000bp, fasta_ignored_small_percent_20000bp, fasta_ignored_total_percent_20000bp, mean_20000bp, sd_20000bp, range_low_20000bp, range_high_20000bp, median_gc_20000bp, contigs_analyzed_20000bp, fasta_recovered_percent_80000bp, fasta_ignored_missing_percent_80000bp, fasta_ignored_small_percent_80000bp, fasta_ignored_total_percent_80000bp, mean_80000bp, sd_80000bp, range_low_80000bp, range_high_80000bp, median_gc_80000bp, contigs_analyzed_80000bp, fasta_recovered_percent_320000bp, fasta_ignored_missing_percent_320000bp, fasta_ignored_small_percent_320000bp, fasta_ignored_total_percent_320000bp, mean_320000bp, sd_320000bp, range_low_320000bp, range_high_320000bp, median_gc_320000bp, contigs_analyzed_320000bp = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
	#Defaults these to "None" in case the analysis did not produce a log file output
	file_1000bp, file_3000bp, file_5000bp, file_20000bp, file_80000bp, file_320000bp = None, None, None, None, None, None
	#Iterative variable:
	Sample = key
	try:  ##Creates an empty list [] if it can't find anything.
		file_1000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_1000_*/log*')
	except Exception:
		pass
	try:
		file_3000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_3000_*/log*')
	except Exception:
		pass
	try:
		file_5000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_5000_*/log*')
	except Exception:
		pass
	try:
		file_20000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_20000_*/log*')
	except Exception:
		pass
	try:
		file_80000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_80000_*/log*')
	except Exception:
		pass
	try:
		file_320000bp = glob(inputpath+'/'+'Run1.831_GC_output_'+key+'_320000_*/log*')
	except Exception:
		pass
	if len(file_1000bp)==0:  ##Checks to see that a file was found. If there was no log file found, the result of len(file_1000bp) will be 0 and this loop will be skipped.
		print('No log file for', key, 'at 1000bp. Skipping...')
		pass
	else:
		with open(file_1000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_1000bp = line7[7]
			fasta_ignored_missing_percent_1000bp= line8[7]
			fasta_ignored_small_percent_1000bp = line9[7]
			fasta_ignored_total_percent_1000bp = line11[0]
			mean_1000bp = line15[3]
			sd_1000bp = line16[5]
			range_low_1000bp = line17[4]
			range_high_1000bp = line17[7]
			median_gc_1000bp = line18[3]
			contigs_analyzed_1000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 1000bp')
	if len(file_3000bp)==0:
		print('No log file for', key, 'at 3000bp. Skipping...')
		pass
	else:
		with open(file_3000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_3000bp = line7[7]
			fasta_ignored_missing_percent_3000bp= line8[7]
			fasta_ignored_small_percent_3000bp = line9[7]
			fasta_ignored_total_percent_3000bp = line11[0]
			mean_3000bp = line15[3]
			sd_3000bp = line16[5]
			range_low_3000bp = line17[4]
			range_high_3000bp = line17[7]
			median_gc_3000bp = line18[3]
			contigs_analyzed_3000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 3000bp')
	if len(file_5000bp)==0:
		print('No log file for', key, 'at 5000bp. Skipping...')
		pass
	else:
		with open(file_5000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_5000bp = line7[7]
			fasta_ignored_missing_percent_5000bp= line8[7]
			fasta_ignored_small_percent_5000bp = line9[7]
			fasta_ignored_total_percent_5000bp = line11[0]
			mean_5000bp = line15[3]
			sd_5000bp = line16[5]
			range_low_5000bp = line17[4]
			range_high_5000bp = line17[7]
			median_gc_5000bp = line18[3]
			contigs_analyzed_5000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 5000bp')
	if len(file_20000bp)==0:
		print('No log file for', key, 'at 20000bp. Skipping...')
		pass
	else:
		with open(file_20000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_20000bp = line7[7]
			fasta_ignored_missing_percent_20000bp= line8[7]
			fasta_ignored_small_percent_20000bp = line9[7]
			fasta_ignored_total_percent_20000bp = line11[0]
			mean_20000bp = line15[3]
			sd_20000bp = line16[5]
			range_low_20000bp = line17[4]
			range_high_20000bp = line17[7]
			median_gc_20000bp = line18[3]
			contigs_analyzed_20000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 20000bp')
	if len(file_80000bp)==0:
		print('No log file for', key, 'at 80000bp. Skipping...')
		pass
	else:
		with open(file_80000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_80000bp = line7[7]
			fasta_ignored_missing_percent_80000bp= line8[7]
			fasta_ignored_small_percent_80000bp = line9[7]
			fasta_ignored_total_percent_80000bp = line11[0]
			mean_80000bp = line15[3]
			sd_80000bp = line16[5]
			range_low_80000bp = line17[4]
			range_high_80000bp = line17[7]
			median_gc_80000bp = line18[3]
			contigs_analyzed_80000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 80000bp')
	if len(file_320000bp)==0:
		print('No log file for', key, 'at 320000bp. Skipping...')
		pass
	else:
		with open(file_320000bp[0], 'r') as f:
			content = [x.strip() for x in f.readlines()] #from http://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list-with-python#comment41651423_3277516
			if content[11].startswith("Warning"):  #Checks to see if there's a warning for having a lot of contigs not included which adds a line to the log file.
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split(), content[18].split()
			else:
				line6, line7, line8, line9, line11, line15, line16, line17, line18, line19 = content[5].split(), content[6].split(), content[7].split(), content[8].split(), content[10].split(), content[13].split(), content[14].split(), content[15].split(), content[16].split(), content[17].split()
			fasta_number = line6[0] #only need to define once, but in each call since a file may be missing
			fasta_recovered_percent_320000bp = line7[7]
			fasta_ignored_missing_percent_320000bp= line8[7]
			fasta_ignored_small_percent_320000bp = line9[7]
			fasta_ignored_total_percent_320000bp = line11[0]
			mean_320000bp = line15[3]
			sd_320000bp = line16[5]
			range_low_320000bp = line17[4]
			range_high_320000bp = line17[7]
			median_gc_320000bp = line18[3]
			contigs_analyzed_320000bp = line19[3]
		f.close()
		print('Wrote data for', key, 'at 320000bp')
	#Outputs all of these variables to a new line in the output file:
	outfile.write(Sample + '\t' + Class + '\t' + fasta_number + '\t' + fasta_recovered_percent_1000bp + '\t' + fasta_ignored_missing_percent_1000bp + '\t' + fasta_ignored_small_percent_1000bp + '\t' + fasta_ignored_total_percent_1000bp + '\t' + mean_1000bp + '\t' + sd_1000bp + '\t' + range_low_1000bp + '\t' + range_high_1000bp + '\t' + median_gc_1000bp + '\t' + contigs_analyzed_1000bp + '\t' + fasta_recovered_percent_3000bp + '\t' + fasta_ignored_missing_percent_3000bp + '\t' + fasta_ignored_small_percent_3000bp + '\t' + fasta_ignored_total_percent_3000bp + '\t' + mean_3000bp + '\t' + sd_3000bp + '\t' + range_low_3000bp + '\t' + range_high_3000bp + '\t' + median_gc_3000bp + '\t' + contigs_analyzed_3000bp + '\t' + fasta_recovered_percent_5000bp + '\t' + fasta_ignored_missing_percent_5000bp + '\t' + fasta_ignored_small_percent_5000bp + '\t' + fasta_ignored_total_percent_5000bp + '\t' + mean_5000bp + '\t' + sd_5000bp + '\t' + range_low_5000bp + '\t' + range_high_5000bp + '\t' + median_gc_5000bp + '\t' + contigs_analyzed_5000bp + '\t' + fasta_recovered_percent_20000bp + '\t' + fasta_ignored_missing_percent_20000bp + '\t' + fasta_ignored_small_percent_20000bp + '\t' + fasta_ignored_total_percent_20000bp + '\t' + mean_20000bp + '\t' + sd_20000bp + '\t' + range_low_20000bp + '\t' + range_high_20000bp + '\t' + median_gc_20000bp + '\t' + contigs_analyzed_20000bp + '\t' + fasta_recovered_percent_80000bp + '\t' + fasta_ignored_missing_percent_80000bp + '\t' + fasta_ignored_small_percent_80000bp + '\t' + fasta_ignored_total_percent_80000bp + '\t' + mean_80000bp + '\t' + sd_80000bp + '\t' + range_low_80000bp + '\t' + range_high_80000bp + '\t' + median_gc_80000bp + '\t' + contigs_analyzed_80000bp + '\t' + fasta_recovered_percent_320000bp + '\t' + fasta_ignored_missing_percent_320000bp + '\t' + fasta_ignored_small_percent_320000bp + '\t' + fasta_ignored_total_percent_320000bp + '\t' + mean_320000bp + '\t' + sd_320000bp + '\t' + range_low_320000bp + '\t' + range_high_320000bp + '\t' + median_gc_320000bp + '\t' + contigs_analyzed_320000bp + '\n')
	#Outputs the same information in a format better for feeding into the ggplot2 package in R:
	outfile_ggplot.write(Sample + '\t' + Class + '\t' + fasta_number + '\t' + '1000' + '\t' + fasta_recovered_percent_1000bp + '\t' + fasta_ignored_missing_percent_1000bp + '\t' + fasta_ignored_small_percent_1000bp + '\t' + fasta_ignored_total_percent_1000bp + '\t' + mean_1000bp + '\t' + sd_1000bp + '\t' + range_low_1000bp + '\t' + range_high_1000bp + '\t' + median_gc_1000bp + '\t' + contigs_analyzed_1000bp + '\n' + Sample + '\t' + Class + '\t' + fasta_number + '\t' + '3000' + '\t' + fasta_recovered_percent_3000bp + '\t' + fasta_ignored_missing_percent_3000bp + '\t' + fasta_ignored_small_percent_3000bp + '\t' + fasta_ignored_total_percent_3000bp + '\t' + mean_3000bp + '\t' + sd_3000bp + '\t' + range_low_3000bp + '\t' + range_high_3000bp + '\t' + median_gc_3000bp + '\t' + contigs_analyzed_3000bp + '\n' + Sample + '\t' + Class + '\t' + fasta_number + '\t' + '5000' + '\t' + fasta_recovered_percent_5000bp + '\t' + fasta_ignored_missing_percent_5000bp + '\t' + fasta_ignored_small_percent_5000bp + '\t' + fasta_ignored_total_percent_5000bp + '\t' + mean_5000bp + '\t' + sd_5000bp + '\t' + range_low_5000bp + '\t' + range_high_5000bp + '\t' + median_gc_5000bp + '\t' + contigs_analyzed_5000bp + '\n' + Sample + '\t' + Class + '\t' + fasta_number + '\t' + '20000' + '\t' + fasta_recovered_percent_20000bp + '\t' + fasta_ignored_missing_percent_20000bp + '\t' + fasta_ignored_small_percent_20000bp + '\t' + fasta_ignored_total_percent_20000bp + '\t' + mean_20000bp + '\t' + sd_20000bp + '\t' + range_low_20000bp + '\t' + range_high_20000bp + '\t' + median_gc_20000bp + '\t' + contigs_analyzed_20000bp + '\n' + Sample + '\t' + Class + '\t' + fasta_number + '\t' + '80000' + '\t' + fasta_recovered_percent_80000bp + '\t' + fasta_ignored_missing_percent_80000bp + '\t' + fasta_ignored_small_percent_80000bp + '\t' + fasta_ignored_total_percent_80000bp + '\t' + mean_80000bp + '\t' + sd_80000bp + '\t' + range_low_80000bp + '\t' + range_high_80000bp + '\t' + median_gc_80000bp + '\t' + contigs_analyzed_80000bp + '\n' + Sample + '\t' + Class + '\t' + fasta_number + '\t' + '320000' + '\t' + fasta_recovered_percent_320000bp + '\t' + fasta_ignored_missing_percent_320000bp + '\t' + fasta_ignored_small_percent_320000bp + '\t' + fasta_ignored_total_percent_320000bp + '\t' + mean_320000bp + '\t' + sd_320000bp + '\t' + range_low_320000bp + '\t' + range_high_320000bp + '\t' + median_gc_320000bp + '\t' + contigs_analyzed_320000bp + '\n')
	
outfile.close()
outfile_ggplot.close()
print('Regular output:', outfile_path)
print('ggplot2 output:', outfile_ggplot_path)

##  Record the endtime time:
then = datetime.datetime.now()
print('Time at the end of the script:', then)
endtime = time.clock()
elapsedtime = endtime - starttime
print(elapsedtime, "seconds elapsed running this script (i.e., ", elapsedtime/60, "minutes).\n\n======================================\n")
