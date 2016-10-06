#!/usr/bin/python

##  Isochore analysis pipeline. The goal is to generate a histogram with the abundance of GC as a percent of 5,000 base pair fragments (or fragments of another size).
##  A script by Alex Hall
##  allopatry@gmail.com
##  June 8, 2016
##  Version 1.831


##  Usage:
##    >>> isochore_analyzer.py filename.fasta [organism] [group] [5000]
##  Reads input file 'filename.fasta' and uses the second value (if provided) as the non-overlapping window size to analyze (default is 5000). The third value is interpreted as a string and can be the species name of the genome being studied. This helps keep output separate.
##  Output are three files in a folder created in the working directory called output_ followed by the system time at the start of the analysis, to prevent overwriting old output:
##    gc_content_timesinceepoch.tsv << this is the main output as a list of lines where each line has a contig and number of the segment (by window size), a tab character, and the GC percentage out of 100
##    unused_contigs_missing_data_timesinceepoch.tsv << These were unused data as a result of having too many ambiguous or missing sequence. Each line has a contig and number of the segment (by window size), a tab character, and the number of times that exact key was discarded for this reason.
##    unused_short_contigs_content_timesinceepoch.tsv << These were unused data as a result of the contig size being smaller than the provided window size (defaults to 5000). Each line has a contig and number of the segment (by window size), a tab character, and the number of times that exact key was discarded for this reason.

##  Setup environment:
from  __future__ import print_function
import sys, re, os, time, datetime, errno, math, numpy
from collections import OrderedDict
scriptfile, inputfile, path = sys.argv[0], sys.argv[1], os.getcwd()
try:
	organism = str(sys.argv[2])
except IndexError:
	organism = ''
try:
	group = str(sys.argv[3])
except IndexError:
	group = ''
try:
	windowsize_commandline = int(sys.argv[4])
except IndexError:
	windowsize_commandline = 5000

##  Record the start time:
starttime = time.clock()
timesinceepoch = int(time.time())
now = datetime.datetime.now()
print('\n======================================\n\nThe current time is:', now)

##Set up a log file
path_starttime= path+'/'+str(timesinceepoch)
outfolder='Run1.831_GC_output_'+str(organism)+'_'+str(windowsize_commandline)+'_'+str(timesinceepoch)
outfolder_path=path+'/'+outfolder
try:
	os.makedirs(outfolder)
except OSError as exc:  # Python >2.5
	if exc.errno == errno.EEXIST and os.path.isdir(path):
		pass
	else:
		raise
##Check that a valid window size was provided
if windowsize_commandline < 1:
	sys.exit('Error: the window size you provided was less than 1 base.')
print('Window size provided or assumed:', windowsize_commandline, '(i.e.,',windowsize_commandline/1000,'kb).')

##  Check that input is in fasta format:
inputdata = open(inputfile, 'r')
for line in inputdata:
	if line.startswith('#'):
		continue
	else:
		if line.startswith('>'):
			print('Looks like you supplied a fasta-formatted file.')
			break
		else:
			sys.exit('Error: could not verify that the input file was a fasta-formatted file. Check the input data.')

##  Check that file is not too large for the computer to handle:
print('Checking that you have enough physical memory to analyze this file.')
combined = path+'/'+inputfile
filesize = os.path.getsize(combined)
physicalmem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') ##this only works on Linux systems, not on OSX
print('Physical memory =', physicalmem/(1024**2), 'Mb.')
print('File size =', filesize/(1024**2), 'Mb.')
if filesize>physicalmem:
	sys.exit('Error: file size exceeded physical memory.')
else:
	print('The computer should have enough memory to continue.')

##  Setup variables:
print('Setting up variables, initializing dictionaries, and creating output files.')
fragment_length, missingno, GCno, ATno, Ano, Tno, Wno, Gno, Cno, Sno, GCPercent, endtime, headernum, dtotal, GCmean, GCsd, GCrange_lower, GCrange_upper, GCmedian, GCnum, contiglength, contignum_round, headerposition, headerprogress, totaltime = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

dgc, dmissing, dshort, ddata, dgc_temp = {}, {}, {}, {}, {}
sd_list = []

##  Create output files and folder:
dgc_ = 'gc_content_'+str(organism)+'_'+str(windowsize_commandline)+'_'+str(timesinceepoch)+'.tsv'
dmissing_ = 'unused_contigs_missing_data_'+str(organism)+'_'+str(windowsize_commandline)+'_'+str(timesinceepoch)+'.tsv'
dshort_ = 'unused_short_contigs_content_'+str(organism)+'_'+str(windowsize_commandline)+'_'+str(timesinceepoch)+'.tsv'
dgc_output = open(os.path.join(outfolder_path, dgc_), 'w')
dmissing_output = open(os.path.join(outfolder_path, dmissing_), 'w')
dshort_output = open(os.path.join(outfolder_path, dshort_), 'w')
log_ = 'log_'+str(organism)+'_'+str(windowsize_commandline)+'_'+str(timesinceepoch)+'.txt'
log_output = open(os.path.join(outfolder_path, log_), 'a')

##Initialize an input data dictionary [ddata]:
inputdata.seek(0)
i_inputdata = iter(inputdata)
for line in i_inputdata:
	if line.startswith('>'):
		headernum = headernum+1
		header= line.replace(">","")
		header2= header.strip()
		header3= header2.replace(" ","_")
##		Insert code here later to deal with interleaved fasta files
		ddata[header3] = next(i_inputdata).strip()
print(headernum,'fasta headers (i.e., contigs) were detected.')

##  Populate three dictionaries with the ddata dictionary and window size equal to windowsize
print('\nBeginning GC content analysis using a non-overlapping window of', windowsize_commandline, 'base pairs (i.e.,', windowsize_commandline/1000, 'kb or', windowsize_commandline/1000000, 'Mb).')
for key in sorted(ddata.iterkeys()):	##in python 3.X use dgc_temp.keys()
	dgc_temp = {}
	contig_temp=ddata[key]
	start = 0
	end = windowsize_commandline-1
	contiglength= len(contig_temp)
	if contiglength >= windowsize_commandline == False:
		if not key in dshort.keys():
			dshort[key] = 1
		else:
			dshort[key] = dshort[key]+1
	contignum_round=int(round((contiglength/windowsize_commandline)))
	for i in range(0,contignum_round+1):
		key_name=str(key)+'_'+str(i)
		dgc_temp[key_name] = ''
	for key in sorted(dgc_temp.iterkeys()):	##in python 3.X use dgc_temp.keys()
		contig_temp2 = contig_temp[start:end+1]
		dgc_temp[key] = contig_temp2
		start = start+windowsize_commandline
		if ((len(contig_temp)-start)>windowsize_commandline)==True:
			end = end+windowsize_commandline
		else:
			end = (len(contig_temp)-1)
	for key in sorted(dgc_temp.iterkeys()): ##in python 3.X use dgc_temp.keys()
		contig_temp3=dgc_temp[key]
		fragment_length=len(contig_temp3)
		fragment_length_float=float(len(contig_temp3))
		if not fragment_length == windowsize_commandline:
			if not key in dshort.keys():
				dshort[key] = 1
				break
			else:
				dshort[key] = dshort[key]+1
				break
		hyphenno = len(re.findall('\-', contig_temp3, re.IGNORECASE))
		questionmarkno = len(re.findall('\?', contig_temp3, re.IGNORECASE))
		Ano = len(re.findall('A', contig_temp3, re.IGNORECASE))
		Cno = len(re.findall('C', contig_temp3, re.IGNORECASE))
		Tno = len(re.findall('T', contig_temp3, re.IGNORECASE))
		Gno = len(re.findall('G', contig_temp3, re.IGNORECASE))
		Wno = len(re.findall('W', contig_temp3, re.IGNORECASE))
		Sno = len(re.findall('S', contig_temp3, re.IGNORECASE))
		ATno = Ano+Tno+Wno
		GCno = Gno+Cno+Sno
		if (fragment_length_float-(ATno+GCno))/fragment_length_float >= 0.2:
			if not key in dmissing.keys():
				dmissing[key] = 1
				break
			else:
				dmissing[key] = dmissing[key]+1
				break
		else:
			GCPercent = float(GCno/float((ATno+GCno)))*100
			if not key in dgc.keys():
				dgc[key] = GCPercent
			else:
				keyd=key+'_duplicate'
				dgc[keyd] = GCPercent

##  Count number of missing values as percentage and report to terminal:
dgc_num = float(len(dgc.keys()))
dmissing_num = float(len(dmissing.keys()))
dshort_num = float(len(dshort.keys()))
dtotal_num = float(dgc_num+dmissing_num+dshort_num)
dgc_perc = float((dgc_num/dtotal_num)*100)
dmissing_perc = float((dmissing_num/dtotal_num)*100)
dshort_perc = float((dshort_num/dtotal_num)*100)
print('\nAs a percent of all contigs recovered,',dgc_perc,'percent could have their GC content analyzed.')
print('As a percent of all contigs recovered,',dmissing_perc,'percent could not have their GC content analyzed because there was too much ambiguous or missing data.')
print('As a percent of all contigs recovered,',dshort_perc,'percent could not have their GC content analyzed because the contig was not larger than or equal to',windowsize_commandline,'base pairs.')
if (dmissing_perc+dshort_perc) >= 10:
	print(dmissing_perc+dshort_perc,'percent of contigs could not be used to calculate GC content.'	)
	print('Warning! The amount of missing data exceeded 10 percent as measured by the percentage of contigs that could be used. Be sure this is okay for the interpretation of your results.')
else:
	print(dmissing_perc+dshort_perc,'percent of contigs could not be used to calculate GC content.')

##  Summarize the results:
GCnum = len(dgc.keys())
GCsort = sorted(dgc.values())
if GCnum > 0.0:
	GCmean = float(sum(GCsort)/GCnum)
	GCsd = numpy.std(GCsort) #you can replace the numpy call with your own standard deviation equation
	GCrange_lower = GCsort[0]
	GCrange_upper = GCsort[GCnum-1]
	if not GCnum % 2:
		GCmedian= float((GCsort[GCnum/2]+GCsort[(GCnum/2)-1])/2.0)
	GCmedian= float(GCsort[GCnum/2])
else:
	GCmean, GCsd, GCrange_lower, GCrange_upper, GCmedian = 0, 0, 0, 0, 0
print('\nSummary of results from isochore analysis by non-overlapping',windowsize_commandline,'base pair windows:')
print('Mean GC content:',GCmean,'percent.')
print('Standard deviation of GC content:',GCsd, 'percent.')
print('GC content ranges from',GCrange_lower,'percent to',GCrange_upper,'percent.')
print('Median GC content:',GCmedian,'percent.')
print('Overall, there were',GCnum,'contigs analyzed.')

##  Populate output: tsv of dgc, dmissing, dshort
print('\nWriting output files to', outfolder_path)
with dgc_output as d:
	for key, value in dgc.items():
		d.write("{}\t{}\t{}\t{}\n".format(key, value, organism, group))
print('Successfully created', dgc_)
with dmissing_output as m:
	for key, value in dmissing.items():
		m.write("{}\t{}\t{}\t{}\n".format(key, value, organism, group))
print('Successfully created', dmissing_)
with dshort_output as s:
	for key, value in dshort.items():
		s.write("{}\t{}\t{}\t{}\n".format(key, value, organism, group))
print('Successfully created', dshort_)

##  Record the endtime time:
then = datetime.datetime.now()
print('The current time is:', then)
endtime = time.clock()
elapsedtime = endtime - starttime
print(elapsedtime, "seconds elapsed running this script (i.e., ", elapsedtime/60, "minutes).")
print('Now writing a log file to:', log_)

with log_output as l:
	print('Script started at:', now, '\nThe file analyzed was:', combined, '\nThe file was a valid fasta file format.', '\nPassed physical memory check.', '\nWindow size provided or assumed:', windowsize_commandline, '(i.e.,', windowsize_commandline/1000, 'kb).\n',headernum, 'fasta headers (i.e., contigs) were detected.\n','As a percent of all contigs recovered,' ,dgc_perc, 'percent could have their GC content analyzed.\n', 'As a percent of all contigs recovered,', dmissing_perc, 'percent could not have their GC content analyzed because there was too much ambiguous or missing data.', '\nAs a percent of all contigs recovered,' ,dshort_perc, 'percent could not have their GC content analyzed because the contig was not larger than or equal to', windowsize_commandline, 'base pairs.\n', file=l)
	if (dmissing_perc+dshort_perc) >= 10:
		print(dmissing_perc+dshort_perc, 'percent of contigs could not be used to calculate GC content.\n', 'Warning! The amount of missing data exceeded 10 percent as measured by the percentage of contigs that could be used. Be sure this is okay for the interpretation of your results.\n', file=l)
	else:
		print(float(dmissing_perc+dshort_perc), 'percent of contigs could not be used to calculate GC content.\n', file=l)
	print('Summary of results from isochore analysis by non-overlapping', windowsize_commandline, 'base pair windows:\n', 'Mean GC content:', GCmean, 'percent.\n', 'Standard deviation of GC content:', GCsd, 'percent.\n', 'GC content ranges from', GCrange_lower, 'percent to', GCrange_upper, 'percent.\n', 'Median GC content:', GCmedian, 'percent.\n', 'Overall, there were', GCnum, 'contigs analyzed.\n', 'Script finished at:', then,'\n', elapsedtime, 'seconds elapsed running this script (i.e.,', round(elapsedtime/60, 2), 'minutes).', file=l)
print('Successfully created', log_)


##  Close all files.
dgc_output.close()
dmissing_output.close()
dshort_output.close()
log_output.close()
sys.exit('All finished!')
