# isochores
Code used to analyze GC content in genomes downloaded from GenBank.

This project consists of three scripts:
  1. deinterleave.py
  2. isochore_analyzer1.831.py
  3. isochore_summarizer_0.56.py

1. deinterleave.py - After downloading DNA FASTA-formatted files, this is the first script to use. It removes newline characters within a sequence, thus allowing the other two scripts to run properly.
2. isochore_analyzer1.831.py - This script analyzes individual FASTA-formatted files for GC content in a specified window size. You must also provide an organism name and group name. If no window size is provided, the default is 5000 base pairs. This script creates several outputs, all placed within one folder. View the commented text in the script to understand the format of these outputs.
3. isochore_summarizer_0.56.py - This script assumes the isochore_analyzer script was already run. It scans the log files created by the isochore_analyzer script and summarizes their contents in a tab-separated text file with a header row. The commented text within this script describes the exact file hierarchy and log file format necessary for it to work appropriately.

These scripts were originally written for personal use and largely provided as-is. If you have questions about getting the scripts to work properly, contact me via email: allopatry@gmail.com

-Alex
