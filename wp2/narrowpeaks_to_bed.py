#!/usr/bin/env python 

import os
import re
import sys

files = os.listdir(sys.argv[1])
# print all files in folder
print(files)

for file in files:
	# parse narrowPeak files in folder
        if re.match(".*narrowPeak$",file):
	    #create .bed file
            new_file_name = ".".join(file.split(".")[:-1])+".bed"
            with open(sys.argv[1]+new_file_name, 'w') as outfile, open(sys.argv[1]+file, 'r') as infile:
		# only add chr start and stop
                for line in infile:
                    start = line.split()[1]
                    stop = line.split()[2]
                    chr = re.findall("chr[0-9]*",line.split()[0])[0]
		    # only add autosomes
                    if chr == "chr" or chr == "chrY" or chr == "chrX" or chr ==  "chrM":
			continue
                    newline = chr+"\t"+start+"\t"+stop+"\n"
                    outfile.write(newline)
