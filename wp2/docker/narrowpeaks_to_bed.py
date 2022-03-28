#!/usr/bin/env python 

import os
import re
import sys

files = os.listdir(sys.argv[1])
print(files)
for file in files:
        if re.match(".*narrowPeak$",file):
            new_file_name = ".".join(file.split(".")[:-1])+".bed"
            with open(sys.argv[1]+new_file_name, 'w') as outfile, open(sys.argv[1]+file, 'r') as infile:
                for line in infile:
                    start = line.split()[1]
                    stop = line.split()[2]
                    chr = re.findall("chr[0-9]*",line.split()[0])[0]
                    if chr == "chr" or chr == "chrY" or chr == "chrX" or chr ==  "chrM":
			continue
                    newline = chr+"\t"+start+"\t"+stop+"\n"
                    outfile.write(newline)
