#!/usr/bin/env python 

print("start")

import pysam
import re
import sys
BAMFILE=sys.argv[1]
samfile = pysam.AlignmentFile(BAMFILE,"rb",ignore_truncation=True)
output = pysam.AlignmentFile(BAMFILE+"_parsed.bam","wb",template=samfile)

# read each line
for read in samfile:
	query = str(re.findall(":[:A-Z0-9]+",read.query_name)[0])[1:]
	barcode = str(re.findall("^[A-Z]+",read.query_name)[0])
	read_tmp = read
	read_tmp.tags += [("CB",barcode)]
	read_tmp.query_name=query
	output.write(read_tmp)

samfile.close()
output.close()
print("parsing finished")
