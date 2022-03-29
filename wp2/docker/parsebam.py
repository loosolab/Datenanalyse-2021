#!/usr/bin/env python 

print("start")

import pysam
import re
import sys
# Bamfile to parse
BAMFILE=sys.argv[1]
# read Bamfile
samfile = pysam.AlignmentFile(BAMFILE+".bam","rb",ignore_truncation=True)
# write Bamfile
output = pysam.AlignmentFile(BAMFILE+"_parsed.bam","wb",template=samfile)

# iterate over lines in Bamfile
for read in samfile:
		# save query without barcode
		query = str(re.findall(":[:A-Z0-9]+",read.query_name)[0])[1:]
		# read barcode from query
		barcode = str(re.findall("^[A-Z]+",read.query_name)[0])
		read_tmp = read
		# add barcode to tags
		read_tmp.tags += [("CB",barcode)]
		# override query
		read_tmp.query_name=query
		# write line
		output.write(read_tmp)

samfile.close()
output.close()
