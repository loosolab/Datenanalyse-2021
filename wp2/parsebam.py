#!/usr/bin/env python 

print("start")

import pysam
import re
import sys
BAMFILE=sys.argv[1]
samfile = pysam.AlignmentFile(BAMFILE+".bam","rb",ignore_truncation=True)
output = pysam.AlignmentFile(BAMFILE+"_parsed.bam","wb",template=samfile)
cnt = 0
print("test")
for read in samfile:
	if cnt == 0:
		#print(read)
		#print("----------------")
		#print(read.query_name)
		#print("CB:"+read.query_name)
		query = str(re.findall(":[:A-Z0-9]+",read.query_name)[0])[1:]
		#print(query)
		barcode = str(re.findall("^[A-Z]+",read.query_name)[0])
		#print(barcode)
		#print("---------------")
		#cnt += 1
		read_tmp = read
		#read_tmp.query_name = 
		#print("query_name:"+read.query_name)
		#read_tmp.query_sequence = 
		#print("query_sequence:"+read.query_sequence)
		#read_tmp.flag = 
		#print("flag:"+str(read.flag))
		#read_tmp.reference_id = 
		#print("reference_id:"+str(read.reference_id))
		#read_tmp.mapping_quality = 
		#print("mapping_quality:"+str(read.mapping_quality))
		#read_tmp.cigar = 
		#print("cigar:"+str(read.cigar))
		#read_tmp.next_reference_id = 
		#print("next_reference_id:"+str(read.next_reference_id))
		#read_tmp.next_reference_start = 
		#print("next_reference_start:"+str(read.next_reference_start))
		#read_tmp.template_lenght = 


		#print("tags:"+str(read_tmp.tags[1]))
		read_tmp.tags += [("CB",barcode)]
		#print("tags:"+str(read_tmp.tags))
		read_tmp.query_name=query
		






		#print("xxx")
		output.write(read_tmp)
		#print("xxx2")

#print("test4")
samfile.close()
output.close()
