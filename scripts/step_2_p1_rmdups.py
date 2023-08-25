from Bio import SeqIO
import sys
import os

input_file=sys.argv[1]
path=sys.argv[2]
outname="temp_out"
os.chdir(path)

seen = []
records = []

for record in SeqIO.parse(input_file, "fasta"):  
  if str(record.seq) not in seen:
    seen.append(str(record.seq))
    records.append(record)

#writing to a fasta file
SeqIO.write(records, outname, "fasta")
