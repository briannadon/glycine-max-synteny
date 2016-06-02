from te import *
import re
from sys import argv

script, align, bed = argv

with open(bed) as b:
	gc = {}
	for line in b.readlines():
		line = line.split()
		gc[line[3]] = [line[0],int(line[1]),int(line[2])]

with open(align) as a:
	alignment = {}
	for line in a.readlines():
		add_align(alignment,line)

print "Alignment\tside\tChr\tStart\tEnd"

for key in sorted(alignment.keys()):
	aargs = (alignment, gc) + key
	s = get_start(*aargs)
	e = get_end(*aargs)
	ch = gc[alignment[key][0]][0]
	print('\t'.join(map(str,key)) + "\t%s\t%d\t%d" % (ch, s, e))
