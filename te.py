from sys import argv
import re

def te_length(masked, chrom, start, end): ## 'te' is a dict of lines of a number:[ChrXX,Start,End]
	block = masked[chrom][(start-1):end]
	length = block.count('N')
	return length

def get_chrom(alignment, genecoords, x, side):
	chrom = genecoords[alignment[(x, side)][0]][0]
	return chrom

def get_start(alignment, genecoords, x, side):
	start = min(genecoords[alignment[(x,side)][0]][1:3] + genecoords[alignment[(x,side)][-1]][1:3])
	return start

def get_end(alignment, genecoords, x, side):
	end = max(genecoords[alignment[(x,side)][0]][1:3] + genecoords[alignment[(x,side)][-1]][1:3])
	return end

def add_align(alignment, line):
	line = line.strip().split()
	if ((int(line[0]), 'L') in alignment or (int(line[0]), 'R') in alignment):
		alignment[(int(line[0]), 'L')].append(line[2])
		alignment[(int(line[0]), 'R')].append(line[3])
	else:
		alignment[(int(line[0]), 'L')] = [line[2]]
		alignment[(int(line[0]), 'R')] = [line[3]]

def main():

	script, aligns, mask, bed, ages = argv

	with open(bed, 'r') as b:
		genecoords = {} ##genecoords is a dict of GeneName:[Chr,start,end] 
		blines = b.readlines()
		for line in blines:
			line = line.split()
			genecoords[line[3]] = [line[0], int(line[1]), int(line[2])]

	with open(ages, 'r') as a:
		alines = a.readlines()
		al = {}
		for line in alines:
			al[int(line.split()[0])] = line.split()[2]

	with open(mask) as m:
		c = re.compile('Chr\d\d')
		masked = {}
		lines = m.readlines()
		for line in lines:
			if ">Chr" in line:
				curchrom = c.search(line.strip()).group()
				masked[curchrom] = []
			elif "scaffold" in line:
				break
			else:
				seq = line.strip()
				masked[curchrom] += list(seq)

	with open(aligns, 'r') as f:
		lines = f.readlines()
		alignment = {} ## alignment is a dict of {(id, L/R) : [Gene,gene,gene...] } 
		for line in lines:
			add_align(alignment, line)
		print "ID\tTE Length\tAlignment Length\tTE Proportion\tAge\tSide"
		for i in xrange(0,540):
			for side in ('L','R'):
				aargs = (alignment, genecoords, i, side)
				astart = get_start(*aargs)
				aend = get_end(*aargs)
				tel = te_length(masked, get_chrom(*aargs), astart, aend)
				alen = aend - astart
				prop = float(tel) / float(alen)
				print "%d\t%d\t%d\t%f\t%s\t%s" % (i,tel,alen,prop, al[i], side)

if __name__ == "__main__":
	main()

