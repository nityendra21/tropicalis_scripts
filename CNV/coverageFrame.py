#!/usr/bin/python
import re, math, sys, errno, fileinput

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-t", "--louchetreshold", dest="loucheTreshold", type="int", default=500)
parser.add_option("-f", "--framesize", dest="frameSize", type="int", default=1000)
parser.add_option("-s", "--framestep", dest="frameStep", type="int", default=100)
parser.add_option("-l", "--readlength", dest="readLength", type="int", default=300)
parser.add_option("-p", "--printUnaligned", dest="printUnaligned", action="store_true", default=False)
(options, args) = parser.parse_args()
loucheTreshold 	= options.loucheTreshold
frameSize 		= options.frameSize
frameStep 		= options.frameStep
readLength 		= options.readLength
printUnaligned 	= options.printUnaligned
# loucheTreshold 	= 500
# frameSize		= 1000
# frameStep		= 100
# readLength 	= 200


unalignedReads 	= 0
alignedReads	= 0

def extendedCeil(x, m):  # Note: not actually ceil, because if x is a multiple of m, it will return x+m. This is intended behavior
	return x + m - x % m # Actual ceil implementation: x if x % m == 0 else x + m - x % m
def extendedFloor(x, m):
	return x - (x % m)

# Class representing reads as read from the sam file
class SamRead:
	def __init__(self, line):
		dc,dc,scaf,pos,mapq,cigar,rnext,pnext,dc = line.split('\t',8)
		self.rname 		= scaf
		self.pos		= int(pos)-1 # Convert to 0-based system
		# self.cigar		= SamRead.expandCigar(cigar)
		# self.refcigar	= re.sub("I|S|H|P|N","",self.cigar)
		# self.refcigar = self.cigar.replace('I','').replace('S','').replace('H','').replace('P','').replace('N','')
		self.refcigar 	= SamRead.expandRefCigar(cigar)
		self.mapq		= int(mapq)
		self.length 	= len(self.refcigar)
		self.isLouche 	= (rnext=='=' or rnext==self.rname) and (abs(self.pos-int(pnext))<=loucheTreshold)
		self.isAligned 	= scaf != '*'

	@staticmethod
	def expandCigar(cigar):
		if (not cigar) or cigar == "*":
			return ""
		regex = re.match(r'(\d*)([MIDNSHPX=])(.*)',cigar)
		return regex.group(2)*int(regex.group(1)) + SamRead.expandCigar(regex.group(3))

	@staticmethod # Expand the cigar string but only those operations that are present on the reference string
	def expandRefCigar(cigar):
		rep=0
		i=0
		while i<len(cigar):
			char = cigar[i]
			if char.isdigit():
				rep = rep*10 + int(char)
			else:
				if char in ['M','X','=']:
					return char*rep + SamRead.expandRefCigar(cigar[i+1:])
				else:
					return SamRead.expandRefCigar(cigar[i+1:])
			i=i+1
		# for char in cigar:
		# 	if char.isdigit():
		# 		rep = rep*10 + int(char)
		# 	else:
		# 		return char*rep
		return ""
		

	# Count the number of matches of this read in a specific region on the reference
	def countMatches(self, start, end):
		leftCoord 	= max(self.pos, start)
		rightCoord 	= min(self.pos+self.length, end)
		if rightCoord<leftCoord:
			return 0
		else:
			## Previous implementation: enormously inefficient
			# return len(re.findall("M|=|X",self.refcigar)[start-self.pos:end-self.pos])
			return len((self.refcigar[start-self.pos:end-self.pos]).replace('D',''))

# Class containing read frames. Also handles output and global iteration, i.e. keeps the frames in the vicinity and current position
class Frame:
	currentScaf = 0				# Index of the currently processing scaffold
	currentPos	= 0				# Position of the last processed read
	current 	= set()			# Set containing all the frames in the vicinity of the position of the global iterator
	positions	= set()			# The positions of the sets in the current vicinity

	def __init__(self, pos):
		self.pos 		= pos	# Starting position of the frame
		self.end 		= min(pos+frameSize, scaffoldSizes[Frame.currentScaf])
		self.depth		= 0		# Total read depth in the frame
		self.louches 	= 0		# Total read depth caused by "louche" reads
		self.mapq 		= 0		# Total mapq score of all reads aligned in this frame
		self.amtReads	= 0		# Total amount of reads that cover at least one base in this frame
		Frame.current.add(self)
		Frame.positions.add(pos)

	def addCoverage(self, read):
		addedCov 	 = read.countMatches(self.pos, self.pos+frameSize)
		self.depth    += addedCov
		self.louches  += addedCov*read.isLouche
		self.mapq	  += read.mapq
		self.amtReads += (addedCov > 0)

	def flush(self):
		try:
			print scaffolds[Frame.currentScaf],"\t",self.pos,"\t",self.end,"\t",self.depth,"\t",self.louches,"\t",self.mapq/max(1, self.amtReads),"\t",self.amtReads
		except IOError as e:
			if e.errno == errno.EPIPE: sys.exit()
		# print ', {',self.pos,',',self.depth,',',self.louches,',',self.mapq/max(1,self.amtReads),',',self.amtReads,'}',
		Frame.current.discard(self)
		Frame.positions.discard(self.pos)
		

def updatePosition(read):
	pos 			= read.pos
	prevpos 		= max(0, extendedCeil(Frame.currentPos-frameSize, frameStep)) if read.rname == scaffolds[Frame.currentScaf] else 0
	Frame.currentPos= pos
	criticalStart 	= max(0, extendedCeil(pos-frameSize, frameStep))
	criticalEnd		= min(extendedCeil(pos+readLength, frameStep), scaffoldSizes[Frame.currentScaf])
	inRange			= set(range(min(criticalStart, prevpos), criticalEnd+1, frameStep))
	# Create new frames
	newframes		= inRange - Frame.positions
	for n in newframes:
		Frame(n).addCoverage(read)
	# Flush old frames
	for f in sorted(Frame.current, key=lambda frame: frame.pos):
		if f.pos < criticalStart or f.pos > criticalEnd:
			f.flush()

def processRead(read):
	global unalignedReads, alignedReads
	alignedReads += read.isAligned
	unalignedReads += 1 - read.isAligned
	# Progress to next scaffold if necessary
	if read.rname != scaffolds[Frame.currentScaf]:
		if not read.isAligned: return # Stop when arrived at unaligned reads
		newScaf = scaffolds.index(read.rname)
		if newScaf < Frame.currentScaf: sys.exit("input needs to be an ordered sam file") # Stop when the sam file was unordered
		while Frame.current: Frame.current.pop().flush() # Flush all frames in the current vicinity
		alignedReads 	   += 1
		Frame.currentScaf 	= newScaf
		Frame.currentPos 	= 0
		printScaffold(newScaf)
	# Add the information from this read to the frames
	for f in Frame.current: f.addCoverage(read)
	# Crawl further across scaffold
	# Flush frames that fall out of the vicinity and add new frames if necessary		
	updatePosition(read)



# Ratio: iterate over all reads, add the appropriate amount of bases to each frame in the vicinity.
# Frames in the vicinity include the (frameSize/frameStep) frames which start before the current position (but end after it, due to frame overlam)
# as well as possibly the next (maxReadSize/frameStep) frames if the read spills over its current frame.
# Frames are kept in a list variable and removed once they go out of range, and their values flushed to stdout.

def isHeader(line):
	return line.startswith('@')
def isScaffold(line):
	return line.startswith("@SQ")
def extractScaffoldName(line):
	return re.search('SN:([!-)+-<>-~][!-~]*)',line).group(1)
def extractScaffoldSize(line):
	return int(re.search('LN:(\d*)', line).group(1))
def printScaffold(num):
	pass
	# print "\n",scaffolds[num]

# print SamRead.expandRefCigar('20M3D30X')
# sys.exit()

scaffolds		= []
scaffoldSizes	= []
processHeaders	= True
for line in fileinput.input(args):
	if processHeaders:
		### Extract scaffold information
		if not isHeader(line):
			### Print starting info and header of first scaffold
			processHeaders 	= False;
			print "scaffold\tframeStart\tframeEnd\tdepth\tlouchedepth\tavgmapq\tamtreads"
			# print "FrameSize:"+str(frameSize)+" FrameStep:"+str(frameStep)+" loucheTreshold:"+str(loucheTreshold),
			printScaffold(0)
			continue
		if isScaffold(line):
			scaffolds.append(extractScaffoldName(line))
			scaffoldSizes.append(extractScaffoldSize(line))
	else:
		### Process read info
		currentRead = SamRead(line)
		processRead(currentRead)
# BUG: this method skips the first read!

if printUnaligned:
	print "Aligned:"+str(alignedReads)+" Unaligned:"+str(unalignedReads)
