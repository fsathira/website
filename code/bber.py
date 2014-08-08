# bsub -q short -W 12:00 -J bber -o bberPy.out -e bberPy.err -R "mem > 30000" python bber.py
# bsub -q short -W 12:00 -J bber -o bberPy.out -e bberPy.err -M 32768 -R "rusage[mem=32768]" python bber.py

# Build dictionary of SNP -> LD ID
#  CHR          BP1          BP2           KB  NSNPS SNPS
LDfile = open("lgrc_eqtl_qc.blocks.det")
#LDfile = open("lgrc_qc.blocks.det")
LDdic = {} # map from snp to LDid
ld2snp = {}
numLD = 0
header = True
for line in LDfile:
	if header:
		header = False
		continue
	aline = line.strip().split()
	snps = aline[5].split("|")
	ld2snp[numLD] = snps
	for snp in snps:
		LDdic[snp] = numLD
	numLD += 1

LDfile.close()
numLD

# Build dictionary of Gene -> Cluster ID
# ENSG00000000003 1
clustdic = {}
clustfile = open("lgrc.expr.hclust.txt")
clust2gene = {}
clust = set([])
for line in clustfile:
	aline = line.strip().split()
	clustid = int(aline[1])
	gene = aline[0]
	clustdic[gene] = clustid
	clust.add(clustid)
	if not clust2gene.has_key(clustid):
		clust2gene[clustid] = []
	clust2gene[clustid].append(gene)

clustfile.close()
numClust = len(clust)
numClust
max(clust) == numClust # should be True

# Build dictionary of dictionary of array, mapping LD -> genes cluster -> array
# Ensembl_Gene    CHR     SNP     BP      A1      TEST    NMISS   BETA    STAT    P
# ENSG00000000003    0   kgp22785549          0    C        ADD      159    -0.1206       -2.103      0.03711
eqtlfile = open("plink.ADD.0.001.linear") # plink.ADD.linear, plink.ADD.0.01.linear
unmapSNP = []
unmapGene = []
eqtldic = {}
header = True
nline = 0
for line in eqtlfile:
	if nline % 1000000 == 0:
		print nline
	nline += 1
	if header:
		header = False
		continue
	aline = line.strip().split()
	thisSNP = aline[2]
	thisGene = aline[0]
	thisP = float(aline[9])
	
	# if SNP is not in the list, it is a singleton. Add it to the list.
	if not LDdic.has_key(thisSNP):
		unmapSNP.append(thisSNP)
		LDdic[thisSNP] = numLD
		ld2snp[numLD] = [thisSNP]
		numLD += 1
	LDid = LDdic[thisSNP]
	
	# if gene is not in the list, they are the duplicated genes. This is a problem.
	if not clustdic.has_key(thisGene):
		unmapGene.append(thisGene)
		continue
	clustid = clustdic[thisGene]
	
	if eqtldic.has_key(LDid):
		if eqtldic[LDid].has_key(clustid): 
			eqtldic[LDid][clustid].append(thisP)
		else:
			eqtldic[LDid][clustid] = [thisP]
	else:
		eqtldic[LDid] = {clustid:[thisP]}

eqtlfile.close()

# number of singleton LD, should be 216109, it might be less here
len(unmapSNP)
# total to a new number of LD
# this is an underestimation as some singleton SNP clusters have insignificant p
numLD

# this should be zero
len(unmapGene)

# BBER
from scipy import stats

def ranks (v):
	# from http://stackoverflow.com/questions/5284646/rank-items-in-an-array-using-python-numpy
	import numpy as np
	t = np.argsort(v)
	r = np.empty(len(v),int)
	r[t] = np.arange(len(v))
	for i in xrange(1, len(r)):
		if v[t[i]] <= v[t[i-1]]: r[t[i]] = r[t[i-1]]
	return r

# Pick 90th percentile as representative, store in an array for sorting
#bberdic = {}
bberP = []
#bberP = [None]*(numLD*numClust)
for LDid in range(1,numLD+1):
	if LDid % 10000 == 0:
		print LDid
	for clustid in range(1,numClust+1):
		# If the block is not there, all P values are >= 0.05. We use 0.05 as representative stat
		# It will all be ignored at the end anyway.
		if not eqtldic.has_key(LDid) or not eqtldic[LDid].has_key(clustid):
			bberP.append(0.05)
			continue
		allP = eqtldic[LDid][clustid]
		totalBlockSize = len(ld2snp[LDid])*len(clust2gene[clustid])
		allP += [1]*(totalBlockSize - len(allP))
		thisstat = stats.scoreatpercentile(allP,10) # 10th percentile, not 90th
		#bberdic[LDid] = {clustid:thisstat}
		bberP.append(thisstat)
		#bberP[numClust*(LDid-1) + (clustid-1)] = thisstat

# Sort, tie get the same ranks
bberRank = ranks(bberP) + 1 # ranks start from 0

# write out file
bberfile = open("plink.ADD.0.001.bber.linear", "w")
eqtlfile = open("plink.ADD.0.001.linear") # plink.ADD.linear
eqtldic = {} # dic of dic of array. LD -> genes cluster -> array
header = True
for line in eqtlfile:
	aline = line.strip().split()
	if header:
		header = False
		bberfile.write("\t".join(aline + ["BBER_P"]) + "\n")
		continue
	thisSNP = aline[2]
	thisGene = aline[0]
	thisP = float(aline[9])
	LDid = LDdic[thisSNP]
	clustid = clustdic[thisGene]
	# This is the BBER adjustment
	adjP = thisP*(numLD*numClust)/bberRank[numClust*(LDid-1) + (clustid-1)]
	aline.append(str(adjP))
	bberfile.write("\t".join(aline) + "\n")

bberfile.close()

