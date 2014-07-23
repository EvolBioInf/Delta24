import sys
import random
import numpy.random

r=0

def rcomp(seq):
	base={'A':'T', 'C':'G','G':'C','T':'A', '':'','AA':'TT', 'CC':'GG','GG':'CC','TT':'AA', '-':'-'}
	rc=[]
	for s in seq:
		rc.insert(0, base[s])	
	return ''.join(rc)
		

def rand_base(b):
	base={'A':['C','G','T', '', 'AA'], 'C':['A','G','T', '', 'CC'], 'G':['A','C','T', '', 'GG'], 'T':['A','C','G', '', 'TT'], '-':['-']}
#	print random.randint(0, 2)
	e=random.randint(0, 2)
	return base[b][e]
	
def printread(File1, File2, seq, start, start_mate, I, k, REV, F):
	seq=seq.translate(None, '-')
	seq_p=[]
	for s in seq:
		seq_p.append(s)
	for e in k:
		if (e<len(seq)):
			seq_p[e]=rand_base(seq[e])
	global r
	seq_p=''.join(seq_p)
	if(REV):
		File1.write("@B81180ABBXX:1:6:"+str(1800+r)+":"+str(2000+r)+"#GCCAATAT/"+str(F)+"\n"+seq_p+"\n"+"+\n"+"M"*len(seq_p)+'\n')
		File2.write("B81180ABBXX:1:6:"+str(1800+r)+":"+str(2000+r)+"#GCCAATAT/1\t2\tcontig_1\t"+str(start)+"\t100\t"+str(len(seq_p) )+"M\t=\t"+str(start_mate)+"\t"+str(I)+'\t'+seq_p+"\t"+"M"*len(seq_p)+'\n')
	else:
		File1.write("@B81180ABBXX:1:6:"+str(1800+r)+":"+str(2000+r)+"#GCCAATAT/"+str(F)+"\n"+rcomp(seq_p)+"\n"+"+\n"+"M"*len(seq_p)+'\n')
		File2.write("B81180ABBXX:1:6:"+str(1800+r)+":"+str(2000+r)+"#GCCAATAT/2\t2\tcontig_1\t"+str(start)+"\t100\t"+str(len(seq_p) )+"M\t=\t"+str(start_mate)+"\t"+str(I)+'\t'+seq_p+"\t"+"M"*len(seq_p)+'\n')
	r+=1

E=0.01
L=100
I=50
T=int(sys.argv[2])

Joint=True

seq={}

def dsum (O, T):
	return [O[0]+T[0], O[1]+T[1], O[2]+T[2], O[3]+T[3] ]

infile=open(sys.argv[1])

for line in infile:
	line=line.strip('\n')
	if line[0]!="#":
		if line[0]=='>':
			name=line[1:]
			seq[name]=[]
		else:
			seq[name].append(line)
for name in seq.keys():
	seq[name]=''.join(seq[name])

S1=seq["S1"]
S2=seq["S2"]

#S3=seq["S3"]

File=open("ref.fasta", 'w')
File.write(">contig_1\n")
#File.write(S3+'\n')
File.close()

diff=0
for x in range(0, len(S1) ):
	if S1[x]!=S2[x]:
		diff+=1
print float(diff)/float(len(S1))


K1S=[]
K2S=[]
K1E=[]
K2E=[]

File1=open("seq1.fq", 'w')
File2=open("seq2.fq", 'w')
File3=open("true-align.sam", 'w')

amb={'A':{'A':'A', 'C':'M', 'G':'R','T':'W'}, 'C':{'C':'C', 'A':'M', 'G':'S', 'T':'Y'}, 'G':{'G':'G', 'A':'R','C':'S', 'T':'K'}, 'T':{'T':'T', 'A':'W', 'C':'Y', 'G':'K'} }
File4=open("true.fasta", 'w')
File4.write(">contig_1\n")
out=[]
for x in range(0, len(S1) ):
	if S1[x]!='-' and S2[x]!='-':
		out.append(amb[S1[x]][S2[x]])
File4.write(''.join(out))
File4.close()

File3.write("@HD\tVN:1.0\tSO:unsorted\n")
File3.write("@SQ\tSN:contig_1\tAS:silico_silicae\tLN:"+str(len(S1))+'\n')

for g in range(0, len(S1)/(L*2)*T):
	x=random.randint(0, len(S1)-2*L-I)
	S=random.randint(0, 1)
	if S==0:
		errors=numpy.random.multinomial((L),[E, 1-E])[0]
		REV=(random.randint(0, 2) )
	#	REV=False
		k=[]
		for e in range(0, errors):
			k.append(random.randint(0, L-1) )
		printread(File1, File3, S1[x:x+L], x, x+L+I, I, k, REV, "1")
		errors=numpy.random.multinomial((L),[E, 1-E])[0]

		k=[]
		for e in range(0, errors):
			k.append(random.randint(0, L-1) )
	#	REV=True
		printread(File2, File3, S1[x+L+I:x+2*L+I], x+L+I, x, I, k, REV, "2")

		K1S.append(x)
		K1E.append(x+L)
		K1S.append(x+L+I)
		K1E.append(x+2*L+I)
	else:
		REV=(random.randint(0, 2) )
	#	REV=False
		errors=numpy.random.multinomial((L),[E, 1-E])[0]

		k=[]
		for e in range(0, errors):
			k.append(random.randint(0, L-1) )
		printread(File1, File3, S2[x:x+L], x, x+L+I, I, k, REV, "1")
		errors=numpy.random.multinomial((L),[E, 1-E])[0]

		k=[]
	#	REV=True
		for e in range(0, errors):
			k.append(random.randint(0, L-1) )
		printread(File2, File3, S2[x+L+I:x+2*L+I], x+L+I, x, I, k, REV, "2")

		K2S.append(x)
		K2E.append(x+L)
		K2S.append(x+L+I)
		K2E.append(x+2*L+I)

K1S.sort(reverse=True)
K2S.sort(reverse=True)
K1E.sort(reverse=True)
K2E.sort(reverse=True)

k1s=K1S.pop()
k2s=K2S.pop()
k1e=K1E.pop()
k2e=K2E.pop()

K1=0
K2=0

base={'A':[(1-E), E/3, E/3, E/3],'C':[E/3, (1-E), E/3, E/3], 'G':[E/3, E/3, (1-E), E/3],'T':[E/3, E/3, E/3, (1-E)]}

File=open("true-test.txt", 'w')
File.write(">contig_1\n")
for x in range(0, len(S1) ):
	while x==k1s:
		K1+=1
		try:
			k1s=K1S.pop()
		except:
			k1s=len(S1)
	while x==k2s:
		K2+=1
		try:
			k2s=K2S.pop()
		except:
			k2s=len(S1)
	while x==k1e:
		K1-=1
		try:
			k1e=K1E.pop()
		except:
			k1e=len(S1)
	while x==k2e:
		K2-=1
		try:
			k2e=K2E.pop()
		except:
			k2e=len(S1)
	if (Joint):
		if S1[x]!='-':
			M=numpy.random.multinomial(K1, base[S1[x]])  
		if S2[x]!='-':
		        D=numpy.random.multinomial(K2, base[S2[x]])  
	else:
		if S1[x]!='-':
			M=numpy.random.multinomial((K1+K2)/2.0, base[S1[x]])  
		if S2[x]!='-':
	        	D=numpy.random.multinomial((K1+K2)/2.0, base[S2[x]])  
	K=dsum(M, D)
        File.write(str(x)+'\t'+'\t'.join(map(str, K) )+'\t') 
        File.write('\n')
File.close()
