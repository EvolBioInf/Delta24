from sympy import *
import sys
import random
import numpy
import math
import time
import itertools

def simp_print (sym_eq, name, data, parms):
	eq=str(sym_eq)
	for x in range (0, len(data.names) ):
		eq=eq.replace(data.names[x], "X["+str(x)+"]")
	eq=eq.replace("log(0.333333333333333*E*(-E + 1))", "self.coef[1]")
	eq=eq.replace("log(0.0555555555555556*E**2 + 0.166666666666667*E*(-E + 1))", "self.coef[2]")
	eq=eq.replace("log(0.111111111111111*E**2)", "self.coef[3]")
	eq=eq.replace("log((-E + 1)**2)", "self.coef[4]")
	eq=eq.replace("log(0.166666666666667*E*(-E + 1) + 0.5*(-E + 1)**2)", "self.coef[5]")
	eq=eq.replace("log(0.0555555555555556*pow(parms[1],2)+0.5*pow(1-parms[1],2))", "self.coef[6]")
	eq=eq.replace("log(-0.333333333333333*E + 0.5)", "self.coef[7]")
	eq=eq.replace("log(0.333333333333333*E)", "self.coef[8]")

	eq=eq.replace("P**2", "pow(P, 2)")
	eq=eq.replace("(-P + 1)**2", "pow(1-P, 2)")
	eq=eq.replace("F**2", "pow(F, 2)")
	eq=eq.replace("E**2", "pow(E, 2)")
	eq=eq.replace("(-E + 1)**2", "pow(1-E, 2)")
	eq=eq.replace("1.0*", "")
	for x in range (0, len(parms) ):
		eq=eq.replace(parms[x].name, "parms["+str(x)+"].val")
	eq=eq.replace("log(0.0555555555555556*pow(parms[1].val, 2) + 0.5*pow(1-parms[1].val, 2))", "self.coef[6]")
	eq=eq.replace("(parms[2].val*(-pow(parms[0].val, 2) + parms[0].val) + pow(1-parms[0].val, 2))", "self.coef[10]")
	eq=eq.replace("(-parms[2].val + 1)*(-2*pow(parms[0].val, 2) + 2*parms[0].val)", "self.coef[11]")
	eq=eq.replace("(parms[2].val*(-pow(parms[0].val, 2) + parms[0].val) + pow(parms[0].val, 2))", "self.coef[12]")
	print name+"="+eq

def dif (A, B):
	X=int(A[0]!=B[0])+int(A[1]!=B[1])
	return X

def get (Nn1, Nn2):
	di=["Aa","At","Ac","Ag","Ta","Tt","Tc","Tg","Ca","Ct","Cc","Cg","Ga","Gt","Gc","Gg"]
	A=["A1","T1","C1","G1"]
	B=["A2","T2","C2","G2"]
	monolookup={'a':"A2",'t':"T2",'c':"C2",'g':"G2",'A':"A1",'T':"T1",'C':"C1",'G':"G1"}

	ret=["exp(0"]
	AR=[monolookup[Nn1[0]], monolookup[Nn2[0]]]
	if AR[0]==AR[1]:
		AR=[AR[0]]
	BR=[monolookup[Nn1[1]], monolookup[Nn2[1]]]
	if BR[0]==BR[1]:
		BR=[BR[0]]
	AW=[]
	BW=[]
	for a in A:
		if a not in AR:
			AW.append(a)
	for b in B:
		if b not in BR:
			BW.append(b)
	if len(AR)==1:
		AF="log(1-E)"
	else:
		AF="log((1-E)/2.0+E/6.0)"
	if len(BR)==1:
		BF="log(1-E)"
	else:
		BF="log((1-E)/2.0+E/6.0)"
	
	OmE=[]
	ESQ=[]
	OmEh=[]
	for d in di:
		if d==Nn2 or d==Nn1:
			continue
		elif dif(d, Nn1)==1 and dif(d, Nn2)==1:
			OmE.append(d)
		elif dif(d, Nn1)==1 or dif(d, Nn2)==1:
			OmEh.append(d)
		else:
			ESQ.append(d)
	if len(OmE)!=0:
		ret.append("+("+'+'.join(OmE)+")*log(E/3.0*(1-E))")
	if len(OmEh)!=0:
		ret.append("+("+'+'.join(OmEh)+")*log((1-E)*E/6.0+E**2/18.0)")
	ret.append("+("+'+'.join(ESQ)+")*log(E**2/9.0)")
	
	if dif(Nn1, Nn2)==0:
		ret.append("+("+Nn1+")*log((1-E)**2)")
		ret.append("+("+'+'.join(AR)+")*"+AF+"+("+'+'.join(AW+BW)+")*ln(E/3.0)+("+'+'.join(BR)+")*"+BF+")")
		return  "".join(ret)
	if dif(Nn1, Nn2)==1:
		ret.append("+("+Nn1+"+"+Nn2+")*log((1-E)**2/2.0+(1-E)*E/6.0)")
		ret.append("+("+'+'.join(AR)+")*"+AF+"+("+'+'.join(AW+BW)+")*ln(E/3.0)+("+'+'.join(BR)+")*"+BF+")")
		return  "".join(ret)
	if dif(Nn1, Nn2)==2:
		ret.append("+("+Nn1+"+"+Nn2+")*log((1-E)**2/2.0+E**2/18.0) ")
		ret.append("+("+'+'.join(AR)+")*"+AF+"+("+'+'.join(AW+BW)+")*ln(E/3.0)+("+'+'.join(BR)+")*"+BF+")")
		return  "".join(ret)
	print dif(Nn1, Nn2)
	print "Oh shit"
	quit()

#A class for parameters of the likelihood equation.
class p:
	def __init__ (self, name, LB, UB):
		self.name=name
		self.sym=Symbol(name)
		self.val=0
		self.UB=UB
		self.LB=LB
		self.R=UB-LB
		if self.R<=0:
			print "range error: parameter "+name+" can take no values"
			quit()

#A class for the data/observations.
class d:
	def __init__ (self, D):
		self.D=len(D)
		self.Ds=[]
		self.sym=[]
		self.names=[]
		for ps in D:
			self.Ds.append(p(ps[0], ps[1], ps[2]) )
		for ps in self.Ds:
			self.sym.append(ps.sym)
			self.names.append(ps.name)

	def read (self, filename):
		File=open(filename)
		self.data={}
		self.sets=[]
		for name in self.names:
			self.data[name]=[]
		for line in File:
			if line[0]=='>':
				continue
			datum=[]
			line=line.strip('\n').split('\t')
			for x in range(0, self.D):
				t=float(line[x])
				if t<=self.Ds[x].UB and t>=self.Ds[x].LB:
					datum.append(t)
					self.data[self.names[x]].append(t)
			self.sets.append(datum)
		File.close()
		
	def left_sort (self, breaks):
		breaks.insert(0, 0)
		breaks.append(len(self.sets[0]))
		for x in range (0, len(self.sets) ):
			slices=[]
			for y in range(0, len(breaks)-1):
				slices.append(self.sets[x][breaks[y]:breaks[y+1]] )
			for aslice in slices:
				aslice.sort()
			self.sets[x]=[]
			for aslice in slices:
				self.sets[x]+=aslice
	def compress (self):
		count_keys={}
		self.counts=[]
		for X in self.sets:
			key=','.join(map(str, X) )
			try:
				count_keys[key][1]+=1
			except:
				count_keys[key]=[X, 1]
		self.sets=[]
		for key in count_keys.keys():
			self.sets.append(count_keys[key][0])
			self.counts.append([count_keys[key][1]])
				
#to get the vector which increments the parameters.

class system : 
	def __init__ (self, eqs, parms, data):
		self.eqs=eqs
		self.N=len(eqs)
		self.parms=parms
		self.data=data
		self.R=float("infinity")
		self.diff=[]
		for x in range(0, self.N):
			self.diff.append([])
			for y in range(0, self.N) :
				self.diff[x].append(diff(self.eqs[x], parms[y].sym) )
		if len(parms)!=self.N:
			print "system does not have one equation for each parameter. Matt cannot maximize this."
			quit()
		#TODO figure out if I can check to see if any of the equations have algabaic solutions.
	def simp (self, sym_eq, name):
		eq=str(sym_eq)
		self.coef=[0,0,0,0,0,0,0,0,0,0,0,0,0]
		for x in range (0, len(self.data.names) ):
			eq=eq.replace(self.data.names[x], "X["+str(x)+"]")
                eq=eq.replace("log(0.333333333333333*E*(-E + 1))", "self.coef[1]")
                eq=eq.replace("log(0.0555555555555556*E**2 + 0.166666666666667*E*(-E + 1))", "self.coef[2]")
                eq=eq.replace("log(0.111111111111111*E**2)", "self.coef[3]")
                eq=eq.replace("log((-E + 1)**2)", "self.coef[4]")
                eq=eq.replace("log(0.166666666666667*E*(-E + 1) + 0.5*(-E + 1)**2)", "self.coef[5]")

                eq=eq.replace("log(-0.333333333333333*E + 0.5)", "self.coef[7]")
                eq=eq.replace("log(0.333333333333333*E)", "self.coef[8]")
	        eq=eq.replace("log(-E + 1)", "self.coef[9]")
	        eq=eq.replace("P**2", "pow(P, 2)")
	        eq=eq.replace("(-P + 1)**2", "pow(1-P, 2)")
	        eq=eq.replace("F**2", "pow(F, 2)")
	        eq=eq.replace("E**2", "pow(E, 2)")
	        eq=eq.replace("(-E + 1)**2", "pow(1-E, 2)")
		for x in range (0, len(self.parms) ):
			eq=eq.replace(self.parms[x].name, "parms["+str(x)+"].val")
		print name+"="+eq
		E=self.parms[1].val
	        self.coef[1]=math.log(0.333333333333333*E*(-E + 1)) #done
                self.coef[2]=math.log(0.0555555555555556*E**2 + 0.166666666666667*E*(-E + 1)) #done
                self.coef[3]=math.log(0.111111111111111*E**2) #done
                self.coef[4]=math.log((1-E)**2) #done
                self.coef[5]=math.log(0.166666666666667*E*(-E + 1) + 0.5*(-E + 1)**2)
                self.coef[6]=math.log(0.111111111111111*E**2 + 0.5*(-E + 1)**2) #done
                self.coef[7]=math.log((1-E)**2/2.0+E**2/9.0) #done
                self.coef[8]=math.log(0.333333333333333*E) #done
                self.coef[9]=math.log(-E + 1) #done
		print map(str, self.coef)
		return eq
	def eq_eval (self, str_eq, X):
		E=self.parms[1].val
		exec("val="+str_eq)
		return val
	def inc (self):
		self.Rval=zeros((self.N, 1) )
		self.Req=[]
		self.Jval=zeros(self.N)
		test=zeros(self.N)
		self.Jeq=[]
		for x in range(0, self.N):
			self.Jeq.append([])
			self.Req.append(self.simp(str(self.eqs[x]), "R["+str(x)+"]") )
			for y in range(0, self.N):
				K=self.diff[x][y]
				self.Jeq[x].append(self.simp(str(K), "J["+str(x)+"]["+str(y)+"]" ) )
		print "entering main loop"
		for X, Y in itertools.izip(self.data.sets,self.data.counts):
			for x in range(0, self.N):
				for y in range(0, self.N):
					r=self.eq_eval(self.Jeq[x][y], X)
					self.Jval[x,y]+=r*Y[0]
				try:
					r=self.eq_eval(self.Req[x], X)
					self.Rval[x]+=r*Y[0]
				except:
					self.Rval[x]+=float("infinity")
		IJ=self.Jval.inv()
		self.R=0
		for x in range(0, self.N):
			test[x]=self.parms[x].val
		for x in range(0, self.N):
			for y in range(0, self.N):
				test[x]-=(IJ[x,y]*self.Rval[y])
			self.R+=abs(self.Rval[x])
		for x in range(0, self.N):
			if test[x]>self.parms[x].LB and test[x]<self.parms[x].UB:
				self.parms[x].val=test[x]
			else:
				print "MLE fail to converge : setting "+self.parms[x].name+" to "+str(test[x])

#This is a general class for an equation. It associates an equation with data and the parameters of the equation.
class eq :
	def __init__ (self, eq, parms, data):
		self.string=eq
		for x in range (0, len(parms) ):
			eq=eq.replace(parms[x].name, "parms["+str(x)+"].sym")
		for x in range (0, len(data.names) ):
			eq=eq.replace(data.names[x], "data.sym["+str(x)+"]")
		self.parms=parms
		self.data=data
		self.N=len(parms)
		self.T=""
		exec("self.eq="+eq)
	def setfd(self):
		self.T=self.eq
		for parm in self.parms:
			self.T=self.T.subs(parm.sym, parm.val)
		self.fd=lambdify(self.data.sym, self.T)
		

#A class for the Log likelihood equation. Bad things will happen if you give it the likelihood equation itself.
class Leq :
	def __init__ (self, eq, parms, data, surf_size):
		print "initializing likelihood maximization equations, this may take some time."
		self.string=eq
		for x in range (0, len(parms) ):
			eq=eq.replace(parms[x].name, "parms["+str(x)+"].sym")
		for x in range (0, len(data.names) ):
			eq=eq.replace(data.names[x], "data.sym["+str(x)+"]")
		self.parms=parms
		self.data=data
		self.xeq=[]
		self.N=len(parms)
		exec("self.eq="+eq)
		for x in range(0, self.N):
			deq=diff(self.eq, parms[x].sym)
			self.xeq.append(deq)
		self.S=system(self.xeq, parms, data)
		self.R=self.S.R

		self.L=0
		self.I=0
		self.E=0
		print "done."
	#this is unfortunately hand typed, and will not be automatically updated with the program.
	def set_ln (self):
		eq=self.string
		self.coef=[0,0,0,0,0,0,0]
		for x in range (0, len(self.parms) ):
			eq=eq.replace(self.parms[x].name, "parms["+str(x)+"].val")
		for x in range (0, len(self.data.names) ):
			eq=eq.replace(self.data.names[x], "X["+str(x)+"]")
		eq=eq.replace("log(E/3.0*(1-E))", "self.coef[1]")
		eq=eq.replace("log((1-E)*E/6.0+E**2/18.0)", "self.coef[2]")
		eq=eq.replace("log(E**2/9.0)", "self.coef[3]")
		eq=eq.replace("log((1-E)**2)", "self.coef[4]")
		eq=eq.replace("log((1-E)**2/2.0+(1-E)*E/6.0)", "self.coef[5]")
		eq=eq.replace("log((1-E)**2/2.0+E**2/9.0)", "self.coef[6]")
		E=self.parms[1].val
		self.coef[1]=log(E/3.0*(1-E))
		self.coef[2]=log((1-E)*E/6.0+E**2/18.0)
		self.coef[3]=log(E**2/9.0)
		self.coef[4]=((1-E)**2)
		self.coef[5]=log((1-E)**2/2.0+(1-E)*E/6.0)
		self.coef[6]=log((1-E)**2/2.0+E**2/9.0)
		self.eval_eq=eq
	def ln_eval (self, X):
#		self.coef[0]=lg(1+X[0]+X[1]+X[2]+X[3]+X[4]+X[5]+X[6]+X[7])-lg(1+X[0])-lg(1+X[1])-lg(1+X[2])-lg(1+X[3])-lg(1+X[4])-lg(1+X[5])-lg(1+X[6])-lg(1+X[7])
		exec("val="+self.eval_eq)
		return val

	def get_L (self):
		self.L=0
		self.set_ln()
		for X, Y in itertools.izip(self.data.sets,self.data.counts):
			self.L+=self.ln_eval(X)*Y[0]
		return self.L

	def get_Ls (self):
		self.L=0
		self=self.eq
		for parm in self.parms:
			T=T.subs(parm.sym, parm.val)
		self.fd=lambdify(self.data.sym, T)
		for x in self.data.sets[0:100]:
			try:
				self.L+=self.fd(*x)
			except:
				self.L+=float("-infinity")
		return self.L

	def get_f (self):
		T=self.eq
		for parm in self.parms:
			T=T.subs(parm.sym, parm.val)
		self.fd=lambdify(self.data.sym, T)

	def f (self, x):
		return self.fd(x)
	
	def inc(self):
		self.S.inc()
		self.R=self.S.R
		ret=[]
		for p in self.parms:
			ret.append(str(p.name)+":"+str(p.val))
		return '\t'.join(ret)

#read in the data from sys.argv[1], 1 datum per line, 0 and "infinity" are the upper and lower bounds for the data
#right now it is set up to read dinucleotide data. 

data=d([["A1", 0, float("infinity")],["C1", 0, float("infinity")],["G1", 0, float("infinity")],["T1", 0, float("infinity")],["A2", 0, float("infinity")],["C2", 0, float("infinity")],["G2", 0, float("infinity")],["T2", 0, float("infinity")],["Aa", 0, float("infinity")],["Ac", 0, float("infinity")],["Ag", 0, float("infinity")],["At", 0, float("infinity")],["Ca", 0, float("infinity")],["Cc", 0, float("infinity")],["Cg", 0, float("infinity")],["Ct", 0, float("infinity")],["Ga", 0, float("infinity")],["Gc", 0, float("infinity")],["Gg", 0, float("infinity")],["Gt", 0, float("infinity")],["Ta", 0, float("infinity")],["Tc", 0, float("infinity")],["Tg", 0, float("infinity")],["Tt", 0, float("infinity")]])

#data.read("test.txt")
#data.compress()

#this would speed things up, but right now it messes up all our calculations.

#declare all the parameters we want the ML procedure to maximize.
#here P is Pi, E is Epsilon and F is delta. 

P=p("P", 0.00001, 1.0 )
E=p("E", 0.00001, 1.0 )
F=p("F", 0.0, 1.0)

#The paramaters have to be put into a list to pass it to the main maximization procedure.
parms=[P, E, F]

#And we have to set reasonable initial values for the parameters or our estimates may not converge. Which is bad.
P.val=0.01
E.val=0.01
F.val=0.01

DI={}
H00_DI=[]
H01_DI=[]
H11_DI=[]
N=["Aa","At","Ac","Ag","Ta","Tt","Tc","Tg","Ca","Ct","Cc","Cg","Ga","Gt","Gc","Gg"]
E00_DI=[]
E01_DI=[]
E11_DI=[]


for x in range(0, 16):
	for y in range(x, 16):
		Nn1=N[x]
		Nn2=N[y]
		DI[Nn1+Nn2]=get(Nn1, Nn2)
		
		if Nn1==Nn2:
			H00_DI.append(get(Nn1, Nn2) )
		elif Nn1[0]==Nn2[0] or Nn1[1]==Nn2[1]:
			H01_DI.append(get(Nn1, Nn2) )
		else:
			H11_DI.append(get(Nn1, Nn2) )

#Our Likelihood equation. It should apper as python wants to see it, except in quotes.
#The equation itslf is broken into three sections, section one (proceded by (1-P)**2 ) calculates the likelihood of the observation given that the observation is made at a double homozygous site,
#section two (proceded by (1-F)**2 ) calculates the likelihood of the observation given that the observation is made at a single heterozygote site,

H00=eq( "((P-P**2)*F+(1-P)**2)*("+'+'.join(H00_DI)+")/16.0", parms, data)

H01=eq( "2*(P-P**2)*(1-F)*("+'+'.join(H01_DI)+")/48.0", parms, data)

H11=eq( "((P-P**2)*F+P**2)*("+'+'.join(H11_DI)+")/72.0", parms, data)

Eeq=eq( "((P-P**2)*F+(1-P)**2)*("+'+'.join(E00_DI)+")/16.0+2*(P-P**2)*(1-F)*("+'+'.join(E01_DI)+")/48.0+((P-P**2)*F+P**2)*("+'+'.join(E11_DI)+")/72.0", parms, data)

MLRoe=Leq( "ln("+H00.string+"+"+H01.string+"+"+H11.string+")", parms, data, 0)

simp_print (diff(MLRoe.eq, parms[2].sym), "F0", data, parms)

quit()

