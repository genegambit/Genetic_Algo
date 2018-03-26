import os, sys
from numpy import *
from numpy import equal
import blosum
from blosum import BLOSUM50 as B50
from copy import *
import random         


def BlosumScore( mat, abet, s1, s2, gap=-8 ):
    sc = 0
    n = min( [len(s1), len(s2)] )
    for i in range( n ):
        if s1[i] == '-' or s2[i] == '-' and s1[i] != s2[i]:
            sc += gap
        elif s1[i] == '.' or s2[i] == '.':
            pass
        else:
            n1 = abet.index( s1[i] )
            n2 = abet.index( s2[i] )
            sc += mat[n1,n2]
    return sc
    
def EstablishTargets():
    tgsc1 = BlosumScore( blosum.BLOSUM50, blosum.PBET, seq1, seq1 )
    tgsc2 = BlosumScore( blosum.BLOSUM50, blosum.PBET, seq2, seq2 )
    cxsc2 = BlosumScore( blosum.BLOSUM50, blosum.PBET, seq1, seq2 )
    bestsc = max( (tgsc1, tgsc2) ) + cxsc2
    return bestsc
    
def Jumble( PBET, ngenes ):
    folks = []
    ape = (deepcopy( PBET ))*5
    for i in range( ngenes ):
	ape1 = list(ape)
        random.shuffle( ape1 )
        folks.append( deepcopy( ape1 ))
    return folks

def SingleCost( folk, seq1, seq2 ):
    number1 = BlosumScore( blosum.BLOSUM50, blosum.PBET, seq1, folk )
    number2 = BlosumScore( blosum.BLOSUM50, blosum.PBET, seq2, folk )
    singlecost = Bestsc - (number1 + number2)
    return singlecost
    
def CostFunction (seq):
    Cost=[]
    for count in range(len(seq)):
	CostEach = SingleCost(seq[count],seq1,seq2)
	Cost.append(float(CostEach))
    return array(Cost)

def CrossOver( folks, Cost ):
    # convert costs to probabilities
    dim = len( folks[0] )
    prob = Cost + 0.0
    mx = prob.max()
    prob = mx - prob	# lowest cost is now highest numbr.
    mx = prob.max()
    prob = prob / mx	# makes sure numbers aren't too high
    prob = prob / prob.sum()	# normalized.  sum(prob) = 1.0
    # make new kids
    kids = []
    NG = len( folks )
    for i in range( NG/2 ):
        rdad = random.random()
        rmom = random.random()
        # find which vectors to use
        sm = 0.0
        idad = 0
        while rdad > sm:
            sm = sm + prob[idad]
            idad = idad + 1
        sm = 0.0
        imom = 0
        while rmom > sm:
            sm = sm + prob[imom]
            imom = imom+1
        idad,imom = idad-1,imom-1
        # make babies
        x = int(random.random()*(dim-2))+1	# crossover
        kids.append(concatenate((folks[idad][:x],folks[imom][x:])))
        kids.append(concatenate((folks[imom][:x],folks[idad][x:])))
    return kids

def Mutate(Kids):
    SecondKids = deepcopy(Kids)
    mutationlength = 20
    l = len (SecondKids[0])
    lim = random.randint(0,mutationlength+1)
    startposition = random.randint(0,l- lim+1)
    for i in range (len(SecondKids)):
	#kiddd = deepcopy(FirstKids[1])
	aa = SecondKids[i][(startposition) : (startposition + lim)]
	random.shuffle(aa)
	SecondKids[i][(startposition) : (startposition + lim)] = aa
	return SecondKids


def Feud( folks, kids, fcost, kcost ):

    for i in range( 0, len(kids) ):

        if kcost[i] < fcost[i]:

            folks[i] = kids[i]

            fcost[i] = kcost[i]
    return folks
     
######################################################################################################################################   
seq1 = 'MNFTSLLQDGIYEVGNGAIVTDQSPYLGITPDYQGAYGFPTHPWGIFNKAKAKAAGFQVVGAILVFGAYLPAVIKVLISKRTENLAIGMWIISIAGLGLL'
seq2 = 'AIFAWLGVSVNPGGFILVALSETLSCIASIIVFALKIANKAKAKAAGMTELEYCNLNKAKAKAAGHYPIVKKLPKRDGIYEVGNGAIVTDQSPYLGDGIY'

folks = Jumble(blosum.PBET, 100)
sc = BlosumScore(B50,blosum.PBET,seq1,seq2)
Bestsc = EstablishTargets()
B_Score = Bestsc
for i in range (0, 10000):
    CostFolks = CostFunction(folks)
    #print len(CostFolks)
    mc = CostFolks.min()
    if mc < B_Score:
	mc_index = [x for x in range(len(CostFolks)) if CostFolks[x] == mc][0]
	Best = folks[mc_index]
	B_Score = mc
    FirstKids = CrossOver( folks , CostFolks )
    SecondKids = Mutate(FirstKids)
    folks = Mutate(FirstKids)
    print i
    CostSecondKids = CostFunction(SecondKids)
    folks = Feud (folks, SecondKids, CostFolks, CostSecondKids)
    

FCost = CostFunction (folks)
#print FCost
print "BEST SCORE: ", `B_Score`
print "BEST SEQUENCE: ", `Best`