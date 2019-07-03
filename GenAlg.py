import os, sys
from numpy import *
from numpy import equal
import blosum
from blosum import BLOSUM50 as B50
from copy import *
import random 
import time
start_time = time.time()


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
    tgsc1 = BlosumScore (blosum.BLOSUM50, blosum.PBET, seq1, seq1)
    tgsc2 = BlosumScore (blosum.BLOSUM50, blosum.PBET, seq2, seq2)
    cxsc2 = BlosumScore (blosum.BLOSUM50, blosum.PBET, seq1, seq2)
    bestsc = max ((tgsc1, tgsc2) ) + cxsc2
    return bestsc
    
def Jumble( PBET, ngenes ):
    folks = []
    ape = (deepcopy (PBET))*5
    for i in range (ngenes):
        ape1 = list (ape)
        random.shuffle (ape1)
        folks.append (deepcopy(ape1))
    return folks

def SingleScore( folk, seq1, seq2 ):
    number1 = BlosumScore (blosum.BLOSUM50, blosum.PBET, seq1, folk)
    number2 = BlosumScore (blosum.BLOSUM50, blosum.PBET, seq2, folk)
    singlecost = Bestsc - (number1 + number2)
    return singlecost
    
def ScoreFunction (seq):
    Score = []
    for count in range (len(seq)):
        ScoreEach = SingleScore (seq[count],seq1,seq2)
        Score.append (float (ScoreEach))
    return array (Score)

def CrossOver (folks, Cost):
    # convert costs to probabilities
    dim = len (folks[0])
    prob = Cost + 0.0
    mx = prob.max()
    prob = mx - prob 
    mx = prob.max()
    prob = prob / mx   
    prob = prob / prob.sum() 
    # make new kids
    kids = []
    NG = len(folks)
    NG_half = int(NG/2)
    for i in range (NG_half):
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
        x = int(random.random()*(dim-2))+1  # crossover
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
iteration = 10000
for i in range (0, iteration):
    ScoreFolks = ScoreFunction(folks)
    mc = ScoreFolks.min()
    if mc < B_Score:
        mc_index = [x for x in range(len(ScoreFolks)) if ScoreFolks[x] == mc][0]
        Best = folks[mc_index]
        B_Score = mc
    FirstKids = CrossOver (folks , ScoreFolks)
    SecondKids = Mutate (FirstKids)
    folks = Mutate (FirstKids)
    print (i)
    ScoreSecondKids = ScoreFunction (SecondKids)
    folks = Feud (folks, SecondKids, ScoreFolks, ScoreSecondKids)
    

FCost = ScoreFunction (folks)
print ("Sequence A: ", seq1)
print ("Sequence B: ", seq2)
print ("Best Score: ", B_Score)
print ("Best Match: ", ''.join(Best))
print ("Iterations: ", iteration)
end_time = round((time.time() - start_time)/60, 2)
print("Run Time: ", str(end_time) + " Minutes.")
