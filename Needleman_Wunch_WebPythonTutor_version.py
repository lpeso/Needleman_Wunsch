#! /usr/bin/python3
# Luis del Peso
# November 2014
# Implemmentation of the Needleman-Wunch for global sequence alignment
# Based on the pseudocode in:
# http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

def Fmatrix(Seq1,Seq2,SubMat,GapPen):# Generates the Score Matrix (F matrix)
    ScoreMatrix=[] #Initializes empty Score matrix
    for i in list(range(len(Seq1)+1)):#Initializes first row (Gaps)
        tmp=[]
        if i==0:
            for j in list(range(len(Seq2)+1)):#Initializes first column (Gaps)
                tmp.append(GapPen*j)#adds a number
            ScoreMatrix.append(tmp)#adds a vector of length=length(Seq2)
        else:
            tmp.append(GapPen*i)#adds firt value of column (row values)
            for j in list(range(1,len(Seq2)+1)):#fills cells with "*"
                tmp.append("*")
            ScoreMatrix.append(tmp)
    #Calculates cell values. Note that the range goes from +1 to length+1
    for i in list(range(1,len(Seq1)+1)):
        for j in list(range(1,len(Seq2)+1)):            
            Match=ScoreMatrix[i-1][j-1]+SubMat[Seq1[i-1]][Seq2[j-1]]
            Delete=ScoreMatrix[i-1][j]+GapPen
            Insert=ScoreMatrix[i][j-1]+GapPen
            ScoreMatrix[i][j]=max(Match,Delete,Insert)
    return(ScoreMatrix)

def Alignment(ScoreMatrix,Seq1,Seq2,SubMat,GapPen):
    AlSeq1=""
    AlSeq2=""
    i=len(Seq1)
    j=len(Seq2)
    while (i>0 or j>0):
        Match=ScoreMatrix[i-1][j-1]+SubMat[Seq1[i-1]][Seq2[j-1]]
        Delete=ScoreMatrix[i-1][j]+GapPen
        Insert=ScoreMatrix[i][j-1]+GapPen
        if (i>0 and j>0 and ScoreMatrix[i][j]==Match):
            AlSeq1=Seq1[i-1]+AlSeq1
            AlSeq2=Seq2[j-1]+AlSeq2
            i=i-1
            j=j-1
        elif (i>0 and ScoreMatrix[i][j]==Delete):
            AlSeq1=Seq1[i-1]+AlSeq1
            AlSeq2="-"+AlSeq2
            i=i-1
        elif (j>0 and ScoreMatrix[i][j]==Insert):
            AlSeq1="-"+AlSeq1
            AlSeq2=Seq2[j-1]+AlSeq2
            j=j-1
    return(AlSeq1,AlSeq2)
    
SubstMatrix={'T': {'T': 1, 'G': -1, 'A': -1, 'C': -1}, 'G': {'T': -1, 'G': 1, 'A': -1, 'C': -1}, 'A': {'T': -1, 'G': -1, 'A': 1, 'C': -1}, 'C': {'T': -1, 'G': -1, 'A': -1, 'C': 1}}
GapPenalty=-1
Seq1="ACAGT"
Seq2="ACGT"
resFmatrix=Fmatrix(Seq1,Seq2,SubstMatrix,GapPenalty)
AlignSeq1,AlignSeq2=Alignment(resFmatrix,Seq1,Seq2,SubstMatrix,GapPenalty)
print("Seq1\t",AlignSeq1)
print("Seq2\t",AlignSeq2)