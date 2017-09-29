#! /usr/bin/python3
# Luis del Peso
# Dec 2015
# Implemmentation of the Needleman-Wunch for global sequence alignment
# Based on the pseudocode on:
# http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
# run as:
# python3 Needleman_Wunsch_v4.py Scoring_Matrix Sequnece1 Sequence2 gap_penalty

def Fmatrix(Seq1,Seq2,SubMat,GapPen):# Generates the Score Matrix (F matrix)
    ScoreMatrix=[] #Initializes empty Score matrix
    #Following lines generate an dummy matrix containing only gaps values
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
            Match=ScoreMatrix[i-1][j-1]+SubMat[Seq1[i-1]][Seq2[j-1]] #Score match
            Delete=ScoreMatrix[i-1][j]+GapPen #Score for deletion (gap in seq2)
            Insert=ScoreMatrix[i][j-1]+GapPen #Score for insertion (gap in seq1)
            ScoreMatrix[i][j]=max(Match,Delete,Insert) #cell value is the max
    return(ScoreMatrix)

def GetSM(FileName):#reads and returns Scoring Matrix
    SubstMatrix={}#Initializes empty dictionary (elememts will be dictionaries)
    MyFile=open(FileName,"r")
    header=MyFile.readline()
    alphabet=header.strip().split("\t") #reads alphabet from first line
    for Line in MyFile:#reads all matrix rows
        tmp=Line.strip().split("\t")
        tmpDict={}#dictionary representing each row
        for pos in list(range(1,len(tmp))):#note that skips first position
            tmpDict[alphabet[pos]]=int(tmp[pos])#fills tmp dictionary entries
        SubstMatrix[tmp[0]]=tmpDict#adds element (a dictionary) per row
    return(SubstMatrix)

def GetSeq(FileName):
    Seq=""
    MyFile=open(FileName,"r")
    for Line in MyFile:#reads all matrix rows
        if Line.count(">")>0:
          next
        else:
            Seq=Seq+Line.strip()
    return(Seq.upper())

def PrintMatrix(Matrix,Seq1,Seq2):# just shows the matrix...
    print("*\t\t",end="")
    for residue in Seq1:
        print(residue,"\t",end="")
    print()
    for i in list(range(len(resFmatrix[0]))):
        if i>0:
            print(Seq2[i-1],"\t",end="")
        else:
            print("*\t",end="")
        for j in list(range(len(resFmatrix))):
            print(resFmatrix[j][i],"\t",end="")
        print()

def Alignment(ScoreMatrix,Seq1,Seq2,SubMat,GapPen):#generates the alignment
    AlSeq1=""#initializes aligned sequence 1 as an empty string
    AlSeq2=""#initializes aligned sequence 2 as an empty string
    i=len(Seq1)
    j=len(Seq2)
    while (i>0 or j>0):#goes through the Fmatrix in reverse order
        #3 next lines calculate all three possible paths to current position
        Match=ScoreMatrix[i-1][j-1]+SubMat[Seq1[i-1]][Seq2[j-1]]
        Delete=ScoreMatrix[i-1][j]+GapPen
        Insert=ScoreMatrix[i][j-1]+GapPen
        if (i>0 and j>0 and ScoreMatrix[i][j]==Match):#current pos from diagonal
            AlSeq1=Seq1[i-1]+AlSeq1
            AlSeq2=Seq2[j-1]+AlSeq2
            i=i-1
            j=j-1
        elif (i>0 and ScoreMatrix[i][j]==Delete):#current pos from gap in seq2
            AlSeq1=Seq1[i-1]+AlSeq1
            AlSeq2="-"+AlSeq2
            i=i-1
        elif (j>0 and ScoreMatrix[i][j]==Insert):#current pos from gap in seq1
            AlSeq1="-"+AlSeq1
            AlSeq2=Seq2[j-1]+AlSeq2
            j=j-1
    return(AlSeq1,AlSeq2)
    
import sys
SubstMatrix=GetSM(sys.argv[1])#read Scoring Matrix
print(SubstMatrix)
Seq1=GetSeq(sys.argv[2])
Seq2=GetSeq(sys.argv[3])
GapPenalty=int(sys.argv[4])
resFmatrix=Fmatrix(Seq1,Seq2,SubstMatrix,GapPenalty)
#print(resFmatrix)
PrintMatrix(resFmatrix,Seq1,Seq2)
AlignSeq1,AlignSeq2=Alignment(resFmatrix,Seq1,Seq2,SubstMatrix,GapPenalty)
print("Seq1\t",AlignSeq1)
print("Seq2\t",AlignSeq2)

