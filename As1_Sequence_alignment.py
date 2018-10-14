# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 14:31:33 2018

@author: wangyang
sovle problem
"
Implementation the classical affine-gap local alignment algorithm, 
for sequence over the 20-letter protein alphabet. The Blosum scoring matrix can be found in attachment.  
For the data attached there is a collection of two sets of 50 protein sequences, 
each divided into 25 pairs of putatively homologous sequences; 
one is pairs of sequences from human and mouse, 
while the other is from pairs of sequences from human and fruit flies.
"
the algorithms is dynamic programming .
this BLOSUM62-As1 file name exict Blosum Scoring matrix.  
"""
import numpy as np
class BLOSUM62():
    def __init__(self,compareStrA=None,compareStrB=None,compareSqeFile=None,Scorefilepath=None,ScoringMatrix=None):
        self.__Scorefilepath = Scorefilepath
        self.__ScoringMatrix = ScoringMatrix
        self.__compareStrA = compareStrA
        self.__compareStrB = compareStrB
        self.__compareSqeFile = compareSqeFile
        self.__splitChar = None
        self.CompareMatrix = None
        
    def __getScoringMatrix(self):
        if(self.__ScoringMatrix is None):
           ScoreDict = {}
           with open(self.__Scorefilepath,'r',encoding='utf-8') as ftp:
               line = list(ftp.readlines().replace(' ',''))               
               for line in ftp.readlines():
                   listdata = line.replace('  ',' ').spilt(" ")
                   key1 = listdata[0]
                   ScoreDict[key1] = {}
                   valueScore = listdata[1:-1]+listdata[-1]
                   
                   for i in range(0,len(line)):
                      valueScore[key1][line[i]] = int(valueScore[i])
        self.__ScoringMatrix = valueScore
    
    def __getcompareSequence(self):
        if(self.__compareStrA is None and self.__compareStrB is None):
            if(self.__compareSqeFile is None):
                print("please input compare Sequence")
                exit()
            else:
                with open(self.__compareSqeFile,'r') as ftp:
                    Sequence = ftp.read().split(self.__splitChar)
                    self.__compareStrA = Sequence[0]
                    self.__compareStrB = Sequence[1]
                    
    def setSpiltChar(self,splitChar):
        self.__splitChar = splitChar
        
    def __Initialization(self,lensize):
        CompareMatrix = np.zeros(lensize,dtype='int32')
        CompareMatrix[:,0] = 0
        CompareMatrix[0,:] = 0
        
        return CompareMatrix
        
    def __getCompareMatrix(self):
        lensize = [len(self.__compareStrA)+1,len(self.__compareStrB)+1]
        CompareMatrix = self.__Initialization(lensize)
        
        for i in range(1,lensize(0)):
            for j in range(1,lensize(1)):
                listvalue = []
                listvalue.append(CompareMatrix[i-1,j-1] + self.__getScoringMatrix[self.__compareStrA[i-1]][self.__compareStrB[j-1]])
                listvalue.append(CompareMatrix[i,j-1])
                listvalue.append(CompareMatrix[i-1,j])
                listvalue.append(0)
                CompareMatrix[i][j] = max(listvalue)
        self.CompareMatrix = CompareMatrix

    def __Traceback(self):
        pass
        '''
        不想写了，参考https://github.com/wangyang2014/Local-Alignment-Problem/blob/master/local%20alignment%20problem/LAPClass.cpp
        的void LAPCLASS::Traceback(vector<point> &position) 函数的实现！
        '''
