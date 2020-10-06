import os
from ROOT import *
#hino=[]
#LSP=[]

f = open("higgsino2DFileNames.txt", "r")

#For higgsino, 2D only!
for line in f:
    x = line.split("_")
    hino = int(x[5])
    LSP = int(x[6])
#    os.system("../bin/ALPHABET 1 MC2016 0 %d %d" %(hino,LSP))
#    os.system("../bin/ALPHABET 1 MC2017 0 %d %d" %(hino,LSP))
#    os.system("../bin/ALPHABET 1 MC2018 0 %d %d" %(hino,LSP))
    os.system("python QuickDataCardsABCDNorm_Higgsino.py %d %d" %(hino,LSP))
#os.system("../bin/ALPHABET 1 MC2016 0 400 50")
f.close()

