import os
from ROOT import *
hino=[]
LSP=[]

#For higgsino
for h in range(0, 34):
	hino.append(175+h*25)
hino.append(1100)
hino.append(1200)
for h in hino:
	# os.system("../bin/ALPHABET 1 MC2016 %d 1 0" %(h))
	os.system("python QuickDataCardsABCDNorm_Higgsino.py %d 1" %(h))
