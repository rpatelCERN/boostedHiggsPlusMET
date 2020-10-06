import os
from ROOT import *
hino=[]
LSP=[]

#For higgsino
for h in range(0, 52):
	hino.append(175+h*25)
for h in hino:
	os.system("python QuickDataCardsABCD.py %d 1" %(h))
