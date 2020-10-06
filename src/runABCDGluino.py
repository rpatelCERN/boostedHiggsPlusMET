import os
from ROOT import *
hino=[]

#For Gluino
for h in range(0, 16):
	hino.append(1000+h*100)
for h in hino:
	os.system("python QuickDataCardsABCD_Gluino.py %d" %(h))
