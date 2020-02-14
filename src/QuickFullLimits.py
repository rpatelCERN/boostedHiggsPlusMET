import os
from ROOT import *
hino=[]

for h in range(0, 34):
	hino.append(175+h*25)
for h in hino:
	os.system("python QuickDataCardsABCD.py %d" %(h))
