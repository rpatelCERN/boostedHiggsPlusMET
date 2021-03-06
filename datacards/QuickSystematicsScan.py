import sys
import os
mGluinos=[750, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,2000, 2100]
#mGluinos=[750,1200,1600, 2000]
#mGluinos=[2100]
for g in mGluinos:
	for i in range(0,1):
		os.system("python buildCards.py --model=T5HZ --mGo %d --mu 0.0 " %(g))	
		idir="UnblindingT5HZ%d" %(g)
		cardsToAdd=""
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_signal%d.txt " %(idir,i)
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_signal1H%d.txt " %(idir,i)
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_antitagRegion%d.txt " %(idir,i)
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_sidebandRegion%d.txt " %(idir,i)
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_sidebandRegion1H%d.txt " %(idir,i)
		for i in range(0,3):
			cardsToAdd=cardsToAdd+"%s/card_sidebandATRegion%d.txt " %(idir,i)
		os.system("combineCards.py %s > %s.txt " %(cardsToAdd,idir ))
		#os.system("combine -M ProfileLikelihood  --uncapped 1 --significance %s.txt -n %s -m %d" %(idir,idir,g))
		os.system("combine -M Asymptotic %s.txt -n %s --rMin=0.0001" %(idir,idir))	
