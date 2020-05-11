### 
### Programmed by Antoine BALDO, 2017

import os 
import random
from pprint import pprint

# Var number (nV dimension) setpoint
# You can change the dimension's size
nV = 4

# Sample number setpoint
# You can change the Sample number
nS = 200

# Initialisation:
k=1
# Creation of a list dictionnary
x = {}
x[k] = []

# Loop elements (part1)
for i in range(1,(nV+1)):
	x1 = []

	for j in range(1,(nS+1)):
		a = ((float(j)-1)/nS)
		b = ((float(j))/nS)
		listesample = random.uniform(a,b)
		x1.append(listesample)

	# Select a random number nP times between each Sample and for each Var (part2)
	for k in range(1,nS+1):
		listechoice = random.choice(x1)
		x.setdefault(k, []).append(listechoice)
		x1.remove(listechoice)
	#print 'var %d:' % i
	#x2 = [ '%.6f' % elem for elem in x1]
	# print("%.2f" % x1)
	#print x2
#pprint(x)
print(x[1])
print(x[1][2])
