import time
import random
import matplotlib.pyplot as plt
import math
import sys
import numpy as np
#Configurable parameters
N = 100
inputFileName = "test0.txt"
outputFileName = "out0.txt"

primeNumbers = [393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319]
currPrime = 0
p = primeNumbers[currPrime]
a0 = 7
a1 = 3

# Take name of input file and output file from user.
if(len(sys.argv) == 3):
	inputFileName = sys.argv[1]
	outputFileName = sys.argv[2]
else:
	print("Usage: python linearProb.py <inputFileName> <outputFileName>")
	exit()

Table = []
for i in range(0,N):
	Table.append(-1)

#open files 
inputFile = open(inputFileName, "r")
outputFile = open(outputFileName, "w")

toPrint = "Hash Function: ( ("+ str(a0) + "+"+ str(a1)+ "x) mod "+str(p)+") mod "+ str(N)+"\n"
outputFile.write(toPrint)
print(toPrint)
outputFile.write("i/l/d\tElement\tTime\t1/(1-Load-Factor)\n")

# Maintains the number of elements currently present in the tables
currentNumberOfElemns = 0
# stores current load factor = currentNumberOfElemns/ N
loadFactor = 0
loadX = [] #stores 1/1-load factor for the insertion done
timeY = [] #stores time taken for the insertion done in micro seconds

###########################
###########################
# find the best fit line given two lists x and y
def best_fit(X, Y):

    avgx = sum(X)/len(X)
    avgy = sum(Y)/len(Y)
    n = len(X) 

    num = sum([xi*yi for xi,yi in zip(X, Y)]) - n * avgx * avgy
    den = sum([xi**2 for xi in X]) - n * avgx**2

    m = num / den
    c = avgy - m * avgx

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(c, m))

    return c, m

###########################
###########################
#plots scattered and log-log graph based on the pair (time for insert, load factor at that time)
def plotGraph(x, y):
	w = 15
	h = 12
	d = 70
	plt.figure(figsize=(w, h), dpi=d)

	#plt.plot(x, y, 'ro')
	
	#plt.xlim(0,2.5)
	#plt.ylim(0,20) 

	#plt.scatter(x, y)

	# solution
	#a, b = best_fit(x, y)
	#yfit = [a + b * xi for xi in x]
	#plt.plot(x, yfit)
	
	#plt.show()

	plt.subplot(211)
	plt.plot(x, y, 'o')
	
	
	plt.scatter(x, y)

	# solution
	a, b = best_fit(x, y)
	yfit = [a + b * xi for xi in x]
	plt.plot(x, yfit, '--',label='m=%.2f' % b)
	plt.legend()
	
	plt.title('Scattered plot')
	plt.xlabel('load factor')
	plt.ylabel('Time per Insert Operation (micro-seconds)')
	
	# Numbers between 10 and 100 equally spaced in log scale
	x = np.logspace(1, 2)

  	# Assume ys have a polynomial relationship with xs
	y = x ** 6

  	# Add noise to ys, so it's not a straight line
	y *= 1 + np.random.random(len(x))

  	# Log-log scale
	plt.subplot(212)
	plt.title('Log-log scale plot')
	plt.loglog(x, y, marker='o', linestyle='None')

  	# Regression
	m, b = np.polyfit(np.log10(x), np.log10(y), 1)

  	# Only need two points to define a line
	line_xs = np.array([1, 2])
	line_ys = m * line_xs + b

  	# In order to superimpose plots we need to keep using log-log scale, so we
  	# compensate by exponentiation.
	plt.loglog(10 ** (line_xs), 10 ** (line_ys), '--', label='m=%.2f' % m)
	plt.legend()
	plt.xlabel('load factor')
	plt.ylabel('Time per Insert Operation (micro-seconds)')
	plt.show()
	#plt.savefig("out_25000.png")

###########################
###########################
# resizes the hash table to size and rehashs elements 
def resizeHashTable(size):
	global N, currPrime, currentNumberOfElemns,Table

	N = size

	toPrint = "Resizing table to length " + str(N) + "\n"
	outputFile.write(toPrint)
	print(toPrint)

	#initalize the temporary table with elements from Table which are not -1 and -2
	tempTable = []
	for i in range(0,len(Table)):
		if(Table[i] != -1 and Table[i] != -2):	
			tempTable.append(Table[i])

	
	# clearning current elements present in Table
	for i in range(0,len(Table)):
		Table[i] = -1

	#if( N > len(Table)):
		#print("Adding")
		#i = len(Table)
		#while i < N:
			#print(i)
			#Table.append(-1)
			#i += 1

	# If size is to be increased than add N/2 elements in Table
	if( N > len(Table)):
		for i in range(0,math.ceil(N/2)):
			Table.append(-1)
	else: # else remove last N/2 elements from in Table
		Table = Table[0:N]

	#print(N, "New Size ", len(Table),"len")

	# change hash function if necessary
	#while(size > primeNumbers[currPrime]):
		#currPrime += 1
		#p = primeNumbers[currPrime]

	currentNumberOfElemns = 0

	# insert in new table the elements from temporary table
	for i in range(0, len(tempTable)):
		insertIntoHashTable(tempTable[i])

###########################
###########################
# Given a number/element it computes its index based on hash function and returns the value
def calculateIndex(number):
	global ao, a1, p, N
	temp = ((a0 + a1*int(number)) % p) % N
	return temp 
###########################
###########################
# places element at the index
# and notes down the time and load factor in output file 
def insertElementAndNoteTime(index1, value, currentNumberOfElemns, startTimeStamp):
	global Table, loadFactor, timeY, loadX, N

	Table[index1] = value
	print("Inserting ", index1, "-", value, " ==> Inserted")
	loadFactor = currentNumberOfElemns/N

	#insert difference ( start time, end time ) and loadFactor in output file.
	#timeY.append(round((time.time()- startTimeStamp)*pow(10,6),2))
	#loadX.append(round(loadFactor,2))
	loadX.append(round(1+loadFactor,2))
	timeY.append(round((time.time()- startTimeStamp)*pow(10,6),2))
	toPrint = "i\t" + value + "\t" + str(round((time.time()- startTimeStamp)*pow(10,6),2))  + "  \t"+ str(round(1/(1-loadFactor),2)) + "\n"
	outputFile.write(toPrint)

###########################
###########################
# find the appropriate index at which an element can be inserted based on hash function index
#in case element is not found at expected location does a linear search to find it
# resizes table if occupancy is more than 75%
def insertIntoHashTable(element):
	global currentNumberOfElemns, N
	startTimeStamp = time.time()

	#If number of elements increase more than 75% than resize table to twice its size
	if(currentNumberOfElemns > (3/4)*N):
		resizeHashTable(2*N)

	index = calculateIndex(element)

	print("index -", index, " element -", element )
	if Table[index] == -1 or Table[index] == -2 :	
		currentNumberOfElemns += 1			
		insertElementAndNoteTime(index, element, currentNumberOfElemns, startTimeStamp)		
	else:								
		i = (index + 1) % N
		while(i != index):					
			if Table[i] == -1 or Table[i] == -2 :		
				currentNumberOfElemns += 1
				insertElementAndNoteTime(i, element, currentNumberOfElemns, startTimeStamp)
				break
			
			i = i+1
			i = i%N
			#print(i, "-", N, "this")
###########################
###########################
# looks up the element in the possible location and if not present does linear search till either -1 
# is encountered or we come back to the prev location from which we started.
def lookUpInHashTable(element):
	global currentNumberOfElemns, N

	startTimeStamp = time.time()

	index = calculateIndex(element)
	if Table[index] == -1:
		print("Looking up ", element, " ==> NoMatch")
	elif Table[index] == command[1]:
		print("Looking up ", element, " ==> Match")
	else:
		flagFound = False
		i = (index + 1) % N
		while(i != index and Table[i] != -1):					#loop till end of Table
			if Table[i] == element:	
				flagFound = True		#if free spot found till end of Table then insert and exit loop
				print("Looking up ", element, " ==> Match")
				break
			i += 1
			i %= N

		if not flagFound:
			print("Looking up ", element, " ==> NoMatch")

###########################
###########################
# searches the element just like lookup and if found then puts -2 at its place also rehashes when the
#occupancy of table is less than 25%
def deleteFromHashTable(element):
	global N, currentNumberOfElemns, Table, loadFactor

	index = calculateIndex(element)

	if(Table[index] == element):
		Table[index] = -2
		currentNumberOfElemns -= 1
		loadFactor = currentNumberOfElemns/N
	else:
		i = (index+1) % N
		while(i != index):
			if(Table[i] == -1):
				break
			elif(Table[i] == element):
				Table[i] = -2
				currentNumberOfElemns -= 1
			i += 1
			i %= N

	#If number of elements are less than 25% of the size of tabel than resize table to half its size
	if(4*currentNumberOfElemns < N and N > 10):
		resizeHashTable(int(N/2))
###########################
###########################	
ST = time.time()		
# start reading line from input file one by one. Perform insert with i x and lookup with l x.
line = inputFile.readline()

while(line):
	#print(line)
	command = line.split()
	#insert element
	if(command[0] == "i"): 
		insertIntoHashTable(command[1])
			
	#lookup element
	elif(command[0] == "l"):
		lookUpInHashTable(command[1])
		

	#delete element
	elif(command[0] == "d"):
		deleteFromHashTable(command[1])

	line = inputFile.readline()

print("######TABLE######")
for i in range(0, len(Table)):
	print(i, "-->", Table[i])

ET = time.time()

plotGraph(loadX,timeY)


toPrint = "Total Time in milli-seconds: " + str(round((ET-ST)*pow(10,3),2)) + "\n"
outputFile.write(toPrint)
print(toPrint)

outputFile.close()
inputFile.close()

