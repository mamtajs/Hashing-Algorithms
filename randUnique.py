import random
import sys

if len(sys.argv) != 3:
	print("Usage: python randUnique.py <outputFileName> <n>")
	exit()
else:
	fileName = sys.argv[1]
	n = sys.argv[2]

inputFile = open(fileName, "w");
list = random.sample(range(1, 50000), int(n))

for elem in list:
	rand = random.randint(0,9)
	if( rand <= 5):
		inputFile.write("i " + str(elem) + "\n");
	elif(rand <= 7):
		inputFile.write("l " + str(elem) + "\n");
	elif(rand <= 9):
		inputFile.write("d " + str(elem) + "\n")
	
#print(list)