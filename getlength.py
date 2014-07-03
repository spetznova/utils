import sys

file = sys.argv[1]

f = open(file,'r')
list = f.read().splitlines()
f.close()

print len(list)
