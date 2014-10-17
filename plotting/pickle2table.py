import sys
import pickle

if len(sys.argv)<4:
	print 'Usage:\tpython pickle2table.py picklefile dataOutputFilename weightsOutputFilename'
	sys.exit()

f = open(sys.argv[1])

results = pickle.load(f)

data = results[0][-1]

weights = results[1][-1]

dataf = open(sys.argv[2],'w')

weightsf = open(sys.argv[3],'w')

for x in data:
	for y in x:
		dataf.write(str(y)+'\t')
	dataf.write('\n')

for x in weights:
	weightsf.write(str(x)+'\n')