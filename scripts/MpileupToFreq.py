#!/usr/bin/env python

import sys
import re
inFile = open(sys.argv[1],'r')

print 'contig\tbp\tA\tC\tG\tT\tdel\tins\tinserted\tambiguous'

digitre=re.compile("(\d+).*")

for line in inFile:
        data = line.strip().split('\t')
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()
        contig = data[0]
        types = {'A':0,'C':0,'G':0,'T':0,'-':0,'+':[],'X':[]}

        i = 0

        while i < len(bases):
                base = bases[i]

                if base == '^' or base == '$':
                        i += 1
                elif base == '-':
                        i += 1
                elif base == '*':
                        types['-'] += 1
                elif base == '+':
                        i += 1
                        dm = digitre.match(bases[i:])
                        addNum = int(dm.groups()[0])
                        addSeq = ''
                        for a in range(addNum):
                                i += 1
                                addSeq += bases[i]

                        types['+'].append(addSeq)
                elif base == '.' or base == ',':
                        types[ref] += 1
                else:
                        if types.has_key(base):
                                types[base] += 1
                        else:
                                types['X'].append(base)

                i += 1

        adds = '.'
        if len(types['+']) > 0:
                adds = ','.join(types['+'])

        amb = '.'
        if len(types['X']) > 0:
                amb = ','.join(types['X'])

        out = [contig,int(bp),types['A'],types['C'],types['G'],types['T'],types['-'],len(types['+'])]
        print '\t'.join([str(x) for x in out])

