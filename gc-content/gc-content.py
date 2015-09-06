# Copyright, 2014-2015, Ivan Minevskiy

import collections

size=200
out = []
fin = open('input.txt','r')
list = collections.deque()
for i in range(1,size):
    list.append(fin.read(1))
out.append(float(list.count('C') + list.count('G'))/size)
while True:
    char = fin.read(1)
    if char == "":
        break
    list.popleft()
    list.append(char)
    out.append(float(list.count('C') + list.count('G'))/size)
fin.closed

fout = open('output.txt','w')
for n in out:
    fout.write("%s\n" % n)
fout.closed