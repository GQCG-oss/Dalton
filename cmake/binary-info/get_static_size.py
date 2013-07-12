#!/usr/bin/env python

# get static size of binary
# written by: Radovan Bast <bast@kth.se>
#             Miroslav Ilias <miroslav.ilias@umb.sk>

import subprocess
import sys
import os

nr_top_consumers = 20

if len(sys.argv) == 1:
    print('usage: %s binary' % sys.argv[0])
    sys.exit(1)
else:
    binary = sys.argv[1]

if not os.path.isfile(binary):
    print('binary %s does not exist' % binary)
    sys.exit(1)

command = "size %s 2> /dev/null | grep %s" % (binary, binary)
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
output = process.stdout.read().strip()

print('total static size of %s:' % binary)

total_size = float(output.split()[2])/(1024*1024)
print('%10.2f MB' % total_size)

command = "readelf -s %s 2> /dev/null || echo 'cannot be determined, readelf not present'" % binary
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

d = {}
for line in process.stdout.readlines():
    if 'OBJECT' in line:
        name = line.split()[-1]
        addr = line.split()[1]
        size = int(line.split()[2], 16)
        d[(name, addr)] = float(size)/(1024*1024)

l = []
for entry in d:
    name, addr = entry
    size = d[entry]
    l.append([size, name])

l.sort(reverse=True)

print('\ntop %i objects:' % nr_top_consumers)

mb_sum = 0.0
for i in range(nr_top_consumers):
    size, name = l[i]
    mb_sum += size
    print('%10.2f MB %20s      (sum: %.2f MB)' % (size, name, mb_sum))
