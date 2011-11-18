#!/usr/bin/env python

import subprocess, sys

binary = sys.argv[1]

command = "size %s 2> /dev/null | grep dalton.x" % binary
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
output = process.stdout.read().strip()

total_size = float(output.split()[2])/(1024*1024)
print('total static size of %s: %10.1f MB' % (binary, total_size))
