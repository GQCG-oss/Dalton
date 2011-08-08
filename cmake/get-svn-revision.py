#!/usr/bin/env python

import re
import string
import subprocess

p = subprocess.Popen('svn info', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()

if stderr:
    print 'unknown'
else: # process first three lines of the standard output
    for line in stdout.splitlines():
        if 'Revision:' in line:
            commit = line.split()[-1]
        if 'Last Changed Author:' in line:
            author = line.split(':')[-1]
        if 'Last Changed Date:' in line:
            date = line.split('Last Changed Date:')[-1]
    print '%s   %s  %s' % (commit, author, date)
