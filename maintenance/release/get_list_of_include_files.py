#!/usr/bin/env python

import sys
import glob
import os
import string

if len(sys.argv) != 3:
   sys.exit('usage: %s /sourcedir /builddir' % sys.argv[0])

sourcedir = sys.argv[1]
builddir  = sys.argv[2]

include_files = ''
for root, directories, files in os.walk(sourcedir):
    for directory in directories:
        full_directory = os.path.join(root, directory)
        if not builddir in full_directory:
            for include_file in glob.glob('%s/*.h' % full_directory):
                include_files += '    %s\n' % string.replace(include_file, sourcedir+'/', '')
            for include_file in glob.glob('%s/*.inc' % full_directory):
                include_files += '    %s\n' % string.replace(include_file, sourcedir+'/', '')

s  = '# file generated do not edit; this is not important for compilation, only for the "make release" target\n'
s += '\nset(INCLUDE_FILES\n%s    )' % include_files

f = open(os.path.join(builddir, 'IncludeFiles.cmake'), 'w')
f.write(s)
f.close()
