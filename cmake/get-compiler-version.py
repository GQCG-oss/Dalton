#!/usr/bin/env python

# script to get compiler version

import os
import sys
import string
import subprocess
from optparse import OptionParser, OptionGroup

# define example usage

usage = '''
  Examples:
    ./%prog --compiler=ifort
    '''

# initialize parser

parser = OptionParser(usage)

# define options

parser.add_option('--compiler',
                  type='string',
                  action='store',
                  dest='compiler',
                  default=None,
                  help='set the compiler [default: %default]',
                  metavar='STRING')

# process input

(options, args) = parser.parse_args()

# define commands

command_d = {}

command_d['gfortran'] = 'gfortran --version'
command_d['g95']      = 'g95      --version'
command_d['ifort']    = 'ifort    --version'
command_d['pgf90']    = 'pgf90    -V'
command_d['xlf']      = 'xlf      -V'

command_d['gcc']      = 'gcc      --version'
command_d['icc']      = 'icc      --version'
command_d['pgcc']     = 'pgcc     -V'
command_d['xlc']      = 'xlc      -V'

command_d['g++']      = 'g++      --version'
command_d['iCC']      = 'iCC      --version'
command_d['pgCC']     = 'pgCC     -V'
command_d['xlCC']     = 'xlCC     -V'

def main():
        if len(sys.argv) == 1:
                # user has given no arguments: print help and exit
                print parser.format_help().strip()
                sys.exit()

        compiler_name = options.compiler.split('/')[-1]

        if not compiler_name in command_d:
            print 'unknown version'
        else:
            p = subprocess.Popen(command_d[compiler_name], shell=True, \
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = p.communicate()[0]
            if out == '':
                print 'unknown version'
            else:
                # use only first line, not the copyright stuff
                print out.split('\n')[0]

if __name__ == '__main__':
        main()
