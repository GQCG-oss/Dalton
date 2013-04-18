#!/usr/bin/env python

# you can run this script independently on cmake,
# example: get_compiler_version.py gfortran

import sys
import os
import subprocess

version = 'unknown'

if len(sys.argv) == 1:
    # print unknown version and exit otherwise this does not work on non-cxx projects
    print(version)
    sys.exit()

compiler_name = os.path.split(sys.argv[1])[-1]

command_d = {}

# Fotran compilers
command_d['mpif90']       = 'mpif90       --version'
command_d['gfortran']     = 'gfortran     --version'
command_d['gfortran.exe'] = 'gfortran.exe --version'
command_d['gfortran44']   = 'gfortran44   --version'
command_d['f95']          = 'f95          --version'
command_d['ifort']        = 'ifort        --version'
command_d['pgf90']        = 'pgf90        -V'
command_d['xlf']          = 'xlf          -qversion'

# C compilers
command_d['mpicc']        = 'mpicc        --version'
command_d['gcc']          = 'gcc          --version'
command_d['gcc.exe']      = 'gcc.exe      --version'
command_d['gcc44']        = 'gcc44        --version'
command_d['icc']          = 'icc          --version'
command_d['pgcc']         = 'pgcc         -V'
command_d['xlc']          = 'xlc          -qversion'

# C++ compilers
command_d['mpiCC']        = 'mpiCC        --version'
command_d['g++']          = 'g++          --version'
command_d['g++.exe']      = 'g++.exe      --version'
command_d['c++']          = 'c++          --version'
command_d['g++44']        = 'g++44        --version'
command_d['iCC']          = 'iCC          --version'
command_d['icpc']         = 'icpc         --version'
command_d['pgCC']         = 'pgCC         -V'
command_d['xlCC']         = 'xlCC         -qversion'

if sys.version >= '2.4':
    if compiler_name in command_d:
        p = subprocess.Popen(command_d[compiler_name], shell=True, \
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate()[0]
        if out != '':
            if compiler_name != 'pgf90' and compiler_name != 'pgcc' and  compiler_name != 'pgCC':
                # use only first line, not the copyright stuff
                version = out.split('\n')[0]
            else:
                # for Portland compilers, however, use the second line (the first line is empty)
                version = out.split('\n')[1]
print(version)
