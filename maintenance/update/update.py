#!/usr/bin/env python

import sys
import subprocess
import os

# This script runs git pull. If no merges need to be done it will do so silently,
# otherwise it will expose the git back-end. 
# By Ulf Ekstrom 2013.

package_name = "DALTON/LSDALTON"

def run_git(args=[]):
    """Run git and return (returncode, stdout, stderr)"""
    null = open(os.devnull,'r')
    try:
        P = subprocess.Popen(['git']+args, stdin=null, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = P.communicate()
    except OSError:
        sys.stderr.write("""Cannot find `git'. Please install Git from http://git-scm.com/downloads or ask your system administrator.\n""")
        sys.exit(-1)
    return P.returncode,out,err

def reincarnate():
    """Executes the possibly new self with the same arguments as was used to launch update
    Beware the infiloop.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    # The abspath solution works also if the update script was run with python update.py.
    # Normally the script should not be on the path, but potentially this could 
    # run a different update.py than the original version.
    os.execv(os.path.abspath(sys.argv[0]),sys.argv)

def main():
    if len(sys.argv) == 1:
        pass
    elif sys.argv[1] in ['-h','-help','--help']:
        print "Usage: update.py [--force] [--help]"
        print """Updating script for %s. Run without arguments to check for new updates
and apply them. This script should be run from the main %s directory.""" % (package_name, package_name)
        print "Run with --force to throw away local changes during the update."
        print "If you are a Git expert you can directly do a `git pull' instead of using this script."
        sys.exit(0)
    elif sys.argv[1] != '--force':
        sys.stderr.write("Unknown argument `%s', quitting.\n" % sys.argv[1])
        sys.exit(-1)
    res, dum, dum = run_git(['status']) # Check if we are in a git repo
    if (res != 0) or (not os.path.isdir('.git')): # we also check that we are not accidentaly within some other git repo
        sys.stderr.write("ERROR: The update script must be run from inside the %s directory\n" % package_name)
        sys.exit(-1)
    res = 0
    if len(sys.argv) > 1 and sys.argv[1] == '--force':
        res,stdout,stderr = run_git(['checkout','.'])
    if res == 0:
        res,stdout,stderr = run_git(['pull'])
    if res != 0:
        sys.stderr.write("ERROR: The update script failed.\n")
        sys.stderr.write("       If the error is due to your local changes, resolve the situation, then 'git pull' manually.\n")
        sys.stderr.write("       This is the error message:\n\n")
        sys.stderr.write(stdout)
        sys.stderr.write(stderr)
        sys.exit(-1)
    if stdout == "Already up-to-date.\n":
        print "Up-to-date."
    else:
        print "Applied updates to the %s source code." % package_name
        print "Remember to recompile %s after updating." % package_name
        reincarnate()
    
if __name__== "__main__":
    main()
