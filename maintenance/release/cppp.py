#! /usr/bin/env python

# (c) Ulf Ekstrom

import sys, re, os

class FriendlyFile:
    def __init__(self, path, mode):
        self.f = open(path,mode)
        self.line_nr = 0
        self.path = path
    def readline(self):
        l = self.f.readline()
        if len(l) > 0:
            self.line_nr += 1
        return l
    def locstring(self):
        return "%s:%i" % (self.path,self.line_nr)


re_ifany = re.compile('^\s*#\s*if')
re_ifdef = re.compile('^\s*#\s*ifdef\s+(\w*)')
re_ifndef = re.compile('^\s*#\s*ifndef\s+(\w*)')
re_if = re.compile('^\s*#\s*if\s+(\S*)')
re_else = re.compile('^\s*#\s*else')
re_endif = re.compile('^#\s*endif')
re_elif = re.compile('^\s*#\s*elif\s+(\S*)')
re_else_or_endif = re.compile('^\s*#\s*e[lsnd][lsnd]')
re_eof = re.compile('\b')

defined = []
undefined = []

def skip_until(infile, expr):
    """Skip lines until expr matches, return the matching line"""
    line = infile.readline()
    while line != '' and not expr.match(line):
        if re_ifany.match(line):
            skip_until(infile,re_endif)
        line = infile.readline()
    return line

def process_chunk(infile, dst, until_expr = re_eof):
    """Read (recursively) a chunk assumed to start just below a #ifdef or
    similar. Return at the matching #endif or #else, but do not print that
    line, just return it."""
    line = infile.readline()
    while line != '':
        if until_expr.match(line):
            return line
        elif re_ifdef.match(line):
            m = re_ifdef.match(line)
            name = m.group(1)
            if name in defined:
                line = process_chunk(infile,dst,re_else_or_endif)
                if re_else.match(line):
                    skip_until(infile,re_endif)
            elif name in undefined:
                line = skip_until(infile,re_else_or_endif)
                if re_else.match(line):
                    process_chunk(infile,dst,re_endif)
            else:
                dst.write(line)
        elif re_ifndef.match(line):
            m = re_ifndef.match(line)
            name = m.group(1)
            if name in undefined:
                line = process_chunk(infile,dst,re_else_or_endif)
                if re_else.match(line):
                    skip_until(infile,re_endif)
            elif name in defined:
                line = skip_until(infile,re_else_or_endif)
                if re_else.match(line):
                    process_chunk(infile,dst,re_endif)
            else:
                dst.write(line)
        else:
            dst.write(line)
        line = infile.readline()
    return line

usage = """Usage: cppp [-DMACRO ..] [-UMACRO ...] FILE [FILE ..]
For each FILE perform a partial C preprocessing of #ifdef and #ifndef
directives, completely removing sections that cannot be reached with
the given macros defined (-D) or undefined (-U).

For example, the file

#ifdef A
So A was defined
#ifdef C
#ifndef B
Only if B is not defined
#endif
#endif
#else
A was not defined
#endif

With the arguments -DA -UB will be transformed to

So A was defined
#ifdef C
Only if B is not defined
#endif


NOTE: The program does not process #if statements, although they are reproduced
correctly in the preprocessed files. You may have to manually convert such statements
to chained #ifdef's.

NOTE2: The files are modified IN PLACE! Run this on a copy of the source tree.

By U. Ekstrom, <uekstrom@gmail.com> 2008

"""


def main():
    filenames = []
    target = filenames
    if len(sys.argv) == 1:
        print usage
        return
    for arg in sys.argv[1:]:
        arg = arg.strip()
        if arg.find('-D') == 0:
            arg = arg[2:]
            target = defined
        elif arg.find('-U') == 0:
            arg = arg[2:]
            target = undefined
        if len(arg) > 0:
            target.append(arg)
            if target in [defined, undefined]:
                if '=' in arg:
                    print 'Error: Giving a macro a value (%s) is not supported. Quitting!' % arg
                    return -1
        target = filenames
    for fn in filenames:
        f = FriendlyFile(fn,'r')
        fout = open(fn+'.cppp~','w')
        process_chunk(f,fout)
        fout.close()
        os.rename(fn+'.cppp~',fn)

main()
