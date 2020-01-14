#!/usr/bin/env python3
import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

INSTALL_DIR = pathlib.Path(__file__).resolve().parent
DALEXE = INSTALL_DIR / 'dalton.x'
BASDIR = INSTALL_DIR / 'basis'
os.environ['BASDIR'] = str(BASDIR)


def main():

    args = parse_arguments()

    with tempfile.TemporaryDirectory() as tmp:
        upload_files(args, tmp)
        p = subprocess.run(args.exe, cwd=tmp)
        download_files(args, tmp)

    return p.returncode


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+', help='dalton input')
    parser.add_argument('-f', help='restart file')
    parser.add_argument('-d', action='store_true', help='unused')
    parser.add_argument('-ext', default='out', help='output extension')
    parser.add_argument('-noarch', action='store_true')
    parser.add_argument('-nobackup', action='store_true')
    parser.add_argument('-get', default="")
    parser.add_argument('-put', default="")
    parser.add_argument('-N', type=int, default="1")
    parser.add_argument('-exe', default=[DALEXE])

    args = parser.parse_args()

    if args.N > 1:
        args.exe = ['mpirun', '-np', f'{args.N}'] + args.exe
    return args


def upload_files(args, tmp):

    dal, mol, pot, out = process_args(args.inputs)
    if args.f:
        with tarfile.open(args.f + '.tar.gz', mode='r:gz') as tgz:
            tgz.extractall(path=tmp)

    shutil.copy(dal + '.dal', os.path.join(tmp, 'DALTON.INP'))
    try:
        shutil.copy(mol + '.mol', os.path.join(tmp, 'MOLECULE.INP'))
    except FileNotFoundError:
        pass

    if pot:
        try:
            shutil.copy(pot + '.pot', os.path.join(tmp, 'POTENTIAL.INP'))
        except FileNotFoundError:
            pass

    for f in args.put.split():
        shutil.copy(f, os.path.join(tmp, f))


def download_files(args, tmp):

    dal, mol, pot, out = process_args(args.inputs)
    for f in args.get.split():
        shutil.copy(os.path.join(tmp, f), f"{out}.{f}")

    with tarfile.open(out + ".tar.gz", mode="w:gz") as tgz:
        tgz.add(
            tmp,
            arcname='.',
            filter=lambda f:
                None if re.match(r'.*AO.*|.*\/$', f.name) else f
        )

    shutil.copy(os.path.join(tmp, 'DALTON.OUT'), f"{out}.{args.ext}")


def process_args(inputs):

    assert inputs

    mol = None
    pot = None
    if len(inputs) == 1:
        dal, = inputs
        mol = dal
        pot = dal
        out = dal
    elif len(inputs) == 2:
        dal, mol = inputs
        if dal == mol:
            out = dal
        else:
            out = f'{dal}_{mol}'
    else:
        dal, mol, pot = inputs
        if dal == mol and mol == pot:
            out = dal
        else:
            out = f'{dal}_{mol}_{pot}'
    return dal, mol, pot, out


if __name__ == "__main__":
    sys.exit(main())
