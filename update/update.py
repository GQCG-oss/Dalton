#!/usr/bin/env python
import sys, urllib2, subprocess, os

# Update script 2012 by Ulf Ekstrom

# The url base is read from this file
url_base = 'update/update-url'

# The current version is read from this file
version_file = 'VERSION'
# VERSION contains a number in format X.Y (example: 12.0)
# The script will try to download and apply url_base/X/Y

def full_url():
    # read url
    try:
        f = open(url_base, 'r')
        url = f.read().strip()
    except:
        sys.stderr.write("Cannot find the file `%s', running update from wrong directory?\n" % url_base)
        sys.exit(-1)
    # read version
    try:
        f = open(version_file, 'r')
        version = f.read().strip()
    except:
        sys.stderr.write("Cannot find the file `%s', running update from wrong directory?\n" % version_file)
        sys.exit(-1)
    version_major = int(version.split('.')[0])
    version_patch = int(version.split('.')[1])
    # look for version_patch+1
    version_patch += 1
    new_version = '%s.%s' % (version_major, version_patch)
    # construct full url
    url += '%s/%s' % (version_major, version_patch)
    return new_version, url


def get_patch(new_version, url):
    try:
        u = urllib2.urlopen(url,timeout=10)
    except urllib2.HTTPError:
        print "No updates found on server."
        return None
    except urllib2.URLError:
        sys.stderr.write("Error downloading patch from %s, cannot connect to server.\n" % url)
        return None
    print "Downloading patch %s:" % new_version,
    try:
        patch = u.read()
    except:
        sys.stderr.write("Error downloading patch from %s, cannot connect to server.\n" % url)
        return None
    print "OK"
    return patch


def has_patchprogram():
    with open(os.devnull, "w") as fnull:
        res = subprocess.call(['patch','-v'],stdout = fnull, stderr = fnull)
    if res == 0:
        return True
    else:
        return False


def test_or_apply_patch(patch,test=True):
    with open(os.devnull, "w") as fnull:
        args = ['patch','-p1']
        if test:
            args.append('--dry-run')
        p = subprocess.Popen(args,stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate(patch)
        res = p.returncode
    if res == 0:
        return True
    else:
        f = open('update-error.log','w')
        f.write(out[0])
        f.write(out[1])
        return False


def main():
    if len(sys.argv) > 1:
        print "Usage: update.py"
        print "If run with no arguments this script attempts to download and apply updates to DIRAC."
        print "The updating procedure relies on the UNIX `patch' program, and may fail if you have"
        print "made local modifications to the DIRAC source code."
        sys.exit(-1)
    new_version, url = full_url()
    p = get_patch(new_version, url)
    if p:
        if not has_patchprogram():
            sys.stderr.write("The `patch' program was not found, cannot apply updates.\n")
            sys.exit(-1)
        if test_or_apply_patch(p,test=True):
            print "Verified patch, now applying",
            res = test_or_apply_patch(p,test=False)
            if res:
                print "Update applied successfully."
                new_version, new_url = full_url()
                if new_url == url:
                    sys.stderr.write("After update: Idiotic patch, updated version same as the old. Please report this to the developers.\n")
                    sys.exit(-1)
            else:
                sys.stderr.write("Something went wrong applying the patch even after verification, see update-error.log.\n")
                sys.stderr.write("Please extract the official distribution archive and update directly\nfrom that.\n")
                sys.stderr.write("The source tree may be left in an inconsistent and non-compiling state.\n")
                sys.exit(-1)
        else:
            sys.stderr.write("Patch cannot be applied, please find error log in update-error.log.\n")
            sys.stderr.write("This may be due to your modifications of DIRAC or a developer error.\n")
            sys.stderr.write("Please extract the official distribution archive and update directly\nfrom that.\n")
            sys.stderr.write("The error was detected before any files where modified.\n")
            sys.exit(-1)
    else:
        sys.exit(0)


print 'DIRAC update script, attempting to locate and download updates ...'
print
while True:
    main()
    print "Checking for more updates ..."
