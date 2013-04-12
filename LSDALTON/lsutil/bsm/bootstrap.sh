#! /bin/sh
# bootstrap file to be used when autogen.sh fails.
# consider using --foreign flat to automake
echo "Running aclocal..."
aclocal || exit 1
#echo "Running libtoolize..."
#libtoolize --force || exit 1
#echo "Running autoheader..."
#autoheader || exit 1
echo "Running autoconf..."
autoconf || exit 1
echo "Running automake..."
automake --foreign --add-missing --copy || exit 1
# automake-1.8 does not copy all the required files properly.
test -f config.sub || {
  echo "Your automake is buggy and does not copy
config.sub and config.guess files properly.
I will try to work around it but the recommended solution 
is to upgrade to automake-1.9.1 or more recent."
  libtoolize -c > /dev/null 2>&1 
}
echo "Running configure $* ..."
./configure "$@"
