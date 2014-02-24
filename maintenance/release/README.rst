

General
-------

Normally you want to only edit delete_flags and/or
paths_to_exclude. All the other files in this directory
normally do not need any editing (unless they need bugfixes).


delete_flags
------------

if you add -URABOOF to delete_flags, then "make release"
will filter the following code:

#ifdef RABOOF
    remove this 1
#endif
    keep this 1
#ifdef RABOOF
    remove this 2
#else
    keep this 2
#endif

and produce:

    keep this 1
    keep this 2


paths_to_exclude
----------------

List all source files and directories to exlude in
paths_to_exclude, each entry on a new line.

These entries could look like this:

DALTON/dft/secret_new_functionality.F
DALTON/gp/something_else_that_i_dont_want_to_share.F90
