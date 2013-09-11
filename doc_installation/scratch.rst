

Scratch directory
=================

DALTON and LSDALTON need a scratch directory to write temporary files.
Ideally this should be a fast-access disk.

You should always specify an explicit scratch directory with::

  $ export DALTON_TMPDIR=/full/path/scratch

or by using::

  $ dalton   -t /full/path/scratch [other flags and options]
  $ lsdalton -t /full/path/scratch [other flags and options]

Which overrides ``DALTON_TMPDIR``.
If ``DALTON_TMPDIR`` is neither set nor passed to the run scripts, DALTON
and LSDALTON will search
/global/work/$USER /scratch/$USER /work /scratch /scr /temp /tmp
as candidates for a scratch directory. However you should not let
(LS)DALTON default to those.

Do not point ``DALTON_TMPDIR`` to your home directory!
To prevent loss of data (LS)DALTON 
always appends a directory to ``DALTON_TMPDIR``.
