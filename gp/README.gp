directory gp: G(eneral) P(urpose) routines

Programming advice:

  We define a module as files in a specific sub-directory.

  All general purpose subroutines and functions which
  - are not dependent on variables in common blocks from any other module directory
    (they can have local "private" common blocks only used inside gp/ )
  - are callable from any module
  SHOULD be stored in this directory.

  In this way we can, among other benefits, hope to avoid duplication of code
  for generic tasks.

  Exceptions: Routines in the public domain should be stored in
  pdpack/ instead of in this directory. The routines in this
  directory are covered by the licence agreement.

  -- Hans Joergen Aa. Jensen (last revision 16. Sep. 2006)

