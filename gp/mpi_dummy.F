#define LUPRI_MPI 0
      SUBROUTINE MPI_INIT(IERR)
         IERR = 0
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_INIT called, nothing done'
#ifdef MPI_DEBUG
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_COMM_RANK(I1,I2,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_COMM_RANK called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_COMM_SIZE(I1,I2,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_COMM_SIZE called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_FINALIZE(IERR)
         IERR = 0
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_FINALIZE called, exiting ...'
#ifdef MPI_DEBUG
         CALL QTRACE(LUPRI_MPI)
#endif
         CALL QUIT('dummy MPI_FINALIZE called, exiting ...')
      RETURN
      END
      SUBROUTINE MPI_ABORT()
         CALL QUIT('dummy MPI_ABORT called, exiting ...')
      RETURN
      END
      SUBROUTINE MPI_BCAST(i1,i2,i3,i4,i5,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_BCAST called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_RECV(BUFFER,COUNT,DATATYPE,MPI_ANY_SOURCE,TAG,
     &                    MPI_COMM_WORLD,ISTAT,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_RECV called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_SEND(BUFFER,COUNT,DATATYPE,DEST,TAG,
     &                    MPI_COMM_WORLD,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_SEND called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_REDUCE()
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_REDUCE called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_BARRIER()
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_BARRIER called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_ERROR_CLASS(IERR,IERRCL,IERR2)
        WRITE(LUPRI_MPI,*)
     &      'info: dummy MPI_ERROR_CLASS called, ierr =',IERR
        CALL QUIT('ERROR, dummy MPI_MYFAIL called')
      RETURN
      END
