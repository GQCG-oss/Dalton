      SUBROUTINE MPI_INIT(IERR)
         IERR = 0
         WRITE (0,*) 'info: MPI_INIT called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_COMM_RANK(I1,I2,IERR)
         IERR = 0
         WRITE (0,*) 'info: MPI_COMM_RANK called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_COMM_SIZE(I1,I2,IERR)
         IERR = 0
         WRITE (0,*) 'info: MPI_COMM_SIZE called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_FINALIZE(IERR)
         IERR = 0
         WRITE (0,*) 'info: MPI_FINALIZE called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_ABORT()
         WRITE (0,*) 'info: MPI_ABORT called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_BCAST(i1,i2,i3,i4,i5,IERR)
         IERR = 0
         WRITE (0,*) 'info: MPI_BCAST called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_RECV()
         WRITE (0,*) 'info: MPI_RECV called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_SEND()
         WRITE (0,*) 'info: MPI_SEND called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_REDUCE()
         WRITE (0,*) 'info: MPI_REDUCE called, nothing done'
      RETURN
      END
      SUBROUTINE MPI_BARRIER()
         WRITE (0,*) 'info: MPI_BARRIER called, nothing done'
      RETURN
      END
