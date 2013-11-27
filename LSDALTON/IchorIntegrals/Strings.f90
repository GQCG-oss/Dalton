MODULE stringsMODULE

character(len=132) :: STRING 
integer :: iString

INTERFACE AddToString
  MODULE PROCEDURE AddIntegerToString,AddCharacterToString
END INTERFACE ADDTOSTRING

CONTAINS
subroutine initString(N)
integer :: N,I
iSTRING = 1
do I=1,N
   STRING(iSTRING:iSTRING) = ' '
   iSTRING = iSTRING+1
enddo
end subroutine initString

subroutine writeString(lupri)
integer :: lupri
IF(lupri.EQ.-1)THEN
   WRITE(*,'(A)') STRING(1:iSTRING-1)
ELSE
   WRITE(lupri,'(A)') STRING(1:iSTRING-1)
ENDIF
end subroutine writeString

subroutine AddIntegerToString(I)
integer :: I

if(I.LT.10)THEN
   WRITE(STRING(iSTRING:iSTRING),'(I1)')I
   iSTRING=iSTRING+1
elseif(I.LT.100)THEN
   WRITE(STRING(iSTRING:iSTRING+1),'(I2)')I
   iSTRING=iSTRING+2
elseif(I.LT.1000)THEN
   WRITE(STRING(iSTRING:iSTRING+2),'(I3)')I
   iSTRING=iSTRING+3
elseif(I.LT.10000)THEN
   WRITE(STRING(iSTRING:iSTRING+3),'(I4)')I
   iSTRING=iSTRING+4
elseif(I.LT.100000)THEN
   WRITE(STRING(iSTRING:iSTRING+4),'(I5)')I
   iSTRING=iSTRING+5
else
   stop 'ERROR AddIntegerToString'
endif
end subroutine AddIntegerToString

subroutine AddCharacterToString(B)
character*(*) :: B
character(len=4) :: FORMAT1
character(len=5) :: FORMAT2
integer :: L
L = LEN(B)
IF(L.LT.10)THEN
   WRITE(FORMAT1,'(A2,I1,A1)')'(A',L,')'
   WRITE(STRING(iSTRING:iSTRING+L-1),FORMAT1) B 
ELSE
   WRITE(FORMAT2,'(A2,I2,A1)')'(A',L,')'
   WRITE(STRING(iSTRING:iSTRING+L-1),FORMAT2) B 
ENDIF
iSTRING = iSTRING+L
end subroutine AddCharacterToString

end MODULE STRINGSMODULE

