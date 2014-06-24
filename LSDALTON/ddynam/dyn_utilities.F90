Module dyn_util
Use precision
Use memory_handling
Use ls_util
Use files
Contains
!===================
!   Move_Molecule
!===================
Subroutine Move_Molecule(NAtoms, Coordinates, MoveVec)
!
! Translate molecule
!
  Implicit None
  Real(realk), Dimension(3*NAtoms) :: Coordinates
  Real(realk), Dimension(3)        :: MoveVec
  Integer                               :: NAtoms, I
!
! Adds specified vector to all the atomic coordinates
!
  Do I = 1, NAtoms
    Coordinates((I-1)*3+1) = Coordinates((I-1)*3+1) + MoveVec(1)
    Coordinates((I-1)*3+2) = Coordinates((I-1)*3+2) + MoveVec(2)
    Coordinates((I-1)*3+3) = Coordinates((I-1)*3+3) + MoveVec(3)
  End Do
End Subroutine Move_Molecule

!====================
!   Center_of_Mass
!====================
Subroutine Center_of_Mass(NAtoms, Coordinates, AtomMass, CMVec)
!
! Determine center of mass for the system
!
  Implicit None
  Real(realk), Dimension(3*NAtoms) :: Coordinates
  Real(realk), Dimension(NAtoms)   :: AtomMass
  Real(realk), Dimension(3)        :: CMVec
  Real(realk)                      :: TotMass
  Integer                          :: NAtoms, I
!
  CMVec = 0E0_realk
  TotMass = Sum(AtomMass(1:NAtoms))
  Do I = 1, NAtoms
    CMVec(1) = CMVec(1) + AtomMass(I)*Coordinates((I-1)*3+1)
    CMVec(2) = CMVec(2) + AtomMass(I)*Coordinates((I-1)*3+2)
    CMVec(3) = CMVec(3) + AtomMass(I)*Coordinates((I-1)*3+3)
  End Do
  CMVec = CMVec/TotMass
End Subroutine Center_of_Mass

!===================
!   Bond_Distance
!===================
Function Bond_Distance(NAtoms, Coordinates, Atom1, Atom2) Result(BondDist)
!
! Calculate the bond distance between two atoms
!
  Implicit None
  Real(realk)                      :: BondDist
  Real(realk), Dimension(3*NAtoms) :: Coordinates
  Real(realk), Dimension(3)        :: MoveVec
  Integer                               :: NAtoms, Atom1, Atom2
!
  BondDist = Sqrt( &
       (Coordinates((Atom1-1)*3+1)-Coordinates((Atom2-1)*3+1))**2 + &
       (Coordinates((Atom1-1)*3+2)-Coordinates((Atom2-1)*3+2))**2 + &
       (Coordinates((Atom1-1)*3+3)-Coordinates((Atom2-1)*3+3))**2)
End Function Bond_Distance

!=================
!   Fragment_CM
!=================
Subroutine Fragment_CM(NAtoms, Coordinates, AtomMass, &
                       FragNum, Fragments, CMVec)
!
! Determine center of mass for a fragment
!
  Implicit None
  Real(realk), Dimension(3*NAtoms) :: Coordinates
  Real(realk), Dimension(NAtoms)   :: AtomMass
  Real(realk), Dimension(3)        :: CMVec
  Real(realk)                      :: TotMass
  Integer                               :: NAtoms, FragNum, I
  Integer, Dimension(NAtoms)            :: Fragments
!
  CMVec = 0E0_realk; TotMass = 0E0_realk
  Do I = 1, NAtoms
    If (Fragments(I) == FragNum) Then
      CMVec(1) = CMVec(1) + AtomMass(I)*Coordinates((I-1)*3+1)
      CMVec(2) = CMVec(2) + AtomMass(I)*Coordinates((I-1)*3+2)
      CMVec(3) = CMVec(3) + AtomMass(I)*Coordinates((I-1)*3+3)
      TotMass = TotMass + AtomMass(I)
    End If
  End Do
  If (TotMass < 1E-12_realk) Then
    CMVec = 0E0_realk
  Else
    CMVec = CMVec/TotMass
  End If
End Subroutine Fragment_CM

!===================
!   Fragment_Size
!===================
Function Fragment_Size(NAtoms, Coordinates, FragNum, Fragments) &
         Result(Size_of_Frag)
!
! Determine size of a fragment, i.e. the distance between
! the two atoms in the fragment furthest apart
!
  Implicit None
  Real(realk), Dimension(3*NAtoms) :: Coordinates
  Real(realk)                      :: Size_of_Frag, Dist
  Integer                               :: NAtoms, FragNum, Atm1, Atm2
  Integer, Dimension(NAtoms)            :: Fragments
!
  Size_of_Frag = 0E0_realk
  Do Atm1 = 1, NAtoms - 1
    If (Fragments(Atm1) == FragNum) Then
      Do Atm2 = Atm1 + 1, NAtoms
        If (Fragments(Atm2) == FragNum) Then
	  Dist = Bond_Distance(NAtoms, Coordinates, Atm1, Atm2)
          If (Dist > Size_of_Frag) Size_of_Frag = Dist
        End If
      End Do
    End If
  End Do
End Function Fragment_Size

!=======================
!   Make_RightHandSys
!=======================
Subroutine Make_RightHandSys(Vectors)
!
! Ensure that three vectors constitue a right hand coordinate system
!
  Implicit None
  Real(realk), Dimension(3,3) :: Vectors
!
  Vectors(1,3) = Vectors(2,1)*Vectors(3,2) - Vectors(3,1)*Vectors(2,2)
  Vectors(2,3) = Vectors(3,1)*Vectors(1,2) - Vectors(1,1)*Vectors(3,2)
  Vectors(3,3) = Vectors(1,1)*Vectors(2,2) - Vectors(2,1)*Vectors(1,2)
End Subroutine Make_RightHandSys
!===================
!   Rotate_Matrix
!===================
Subroutine Rotate_Matrix(NAtoms, MatrixIn, MatrixOut, RMat, RotMat, &
                         TempMat, Mode)
!
! Rotates a 3N x 3N matrix based on the 3 x 3 rotational matrix RMat
!
  Implicit None
  Real(realk), Dimension(3*NAtoms,3*NAtoms) :: MatrixIn, MatrixOut, &
                                                    RotMat, TempMat
  Real(realk), Dimension(3,3)      :: RMat
  Integer                               :: NAtoms, Mode, I
!
! The 3x3 RMat is expanded to the full 3Nx3N rotational matrix, RotMat
!
  RotMat = 0E0_realk
  Do I = 1, NAtoms
    RotMat((I-1)*3+1:(I-1)*3+3, (I-1)*3+1:(I-1)*3+3) = RMat
  End Do
!
! A matrix is then rotated forward (Mode = 1) or backward (Mode = -1)
!
  If (Mode > 0) Then
    TempMat   = MatMul(RotMat,MatrixIn)
    MatrixOut = MatMul(TempMat,Transpose(RotMat))
  Else
    TempMat   = MatMul(Transpose(RotMat),MatrixIn)
    MatrixOut = MatMul(TempMat,RotMat)
  End If
End Subroutine Rotate_Matrix

!===================
!   Rotate_Vector
!===================
Subroutine Rotate_Vector(NAtoms, VectorIn, VectorOut, RMat, RotMat, &
                         TempVec, Mode)
!
! Rotates a 3N vector based on the 3 x 3 rotational matrix RMat
!
  Implicit None
  Real(realk), Dimension(3*NAtoms,3*NAtoms) :: RotMat
  Real(realk), Dimension(3,3)      :: RMat
  Real(realk), Dimension(3*NAtoms) :: VectorIn, VectorOut, TempVec
  Integer                               :: NAtoms, Mode, I
!
! The 3x3 RMat is expanded to the full 3Nx3N rotational matrix, RotMat
!
  RotMat = 0E0_realk
  Do I = 1, NAtoms
    RotMat((I-1)*3+1:(I-1)*3+3, (I-1)*3+1:(I-1)*3+3) = RMat
  End Do
!
! A vector is then rotated using RotMat (Mode = 1) or the
! transpose (inverse) of RotMat (Mode = -1)
!
  If (Mode > 0) Then
    TempVec = MatMul(RotMat,VectorIn)
  Else
    TempVec = MatMul(Transpose(RotMat),VectorIn)
  End If
  VectorOut = TempVec
End Subroutine Rotate_Vector

!====================
!   Print_Vector
!====================
Subroutine Print_Vector(FUnit, NAtoms, Labels, Vector)
!
! Prints out vector entities of the system including atomic labels (names)
!
  Implicit None
  Integer                               :: FUnit, NAtoms, I, J
  Character(Len = 4), Dimension(NAtoms) :: Labels
  Real(realk), Dimension(3*NAtoms) :: Vector
!
  Do I = 1, NAtoms
    Write(FUnit,'(1X,A,F17.10,2F24.10)') Labels(I), &
                                        (Vector(3*(I-1)+J), J = 1, 3)
  End Do
  Write(FUnit,'()')
End subroutine Print_Vector
!=================
!   Doublelines
!=================
Subroutine Doublelines(FUnit, Text, Indent)
!
! Puts double lines over and under the given text. If Indent < 0 it is also
! centered.
!
  Implicit None
  Integer            :: FUnit, Indent, TextLength, Indentation, I
  Character(Len = *) :: Text
!
  TextLength = Len(Text)
  Indentation = (72 - TextLength+6)/2 + 1
  If (Indent >= 0) Indentation = Indent + 1
  Write(FUnit, '(/,80A)') (' ', I = 1, Indentation), ('=', I = 1, TextLength+6)
  Write(FUnit, '(80A)') (' ', I = 1, Indentation+3), Text
  Write(FUnit, '(80A)') (' ', I = 1, Indentation), ('=', I = 1, TextLength+6)
  Write(FUnit, '()')
End Subroutine Doublelines

!=================
!   Singlelines
!=================
Subroutine Singlelines(FUnit, Text, Indent)
!
! Puts single lines over and under the given text. If Indent < 0 it is also
! centered.
!
  Implicit None
  Integer            :: FUnit, Indent, TextLength, Indentation, I
  Character(Len = *) :: Text
!
  TextLength = Len(Text)
  Indentation = (72 - TextLength+6)/2 + 1
  If (Indent >= 0) Indentation = Indent + 1
  Write(FUnit, '(/,80A)') (' ', I = 1, Indentation), ('-', I = 1, TextLength+6)
  Write(FUnit, '(80A)') (' ', I = 1, Indentation+3), Text
  Write(FUnit, '(80A)') (' ', I = 1, Indentation), ('-', I = 1, TextLength+6)
  Write(FUnit, '()')
End Subroutine Singlelines

!====================
!   DoublelinesInt
!====================
Subroutine DoublelinesInt(FUnit, Text, Number, LineLength)
!
! Puts singles lines over and under the given text + integer
!
  Implicit None
  Integer            :: FUnit, Number, LineLength, I
  Character(Len = *) :: Text
!
  Write(FUnit, '(/,A,80A)') ' ', ('=', I = 1, LineLength)
  Write(FUnit, '(A,A,I7)') '   ', Text, Number
  Write(FUnit, '(A,80A)') ' ', ('=', I = 1, LineLength)
  Write(FUnit, '()')
End Subroutine DoublelinesInt

!====================
!   SinglelinesInt
!====================
Subroutine SinglelinesInt(FUnit, Text, Number, LineLength)
!
! Puts singles lines over and under the given text + integer
!
  Implicit None
  Integer            :: FUnit, Number, LineLength, I
  Character(Len = *) :: Text
!
  Write(FUnit, '(/,A,80A)') ' ', ('-', I = 1, LineLength)
  Write(FUnit, '(A,A,I7)') '   ', Text, Number
  Write(FUnit, '(A,80A)') ' ', ('-', I = 1, LineLength)
  Write(FUnit, '()')
End Subroutine SinglelinesInt

!====================
!   Singlelines2Int
!====================
Subroutine Singlelines2Int(FUnit, Text, Number, Text2, Number2, LineLength)
!
! Puts singles lines over and under the given text + integer
!
  Implicit None
  Integer            :: FUnit, Number, Number2, LineLength, I
  Character(Len = *) :: Text, Text2
!
  Write(FUnit, '(/,A,80A)') ' ', ('=', I = 1, LineLength)
  Write(FUnit, '(A,A,I7,A,I7)') '   ', Text, Number, Text2, Number2
  Write(FUnit, '(A,80A)') ' ', ('=', I = 1, LineLength)
  Write(FUnit, '()')
End Subroutine Singlelines2Int

!===============
!   Underline
!===============
Subroutine Underline(FUnit, Text, Indent)
!
! Underlines the given text, if Indent < 0 it is also centered.
!
  Implicit None
  Integer            :: FUnit, Indent, TextLength, Indentation, I
  Character(Len = *) :: Text
!
  TextLength = Len(Text)
  Indentation = (72 - TextLength)/2 + 1
  If (Indent >= 0) Indentation = Indent + 1
  Write(FUnit, '(/,80A)') (' ', I = 1, Indentation), Text
  Write(FUnit, '(80A)') (' ', I = 1, Indentation), ('-', I = 1, TextLength)
  Write(FUnit, '()')
End Subroutine Underline

!=====================
!   Generate_Random
!=====================
!Subroutine Generate_Random
!
! Prints out the requested amount of random numbers
!
!  Use File_Units
!  Use Dynamics_Param
!  Use Dummies
!  Implicit None
!  Real(realk) :: Random_Num
!  Integer          :: I
!
!  Call Doublelines(FUPrint, 'Random number generation', 0)
!  Write(FUPrint,'(A,I15,/)') ' Chosen random generator seed  : ', RanSeed
!  Call Init_Random_Gen(0, RanSeed, Dummy, Dummy, Dummy, 0, MaxRanCnt)
!
!  Do I = 1, GetRandom
!    Write(FUPrint,'(5X,I15,F18.15)') I, Random_Num(0)*(UpperRandom-LowerRandom) + LowerRandom
!  End Do
!End Subroutine Generate_Random

!================
!   Random_Gen
!================
!Function Random_Gen(DummyArg) Result(RanGen)
!
!      Random_Gen      - Generates random numbers using a linear
!                        congruental uniform random generator.
!
!      Reference:  D.E. Knuth, "The art of computer programming",
!      Vol. 2 of "seminumerical algorithms" (Addison-Wesley, 1981).
!      2nd edition.(pp.32)
!
!  Use Random_Data
!  Implicit None
!  Real(realk) :: BaseInv, RanGen
!  Integer          :: DummyArg, ITemp
!  Integer          :: I, J
!  GenKey(1:8)   = (/ 45, 127, 149, 76, 45, 244, 81, 88 /)
!  GenKey2(1:16) = (/ (0, i = 1, 16) /)
!  BaseInv = 1.0E0_realk/Base
!  Do I = 2, 9
!    ITemp = 0
!    Do J = 1, I-1
!      ITemp = ITemp + GenKey(J)*RanKey(I-J)
!    End Do
!    ITemp = ITemp + GenKey2(I-1)
!    GenKey2(I-1) = Mod(ITemp,Base)
!    GenKey2(I)   = ITemp/Base
!  End Do
!  RanKey(1:8) = GenKey2(1:8)
!  RanGen = RanKey(1)
!  Do I  = 2, 8
!    RanGen = RanKey(I) + RanGen*BaseInv
!  End Do
!  RanGen = RanGen*BaseInv
!End Function Random_Gen

!================
!   Random_Num
!================
!Function Random_Num(Mode) Result(RanNum)
!
!      Random_Num     returns random number
!         Mode = 0    Return a random number.
!         Mode = 1    Spool ahead to next batch (determined by Count
!                     and MaxCount).
!
!      Reference:  D.E. Knuth, "The art of computer programming",
!      Vol. 2 of "seminumerical algorithms" (Addison-Wesley, 1981).
!      2nd edition.(pp.32)
!
!  Use Random_Data
!  Implicit None
!  Real(realk) :: RanNum, Random_Gen
!  Integer          :: ITemp, Mode
!  Integer          :: I
!  Select Case (Mode)
!    Case (0)
!      Count = Count + 1
!      If (Count >= MaxCount) Then
!        Count = 1
!        Write(*,*) &
!          'WARNING: Exceeded limit for batch of random numbers in Random_Num!'
!      End If
!      ITemp = Int(99E0_realk*RandomData(100))+1
!      RanNum = RandomData(ITemp)
!      RandomData(ITemp) = Random_Gen(0)
!    Case (1)
!      Do I = Count + 1, MaxCount
!        ITemp = Int(99E0_realk*RandomData(100))+1
!        RanNum = RandomData(ITemp)
!        RandomData(ITemp) = Random_Gen(0)
!      End Do
!      Count = 0
!    Case Default
!      Write(*,*) 'ERROR: Illegal mode selected for Random_Num!'
!      Stop
!  End Select
!End Function Random_Num

!====================
!   Element_Symbol
!====================
Subroutine Element_Symbol(AtomNumber, Symbol)
!
! Returns the one- or two-character elemental symbol for a given atomic number
!
  Implicit None
  Integer, Intent(In)                :: AtomNumber
  Character(Len = 2), Intent(Out)    :: Symbol
  Character(Len = 2), Dimension(103) :: AtomSymbol = &
    (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
       'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
       'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
       'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
       'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
       'Md', 'No', 'Lr' /)
  If ((AtomNumber < 1) .or. (AtomNumber > 103)) Then
    Symbol = 'XX'
  Else
    Symbol = AtomSymbol(AtomNumber)
  End If
End Subroutine Element_Symbol

!===================
!   Name_Fragment
!===================
Subroutine Name_Fragment(NAtoms, FragNum, Charges, Fragments, FragName)
!
! Constructs the chemical formula for either a fragment or the whole system
! (FragNum = 0). The formula is returned as a string (FragName).
!
  Implicit None
  Integer, Intent(In)                             :: NAtoms, FragNum
  Integer, Dimension(NAtoms), Intent(In)          :: Fragments
  Real(realk), Dimension(NAtoms), Intent(In) :: Charges
  Character(Len = 80), Intent(Out)                :: FragName
  Integer                                         :: Pos, I, J, Cnt, Digits
!
  Pos = 1
  FragName = '                                        ' &
          // '                                        '
  Do I = 1, 103
    Cnt = 0
    Do J = 1, NAtoms
      If ((NInt(Charges(J)) == I) .and. ((FragNum == 0) &
           .or. (Fragments(J) == FragNum))) Cnt = Cnt + 1
    End Do
    If (Cnt > 0) Then
      Call Element_Symbol(I, FragName(Pos:Pos+1))
      If (FragName(Pos+1:Pos+1) == ' ') Then
        Pos = Pos + 1
      Else
        Pos = Pos + 2
      End If
      If (Cnt > 1) Then
        Digits = Floor(Log10(1E0_realk*Cnt)) + 1
        Select Case(Digits)
          Case(1)
            Write(FragName(Pos:Pos),'(I1)') Cnt
          Case(2)
            Write(FragName(Pos:Pos+1),'(I2)') Cnt
          Case(3)
            Write(FragName(Pos:Pos+2),'(I3)') Cnt
          Case Default
            FragName(Pos:Pos+2) = '***'
            Digits = 3
        End Select
        Pos = Pos + Digits
      End If
    End If
  End Do
End Subroutine Name_Fragment

!====================
!   Analyze_Forces
!====================
!Subroutine Analyze_Forces(NAtm, Geom, Frc, TotFrc, AngFrc, Prnt, PrntTxt, FU)
!
! Analyze forces. Calculates both the total net force and the angular force
!
!  Use VecMat_Functions
!  Implicit None
!  Real(realk), Dimension(3*NAtm) :: Geom, Frc
!  Real(realk), Dimension(3)      :: TotFrc, AngFrc
!  Integer            :: NAtm, FU, I
!  Logical            :: Prnt
!  Character(Len = *) :: PrntTxt
!  TotFrc = 0E0_realk
!  AngFrc = 0E0_realk
!  Do I = 1, NAtm
!    TotFrc = TotFrc + Frc((I-1)*3+1:(I-1)*3+3)
!    AngFrc = AngFrc + Vector_Product(Geom((I-1)*3+1:(I-1)*3+3), &
!                                     Frc((I-1)*3+1:(I-1)*3+3))
!  End Do
!  If (Prnt) Write(FU, '(A,/,A,3E15.6,/,A,3E15.6)') &
!            PrntTxt, 'Net force    : ', TotFrc, 'Angular force: ', AngFrc
!End Subroutine Analyze_Forces

!==================
!   Swap_Vectors
!==================
Subroutine Swap_Vectors(NSwap, Vec1, Vec2)
!
! Swaps the NSwap first elements of vectors Vec1 and Vec2
!
  Implicit None
  Real(realk), Dimension(NSwap) :: Vec1, Vec2
  Real(realk) :: Temp
  Integer          :: NSwap, I
  Do I = 1, NSwap
    Temp = Vec1(I)
    Vec1(I) = Vec2(I)
    Vec2(I) = Temp
  End Do
End Subroutine Swap_Vectors

!====================
!   Schmidt_Orthog
!====================
!Subroutine Schmidt_Orthog(Dim, OrtBasis, TmpVc)
!
! Make othogonal basis through Schmidt orthogonalization
!
!  Use VecMat_Functions
!  Implicit None
!  Integer :: Dim
!  Real(realk), Dimension(Dim,Dim) :: OrtBasis
!  Real(realk), Dimension(Dim)     :: TmpVc
!  Real(realk) :: VecNrm, Comp
!  Integer          :: NFound, NVec, I
!!
!  If (Dim < 1) &
!    Call Quit('Number of dimensions in Schmidt_Orthog must be positive!')
!
!  TmpVc = OrtBasis(:,1)
!  VecNrm = Sqrt(Dot_Product(TmpVc,TmpVc))
!
! If the initial "seed" vector is a zero vector, OrtBasis is simply set equal
! a unit matrix
!  
!  If (VecNrm < 1E-12_realk) Then
!    OrtBasis = Unit_Matrix(Dim)
!    Return
!  Else
!    OrtBasis(:,1) = TmpVc/VecNrm
!  End If
!  NFound = 1
!  NVec = 1
!
! Run through the rest of the dimensions, building up an
! orthogonal basis
!
!  Do While ((NFound < Dim) .and. (NVec <= Dim))
!    TmpVc = 0E0_realk; TmpVc(NVec) = 1E0_realk
!    Do I = 1, NFound
!      Comp = Dot_Product(TmpVc,OrtBasis(:,I))
!      TmpVc = TmpVc - Comp*OrtBasis(:,I)
!    End Do
!    VecNrm = Sqrt(Dot_Product(TmpVc,TmpVc))
!    If (VecNrm > 1E-12_realk) Then
!      NFound = NFound + 1
!      OrtBasis(:,NFound) = TmpVc/VecNrm
!    End If
!    NVec = NVec + 1
!  End Do
!  If (NFound < Dim) Call Quit('Construction of orthogonal basis failed')
!End Subroutine Schmidt_Orthog
!
!====================!
! Calc_Kinetic_Cart  !
!====================!        
! Calculates kinetic energy for Cartesian velocities
Subroutine Calc_Kinetic_Cart(NumCoord,N,Mass,Vels,Kinetic)
!
Integer :: NumCoord, N
Real(realk) Mass(N)
Real(realk) Vels(NumCoord)
Real(realk) :: Kinetic
Real(realk), external :: DDOT
Integer :: i
  Kinetic = 0E-7_realk 
  print*,'Vels',Vels
  print*,'Mass',Mass
  Do i = 1,N
     Kinetic = Kinetic + (DDOT(3,Vels(3*i-2:3*i),1,Vels(3*i-2:3*i),1 )*Mass(i) )/2E0_realk
  Enddo
  print*,'Kinetic',Kinetic
End Subroutine Calc_Kinetic_Cart
!====================!
! Calc_Kinetic       !
!====================!        
! Calculates kinetic energy for mass-weighted Cartesian coordinates
Subroutine Calc_Kinetic(NumCoord,N,Vels,Kinetic)
!
Integer :: NumCoord, N
Real(realk) Vels(NumCoord)
Real(realk) :: Kinetic
Real(realk), external :: DDOT
Integer :: i
  Kinetic = 0E-7_realk  
  print*,'Vels',Vels
  Do i = 1,N
     Kinetic = Kinetic + DDOT(3,Vels(3*i-2:3*i),1,Vels(3*i-2:3*i),1 )/2E0_realk
  Enddo
  print*,'Kinetic',Kinetic
End Subroutine Calc_Kinetic
!===================!
! Calc_AngMom_Cart  !
!===================!
! Calculates angular momentum for Cartesian integration
Subroutine Calc_AngMom_Cart(NumCoord,N,Mass,Coord,Vels,AngMom)
Implicit none
Integer NumCoord,N
Real(realk) Mass(N)
Real(realk) Coord(NumCoord)
Real(realk) Vels(NumCoord)
Real(realk) AngMom(3)
Integer :: i
  AngMom = 0E-7_realk
  Do i = 1,N
     AngMom(1) = AngMom(1) + ( Coord(3*i-1)*Vels(3*i) - Coord(3*i)*Vels(3*i-1) )*Mass(i)
     AngMom(2) = AngMom(2) + ( Coord(3*i-2)*Vels(3*i) - Coord(3*i)*Vels(3*i-2) )*Mass(i)
     AngMom(3) = AngMom(3) + ( Coord(3*i-2)*Vels(3*i-1) - Coord(3*i-1)*Vels(3*i-2))*Mass(i)
  Enddo
End Subroutine Calc_AngMom_Cart 
!===================!
! Calc_AngMom       !
!===================!
! Calculates angular momentum for mass-weighted Cartesian integration
Subroutine Calc_AngMom(NumCoord,N,Coord,Vels,AngMom)
Implicit none
Integer NumCoord,N
Real(realk) Coord(NumCoord)
Real(realk) Vels(NumCoord)
Real(realk) AngMom(3)
Integer :: i
  AngMom = 0E-7_realk
  Do i = 1,N
     AngMom(1) = AngMom(1) + ( Coord(3*i-1)*Vels(3*i) - Coord(3*i)*Vels(3*i-1) )
     AngMom(2) = AngMom(2) + ( Coord(3*i-2)*Vels(3*i) - Coord(3*i)*Vels(3*i-2) )
     AngMom(3) = AngMom(3) + ( Coord(3*i-2)*Vels(3*i-1) - Coord(3*i-1)*Vels(3*i-2))
  Enddo
End Subroutine Calc_AngMom 
!=======================
!   Write_PhaseSpace
!=======================
Subroutine Write_PhaseSpace(FUnit, NAtoms, Charges, Coordinates, &
                          Velocities, Time, Energy, PotEnergy, KinEnergy, &
                          EnergyCons, AngMomCons)
!
! Writes energy and phase space information
!
  Implicit None
  Integer                               :: FUnit, NAtoms, I, J
  Real(realk), Dimension(NAtoms)   :: Charges
  Real(realk), Dimension(3*NAtoms) :: Coordinates, Velocities
  Real(realk) :: Time, Energy, PotEnergy, KinEnergy, EnergyCons, &
                      AngMomCons
  Write(FUnit,'(F12.5,A3,2F15.8,F12.8,A3,F13.10,F17.14)') Time, ' * ', &
                Energy, PotEnergy, KinEnergy, ' * ', EnergyCons, AngMomCons
  Do I = 1, NAtoms
    Write(FUnit,'(I6,F12.6,A3,3F23.15)') I, Charges(I), ' C ', &
                                        (Coordinates(3*(I-1)+J), J = 1, 3)
    Write(FUnit,'(I6,F12.6,A3,3F23.15)') I, Charges(I), ' V ', &
                                        (Velocities(3*(I-1)+J), J = 1, 3)
  End Do
End Subroutine Write_PhaseSpace
!
!====================
!   Final_Analysis
!====================
Subroutine Final_Analysis(NumAtoms,MaxSteps,NumTrajs,FUPrint,Phase,PrintLevel)
!
! back in from the file DALTON.PHS
!
  Implicit None
  Integer :: FUPrint,Phase  ! File units
  Integer :: NumAtoms
  Integer :: PrintLevel
  Real(realk), pointer :: Coords(:,:,:,:)
  Real(realk), pointer :: Time(:,:), TotEnrg(:,:), PotEnrg(:,:), &
                       KinEnrg(:,:), EnrgCons(:,:), AngMomCons(:,:)
  Real(realk), pointer ::  AtmArr(:,:)
  Real(realk)    :: Sumt,SumE,SumEt,Sumt2,Noise,Drift,Free_term   ! Needed for noise and drift
  Integer, Dimension(NumTrajs)                   :: NumSteps
  Integer             :: MaxSteps, NumTrajs, FirstTraj, &
                         NumTraj, Iteration, I, J, K, Stride, BondLim
  Integer             :: MCoord, MBonds, MBndAt, MEnd, Atom1, Atom2, NumBonds
  Logical             :: First, DoIt
  Character(Len = 2)  :: Elem1, Elem2
  Character(Len = 30) :: HeaderText
  Character(Len = 64) :: PrnString, PrnString2, PrnString3
  Character(Len = 80) :: FragTxt
!
  Call Doublelines(FUPrint, 'Final analysis of ab initio trajectories', 0)
!  Call AtmIni(AtmArr, I, .True.)
! Allocate memory
  Call mem_alloc(Coords,NumTrajs,MaxSteps+1,NumAtoms,3)
  Call mem_alloc(Time,NumTrajs,MaxSteps+1)
  Call mem_alloc(TotEnrg,NumTrajs,MaxSteps+1)
  Call mem_alloc(PotEnrg,NumTrajs,MaxSteps+1)
  Call mem_alloc(KinEnrg,NumTrajs,MaxSteps+1)
  Call mem_alloc(EnrgCons,NumTrajs,MaxSteps+1)
  Call mem_alloc(AngMomCons,NumTrajs,MaxSteps+1)
  Call mem_alloc(AtmArr,NumAtoms,8)
!
! Read information
!
  Call Read_PhaseFile(NumAtoms,Phase,FUPrint,PrintLevel,NumSteps,Time,&
                     &TotEnrg,PotEnrg,KinEnrg,EnrgCons,AngMomCons,&
                     &FirstTraj,NumTrajs,Coords)
!
! Print summary of each trajectory
!
  Do I = 1, NumTrajs
    HeaderText = 'Summary of trajectory #       '
    Write(HeaderText(24:30),'(I7)') I+FirstTraj-1
    Call Singlelines(FUPrint,HeaderText,-1)
    Stride = 1 + (NumSteps(I)-1)/30
!
! Calculating energy drift and noise
!   
    ! Setting entities equal to 0 (since differ for each traj.)
    Sumt = 0E0_realk
    Sumt2 = 0E0_realk
    SumEt = 0E0_realk
    SumE = 0E0_realk
    Drift = 0E0_realk
    Noise = 0E0_realk
    Free_term = 0E0_realk
    Do J = 1, NumSteps(I)
       Sumt = Sumt + Time(I,J)
       Sumt2 = Sumt2 + Time(I,J)**2
       SumE = SumE + TotEnrg(I,J) 
       SumEt = SumEt + TotEnrg(I,J)*Time(I,J)
    Enddo
    Drift = (NumSteps(I)*SumEt - SumE*Sumt)/(NumSteps(I)*Sumt2 - Sumt**2)
    Free_term = (SumE*Sumt2 - Sumt*SumEt)/(NumSteps(I)*Sumt2 - Sumt**2)
    ! Calculating noise
    Do J = 1, NumSteps(I)
       Noise = Noise + ( TotEnrg(I,J) - (Time(I,J)*Drift + Free_term) )**2
    Enddo
    Noise = (Noise/NumSteps(I))**0.5E0_realk 
    ! Scaling to mHartree
    Write(*,*) Free_term
    Write(*,*) Drift
    Drift = (Drift*10E0_realk**9E0_realk)
    Write(*,*) Drift
    Noise = Noise*10E0_realk**6E0_realk
    ! 
! If there's more than 30 steps, we print a long and a short version of
! the energy and conservation information
!
    Write(FUPrint,*)
    Write(FUPrint,'(5X,A)') '+----------------------------------------' & 
            // '------------------------------------------------+'
    Write(FUPrint,'(5X,A)') '|                 Energy and conservation' &
            // ' information - Complete version                 |'
    Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
            // '------------------+-----------------------------+'
    Write(FUPrint,'(5X,A)') '|                  |     Total       Potential' &
            // '    Kinetic  |   Energy     Angular mom.   |'
    Write(FUPrint,'(5X,A)') '|   Step    Time   |     energy       energy' &
            // '      energy   |  conserv.      conserv.     |'
    Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
            // '------------------+-----------------------------+'
    If (Stride > 1) Then
      Do J = 1, NumSteps(I)
        Write(FUPrint,'(5X,A,I7,F10.3,A,2F13.6,F12.8,A,F11.8,F17.14,A)') &
                      '|', J-1, Time(I,J), ' |', TotEnrg(I,J), PotEnrg(I,J), &
                      KinEnrg(I,J), ' |', EnrgCons(I,J), AngMomCons(I,J), ' |'
      End Do
      Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
            // '------------------+-----------------------------+'
      Write(FUPrint,*)
      Write(FUPrint,'(5X,A)') '+----------------------------------------' & 
            // '------------------------------------------------+'
      Write(FUPrint,'(5X,A)') '|                  Energy and conservation' &
            // ' information - Short version                   |'
      Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
            // '------------------+-----------------------------+'
      Write(FUPrint,'(5X,A)') '|                  |     Total       ' &
            // 'Potential    Kinetic  |   Energy     Angular mom.   |'
      Write(FUPrint,'(5X,A)') '|   Step    Time   |     energy       energy' &
            // '      energy   |  conserv.      conserv.     |'
      Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
            // '------------------+-----------------------------+'
    End If
    Do J = 1, NumSteps(I), Stride
      Write(FUPrint,'(5X,A,I7,F10.3,A,2F13.6,F12.8,A,F11.8,F17.14,A)') &
                    '|', J-1, Time(I,J), ' |', TotEnrg(I,J), PotEnrg(I,J), &
                    KinEnrg(I,J), ' |', EnrgCons(I,J), AngMomCons(I,J), ' |'
    End Do
    If ((Stride > 1) .and. (Mod(NumSteps(I)-1, Stride) /= 0)) Then
      Write(FUPrint,'(5X,A,I7,F10.3,A,2F13.6,F12.8,A,F11.8,F17.14,A)') &
                    '|', NumSteps(I)-1, Time(I,NumSteps(I)), ' |', &
                    TotEnrg(I,NumSteps(I)), PotEnrg(I,NumSteps(I)), &
                    KinEnrg(I,NumSteps(I)), ' |', EnrgCons(I,NumSteps(I)), &
                    AngMomCons(I,NumSteps(I)), ' |'
    End If
    Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
                      // '------------------+-----------------------------+'
!
! Max and min values are also printed
!
    Write(FUPrint,'(5X,A,2F13.6,F12.8,A,F11.8,F17.14,A)') &	
      '|          Minimum |', MinVal(TotEnrg(I,1:NumSteps(I))), &
      MinVal(PotEnrg(I,1:NumSteps(I))), MinVal(KinEnrg(I,1:NumSteps(I))), &
      ' |', MinVal(EnrgCons(I,2:NumSteps(I))), &
      MinVal(AngMomCons(I,2:NumSteps(I))),' |'
    Write(FUPrint,'(5X,A,2F13.6,F12.8,A,F11.8,F17.14,A)') &	
      '|          Maximum |', MaxVal(TotEnrg(I,1:NumSteps(I))), &
      MaxVal(PotEnrg(I,1:NumSteps(I))), MaxVal(KinEnrg(I,1:NumSteps(I))), &
      ' |', MaxVal(EnrgCons(I,2:NumSteps(I))), &
      MaxVal(AngMomCons(I,2:NumSteps(I))),' |'
    Write(FUPrint,'(5X,A)') '+------------------+---------------------' & 
                      // '------------------+-----------------------------+'
      ! Writing noise and drift
    Write(FUPrint,'(A,F10.3)') ' Energy drift, mHartree/ps  : ', Drift
    Write(FUPrint,'(A,F10.3)') ' Energy noise, mHartree     : ', Noise
      !
!
! Bond analysis
!
!      DoIt = .True.
!      BondLim = 8
!      Atom1 = 1
!      Atom2 = 1
!      Write(FUPrint,*)
!      Write(FUPrint,'(7X,A)') '+----------------------------------------' & 
!            // '--------------------------------------------+'
!      Write(FUPrint,'(7X,A)') '|                          Bond length ' &
!            // 'analysis (Aangstrom)                          |'
!      Write(FUPrint,'(7X,A)') '+------------------+---------------------' &
!            // '--------------------------------------------+'
!      MCoord = 1
!      MBonds = MCoord + NumTrajs*MaxSteps*NumAtoms*3
!      MBndAt = MBonds + (MaxSteps+2)*BondLim
!      MEnd   = MBndAt + 2*BondLim
!      If (MEnd > LWork) Then
!        Write(FUPrint,'(7X,A)') '|               NOTE: Not enough memory ' &
!              // 'to perform full bond analysis!               |'
!        Write(FUPrint,'(7X,A,I12,A,I12,A)') '|                     ' &
!              // 'Available:', LWork, '   Needed:', MEnd, &
!                 '                   |'
!        DoIt = .False.
!      End If
!      First = .True.
!      Do While (DoIt .and. (Atom1*Atom2 < (NumAtoms-1)*NumAtoms))
!        Call Calc_Bonds(NumBonds, BondLim, NumTrajs, MaxSteps, NumAtoms, I, &
!                        Atom1, Atom2, NumSteps, Time, AtmArr(1:NumAtoms,5), &
!                        WorkArray(MCoord), &
!                        WorkArray(MBonds), WorkArray(MBndAt))
!        If (NumBonds > 0) Then
!          PrnString = &
!            '                                                                '
!          PrnString2 = PrnString
!          PrnString3 = PrnString
!          Do K = 1, NumBonds
!            Call Element_Symbol(NInt(AtmArr(Nint(WorkArray(MBndAt+ &
!                                K-1)),1)), Elem1)
!            Call Element_Symbol(NInt(AtmArr(Nint(WorkArray(MBndAt+ &
!                                BondLim+K-1)),1)), Elem2)
!            If (Elem1(2:2) == ' ') Then
!              Write(PrnString((K-1)*8+5:K*8),'(3A)') Elem1(1:1), '-', Elem2
!            Else
!              Write(PrnString((K-1)*8+4:K*8),'(3A)') Elem1, '-', Elem2
!            End If
!            Write(PrnString2((K-1)*8+1:K*8),'(A5,I3)') &
!                  '  Atm', Nint(WorkArray(MBndAt+K-1))
!            Write(PrnString3((K-1)*8+1:K*8),'(A5,I3)') &
!                  '  Atm', Nint(WorkArray(MBndAt+BondLim+K-1))
!          End Do
!          If (.Not. First) Write(FUPrint,'(7X,A)') '+==================+' &
!        // '=================================================================+'
!          Write(FUPrint,'(7X,3A)') '|                  |', PrnString, ' |'
!          Write(FUPrint,'(7X,3A)') '|                  |', PrnString2, ' |'
!          Write(FUPrint,'(7X,3A)') '|   Step    Time   |', PrnString3, ' |'
!          Write(FUPrint,'(7X,A)') '+------------------+---------------------' &
!            // '--------------------------------------------+'
!          PrnString2 = &
!            '                                                                '
!          Do J = 1, NumSteps(I), Stride
!            PrnString = PrnString2
!            Write(PrnString,'(40F8.3)') &
!              (WorkArray(MBonds+(J-1)*BondLim+K-1), K=1,NumBonds)
!            Write(FUPrint,'(7X,A,I7,F10.3,3A)') &
!                    '|', J-1, Time(I,J), ' |', PrnString, ' |'
!          End Do
!          If ((Stride > 1) .and. (Mod(NumSteps(I)-1, Stride) /= 0)) Then
!            PrnString = PrnString2
!            Write(PrnString,'(40F8.3)') &
!              (WorkArray(MBonds+(NumSteps(I)-1)*BondLim+K-1), K=1,NumBonds)
!            Write(FUPrint,'(7X,A,I7,F10.3,3A)') &
!              '|', NumSteps(I)-1, Time(I,NumSteps(I)), ' |', PrnString, ' |'
!          End If
!          Write(FUPrint,'(7X,A)') '+------------------+---------------------' &
!            // '--------------------------------------------+'
!          PrnString = PrnString2
!          Write(PrnString,'(40F8.3)') &
!              (WorkArray(MBonds+MaxSteps*BondLim+K-1), K=1,NumBonds)
!          Write(FUPrint,'(7X,3A)') '|          Minimum |', &
!               PrnString, ' |'
!          PrnString = PrnString2
!          Write(PrnString,'(40F8.3)') &
!              (WorkArray(MBonds+(MaxSteps+1)*BondLim+K-1), K=1,NumBonds)
!          Write(FUPrint,'(7X,3A)') '|          Maximum |', &
!               PrnString, ' |'
!          First = .False.
!      End Do
!      Write(FUPrint,'(7X,A)') '+------------------+---------------------' &
!            // '--------------------------------------------+'
!
! Fragment analysis
!
!      DoIt = .True.
!      Write(FUPrint,*)
!      Write(FUPrint,'(7X,A)') '+----------------------------------------' & 
!            // '--------------------------------------------+'
!      Write(FUPrint,'(7X,A)') '|                                 Fragment ' &
!            // 'analysis                                  |'
!      Write(FUPrint,'(7X,A)') '+----------------------------------------' &
!            // '--------------------------------------------+'
!      MBonds = 1
!      MEnd   = MBonds + NumAtoms*NumAtoms
!      If (MEnd > LWork) Then
!        Write(FUPrint,'(7X,A)') '|               NOTE: Not enough memory ' &
!              // 'to perform fragment analysis!                |'
!        Write(FUPrint,'(7X,A,I12,A,I12,A)') '|                     ' &
!              // 'Available:', LWork, '   Needed:', MEnd, &
!                 '                   |'
!        DoIt = .False.
!      End If
!      If (DoIt) Then
!        Write(FUPrint,'(7X,A)') '|                                   Coming' &
!              // ' soon!                                     |'
!      End If
!      Write(FUPrint,'(7X,A)') '+----------------------------------------' &
!            // '--------------------------------------------+'
!    End If
  End Do
!
! Print summary of all trajectories
!
!  Write(FUPrint,*)
!  Call Print_FullSummary(MaxSteps, NumTrajs, NumSteps, Time, TotEnrg, &
!                      PotEnrg, KinEnrg, EnrgCons, AngMomCons, FirstTraj, NumTraj)
!  Write(FUPrint,*)
!  Write(FUPrint,*)
  Call LSClose(Phase,'KEEP')
! Deallocate memory
  Call mem_dealloc(Coords)
  Call mem_dealloc(Time)
  Call mem_dealloc(TotEnrg)
  Call mem_dealloc(PotEnrg)
  Call mem_dealloc(KinEnrg)
  Call mem_dealloc(EnrgCons)
  Call mem_dealloc(AngMomCons)
  Call mem_dealloc(AtmArr)
!
End Subroutine Final_Analysis
!===============
! Print temp
!===============
Subroutine Print_temp(num_steps,T_array,lupri)
! Prints temperature statistics
Implicit none
Integer :: num_steps, i, lupri
Real(realk), dimension(num_steps+1) :: T_array
!
Call Underline(lupri, 'Temperature variation', -1)
Write(lupri,'(5X,A)') '+--------+-----------+'  
Write(lupri,'(5X,A)') '|  Step  |    T, K   |'  
Write(lupri,'(5X,A)') '+--------+-----------+'  
Do i = 1, num_steps+1 
     Write(lupri,'(5X,A,I7,A,F10.3,A)') '|', i-1,  ' |', T_array(i),' |'
Enddo
Write(lupri,'(5X,A)') '+--------+-----------+'  
!
End subroutine Print_temp
!====================
!   Read_PhaseFile
!====================
Subroutine Read_PhaseFile(NumAtoms,Phase,FUPrint,PrintLevel,NumSteps,Time,&
                          TotEnrg, PotEnrg,KinEnrg,EnrgCons,AngMomCons,&
                          FirstTraj,NumTrajs,Coords)
!
! Read information from the file DALTON.PHS
!
  Implicit None
  Integer :: NumAtoms
  Integer :: FUPrint
  Integer :: Phase  ! File unit for DALTON.PHS
  Real(realk)  Coords(:,:,:,:)
  Real(realk)  Time(:,:), TotEnrg(:,:), PotEnrg(:,:), &
             KinEnrg(:,:), EnrgCons(:,:),AngMomCons(:,:)
  Real(realk)                                :: NewTime
  Integer, Dimension(NumTrajs)                    :: NumSteps
  Integer          :: PrintLevel, NumTrajs, FirstTraj,&
                      Iteration, I, J, K
!
  Call LSOpen(Phase,'DALTON.PHS','OLD ','FORMATTED')
  Read(Phase,'(TR11,I7,TR8,I7)') FirstTraj, NumTrajs
  Write(FUPrint,'(A,I7,A)') ' DALTON.PHS contains ', NumTrajs, ' trajectories'
  Write(FUPrint,'(A,I7)')   ' The first trajectory is trajectory number ', &
                              FirstTraj
  Write(FUPrint,'(A)')      ' Starting to read trajectory information'
!
! Loop over the total number of trajectories
!
  Do I = 1, NumTrajs
    Iteration = 0
    Read(Unit=Phase, Fmt="(F12.5)", Advance="no") NewTime
!
! A negative time indicates the end of the current trajectory
!
    Do While (NewTime > -1E0_realk)
      Iteration = Iteration + 1
      Time(I, Iteration) = NewTime
!
! Read energies and conservation information, then loop over atoms to read
! coordinates and velocities (only kept if there's enough memory available).
!
      Read(Unit=Phase, Fmt="(TR4,2F15.8,F12.8,TR3,F13.10,F17.14)", &
           Advance="yes") TotEnrg(I, Iteration), PotEnrg(I, Iteration), &
           KinEnrg(I, Iteration), EnrgCons(I, Iteration), &
           AngMomCons(I, Iteration)
!      If (BondAnalysis) Then
        Do J =  1, NumAtoms
          Read(Phase,'(TR21,3F23.15)') (Coords(I, Iteration, J, K), K=1,3)
          Read(Phase,*)
        End Do
!      Else
!        Do J =  1, 2*NumAtoms
!          Read(Phase,*)
!        End Do
!      End If
      Read(Unit=Phase, Fmt="(F12.5)", Advance="no") NewTime
    End Do
    If (PrintLevel >= 3) Write(FUPrint,'(A,I7)') &
      ' Finished reading trajectory information for trajectory ', &
      FirstTraj + I - 1
    NumSteps(I) = Iteration
    If (I < NumTrajs) Then
      Read(Phase,*)
      Read(Phase,*)
    End If
  End Do
End Subroutine Read_PhaseFile
!=======================!
! Project_Gradient      !
!=======================!
! Removes components of the mass-weighted gradient which break
! translational-rotational symmetry
Subroutine Project_Gradient(NAtoms,Coordinates,mass,MW_Gradient,print_level,lupri)
Implicit none
Integer :: lupri, NAtoms, print_level
Real(realk) :: Coordinates(3*NAtoms)
Real(realk) :: MW_Gradient(3*NAtoms)
Real(realk) :: mass(NAtoms)
Real(realk) :: Principle_Axes(3,3)
Real(realk) :: Principle_Inertia(3)
Real(realk) :: TraRotVec(NAtoms*3,6)
Real(realk) :: TransformPM(NAtoms*3,NAtoms*3) ! Transformation to principle axes
Real(realk) :: TransformPMT(NAtoms*3,NAtoms*3) ! Transpose of TransformPM
Real(realk) :: Proj_MatPM(NAtoms*3,NAtoms*3)  ! Projection in principle axes frame
Real(realk) :: Proj_Mat(NAtoms*3,NAtoms*3)  ! Projection in general frame
Real(realk) :: TotalMass, VecNorm
Real(realk) :: Scale_Vector(NAtoms*3)
Integer :: TraRotDim, i, j, k, RotCount
Logical :: linear  ! Is molecule linear?
! Calculate total mass
TotalMass = 0E0_realk
Do i = 1,NAtoms
   TotalMass = TotalMass + Mass(i)
Enddo
! Remove mass-weighting (Mass_Weighted coordinates are sent to the subroutine!)
Call Mass_weight_vector(NAtoms,Coordinates,Mass,'REMOVE') 
!
Call Rotational_analysis(NAtoms,Coordinates,mass,Principle_Axes,Principle_Inertia,Linear)
!
TraRotVec = 0E0_realk
TraRotDim = 6
If (linear) TraRotDim = 5
Call Make_RightHandSys(Principle_Axes)
! Set up transformation from regular principle axes
TransformPM = 0E0_realk
Do i=1, NAtoms
   TransformPM((I-1)*3+1:(I-1)*3+3, (I-1)*3+1:(I-1)*3+3) = Principle_Axes
Enddo
!Call Underline(lupri, 'Transformation matrix to principle axes', -1)
!Call Output(TransfrmPM, 1, NumCoord, 1, NumCoord, NumCoord, &
!            NumCoord, 1, lupri)
! Set up the translational vectors
TraRotVec(1:(NAtoms*3):3, 1) = Sqrt(Mass(1:NAtoms))/Sqrt(TotalMass)
TraRotVec(2:(NAtoms*3):3, 2) = Sqrt(Mass(1:NAtoms))/Sqrt(TotalMass)
TraRotVec(3:(NAtoms*3):3, 3) = Sqrt(Mass(1:NAtoms))/Sqrt(TotalMass)

! Mass-weight coordinates
Call Mass_weight_vector(NAtoms,Coordinates,Mass,'WEIGHT') 
! Set up the rotation vectors

!For some strange reason this does not work:
!Scale_Vector = MATMUL(Transpose(TransformPM),Coordinates)
!so we do this instead:
TransformPMT = Transpose(TransformPM)
Scale_Vector = MATMUL(TransformPMT,Coordinates)
RotCount = 0

Do I = 1, 3
  If (Abs(Principle_Inertia(I)) > 1E-12_realk) Then
    RotCount = RotCount + 1
    Do J = 1, NAtoms
      Select Case (I)
        Case (1)
          TraRotVec((J-1)*3+2, RotCount+3) = -Scale_Vector((J-1)*3+3)
          TraRotVec((J-1)*3+3, RotCount+3) = Scale_Vector((J-1)*3+2)
        Case (2)
          TraRotVec((J-1)*3+1, RotCount+3) = Scale_Vector((J-1)*3+3)
          TraRotVec((J-1)*3+3, RotCount+3) = -Scale_Vector((J-1)*3+1)
        Case (3)
          TraRotVec((J-1)*3+1, RotCount+3) = -Scale_Vector((J-1)*3+2)
          TraRotVec((J-1)*3+2, RotCount+3) = Scale_Vector((J-1)*3+1)
      End Select
    End Do
    VecNorm = Sqrt(Dot_Product(TraRotVec(1:NAtoms*3,RotCount+3), &
                              TraRotVec(1:NAtoms*3,RotCount+3)))
    TraRotVec(1:NAtoms*3,RotCount+3) = TraRotVec(1:NAtoms*3,RotCount+3)/VecNorm
  End If
End Do
! Projection matrix in the principle inertia axes
Proj_MatPM = -MatMul(TraRotVec, Transpose(TraRotVec))
Do I = 1, NAtoms*3
  Proj_MatPM(I,I) = Proj_MatPM(I,I) + 1E0_realk
End Do

If (print_level >= 9) Then
!  Call Underline(lupri, 'MW geometry in principal axes', -1)
!  Call Output(Scale_Vector, 1, 1, 1, NAtoms*3, 1, NAtoms*3, 1, lupri)
!  Call Underline(lupri, &
!                'Translation and rotation vectors in MW Cartesians', -1)
!  Call Output(TraRotVec, 1, NAtoms*3, 1, TraRotDim, NAtoms*3, 6, 1, lupri)
!  Call Underline(lupri, 'Projection matrix in principal axes', -1)
!  Call Output(Proj_MatPM, 1, NAtoms*3, 1, NAtoms*3, NAtoms*3, NAtoms*3, &
!              1, lupri)
End If
!
! This projection matrix is then contracted with the transformation matrix
! from regular Cartesians to the principal axes.
! 
Proj_Mat = MatMul(Proj_MatPM, Transpose(TransformPM))
Proj_MatPM = Proj_Mat
Proj_Mat = MatMul(TransformPM, Proj_MatPM)
If (print_level >= 3) Then
!    Call Underline(lupri, 'Projection matrix', -1)
!    Call Output(Proj_Mat, 1, NAtoms*3, 1, NAtoms*3, NAtoms*3, NAtoms*3, &
!                1, lupri)
End If
! Do the projection
MW_Gradient = MATMUL(Proj_Mat,MW_Gradient)
!
End subroutine Project_Gradient
!=======================! 
! Rotational_analysis   !
!=======================!
Subroutine Rotational_analysis(NAtoms,Coordinates,mass,Principle_Axes,Principle_Inertia,Linear)
Implicit none
Integer :: NAtoms, lupri
Real(realk) :: Coordinates(3,NAtoms)
Real(realk) :: Mass(NAtoms)
Real(realk) :: Principle_Inertia(3) ! Principle moments of inertia
Real(realk) :: Principle_Axes(3,3)
Real(realk), pointer :: Work(:)
Real(realk) :: CMO(3)  ! Centre of mass origin
Real(realk) :: Inertia_tensor(3,3)
Integer :: i,j,Info,LWork
Real(realk), parameter :: TestLin = 1E-04_realk
Logical :: Linear
Info=0
! Moving origin to centre of mass
Call Center_of_Mass(NAtoms,Coordinates,Mass,CMO)
Call Move_Molecule(NAtoms,Coordinates,-CMO) 
!
! Build the inertia tensor in CMO frame
!
Inertia_tensor = 0E0_realk
Do i = 1, NAtoms  
   ! Diagonal elements
   Inertia_tensor(1,1) = Inertia_tensor(1,1) + &
   &  mass(i)*(Coordinates(2,i)**2+Coordinates(3,i)**2)
   Inertia_tensor(2,2) = Inertia_tensor(2,2) + &
   &  mass(i)*(Coordinates(1,i)**2+Coordinates(3,i)**2)
   Inertia_tensor(3,3) = Inertia_tensor(3,3) + &
   &  mass(i)*(Coordinates(1,i)**2+Coordinates(2,i)**2)
   ! Off-diagonal elements (upper triangular)
   Inertia_tensor(1,2) = Inertia_tensor(1,2) - &
   &  mass(i)*Coordinates(1,i)*Coordinates(2,i)
   Inertia_tensor(1,3) = Inertia_tensor(1,3) - &
   &  mass(i)*Coordinates(1,i)*Coordinates(3,i)
   Inertia_tensor(2,3) = Inertia_tensor(2,3) - &
   &  mass(i)*Coordinates(2,i)*Coordinates(3,i)
Enddo       
! Copy to lower-triangular part
Inertia_tensor(3,1) = Inertia_tensor(1,3)
Inertia_tensor(3,2) = Inertia_tensor(2,3)
Inertia_tensor(2,1) = Inertia_tensor(1,2)
! Find principle axes of inertia
LWork = 15
Call mem_alloc(Work,LWork)
Call dsyev('V','U',3,Inertia_Tensor,3,Principle_Inertia,Work,LWork,Info)
Call mem_dealloc(Work)
! Sort the eigenvalues and eigenvectors
! NB! Inertia_Tensor after DSYEV contains eigenvectors
Call ls_ORDER2(Inertia_Tensor,Principle_Inertia,3,3)
! Find whether the molecule is linear of not
If (Principle_Inertia(1) .LT. TestLin) then
   Linear = .TRUE.
Else
   Linear = .FALSE.
Endif
!
Principle_Axes = Inertia_Tensor
! Moving origin back
Call Move_Molecule(NAtoms,Coordinates,CMO) 
!
contains
  SUBROUTINE LS_ORDER2(EVEC,EVAL,N,NEVEC)    
!
!     Copied from gphjj.F  _ f90 by TK
!
! Purpose: order the N values in optinfo%EVAL and their associated vectors
!          in EVEC so optinfo%EVAL(i+1) .ge. optinfo%EVAL(i)
!
! Revisions:
!   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
!    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))
!   27-Oct-1984 hjaaj (reduced number of swaps)
!
    implicit none
!    Implicit Real(realk) (A-H,O-Z)
    Real(realk) :: EVEC(*),EVAL(*)
    integer :: N,NEVEC
!
    integer :: IN,I,IMIN,J
    Real(realk) :: EMIN

    IF (N.LE. 1) RETURN
    IN = 1
    DO I=1,N-1
       EMIN = EVAL(I)
       IMIN = I
       DO J=I+1,N
          IF (EVAL(J) .LT. EMIN) THEN
             EMIN = EVAL(J)
             IMIN = J
          ENDIF
       ENDDO
       IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
             CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
       ENDIF
       IN = IN + NEVEC
    ENDDO
  END SUBROUTINE LS_ORDER2

End subroutine Rotational_analysis
!======================!
!  Mass_weight_vector  !
!======================!
! Mass-weights a vector,or removes
! mass-weighting
Subroutine Mass_weight_vector(NAtoms,Vector,Mass,Mode)
Implicit none
Integer,intent(in) :: NAtoms
Real(realk),intent(inout) :: Vector(NAtoms*3)
Real(realk),intent(in) :: Mass(NAtoms)
Character(len=6),intent(in) :: Mode
Integer ::  i
! Do mass-weighting
If (Mode .EQ. 'WEIGHT' ) then
   Do i = 1, NAtoms
      Vector(i*3-2:i*3)=sqrt(Mass(i))*(Vector(i*3-2:i*3))
   Enddo
Endif   
! Remove mass-weighting
If (Mode .EQ. 'REMOVE' ) then
   Do i = 1, NAtoms
      Vector(i*3-2:i*3)=Vector(i*3-2:i*3)/sqrt(Mass(i))
   Enddo
Endif   
End subroutine
!=====================!
!  Init_random_seed   !
!=====================!
Subroutine init_random_seed()
  implicit none
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            i=0
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
End subroutine
!=======================!
! Cubic_equation_root   !
!=======================!
! Picks one real root of the cubic equation:
! x**3+ax**2+bx+c=0
! Trigonometric Viet formula is used.
Function Cubic_equation_root(a,b,c) 
Implicit none
Real(realk) :: Cubic_equation_root
Real(realk) :: a,b,c  ! Equation coefficients
Real(realk) :: Q,R,S,phi ! Needed for solution
!
Q = (a**2-3*b)/9
R = (2*a**2-9*a*b+27*c)/54
S = Q**3-R**2
!
If (S .GT. 0.00D0) then
   phi = (1/3)*acos(R/Q**(3/2))
   Cubic_equation_root = -2*sqrt(Q)*cos(phi) - a/3
Endif
!
!If (S .LT. 0.00D0) then
!   If (Q .GT. 0.00D0) then
!   Else
!Endif
!
!If (S .EQ. 0.00D0) then
!
!Endif
!
End function Cubic_equation_root
!
!=====================!
! Least_squares_extr  !
!=====================!
! Extrapolates smth. by polynomial fitting
Subroutine Least_squares_extr(Poly_Ord,NPoints,table,time_array,time)
Implicit none
Integer :: Poly_Ord, NPoints,s,n,i,j,p,q
Real(realk) :: time
Real(realk) :: table(NPoints+1), Time_array(NPoints+1),Red_time_array(NPoints+1)
Real(realk), dimension(NPoints,Poly_Ord+1) :: A
Real(realk), dimension(Poly_Ord+1,NPoints) :: AAA
Real(realk), dimension(Poly_Ord+1,Poly_Ord+1) :: AA, InvAA
Real(realk), dimension(NPoints) :: d   ! Final coefficients
Integer :: ErrorFlag
!
s = NPoints
n = Poly_Ord
!
! Calculating reduced time
!
Time_array(s+1) = Time
Red_time_array(1) = 1.D0
Red_time_array(s+1) = s + 1
Do j = 2, s
   Red_time_array(j) = 1.D0 + s*(Time_array(j) - Time_array(1))/(Time_array(s+1) - Time_Array(1))
Enddo
!
! Displacing reduced time to center it about zero
!
Do j = 1, (s+1)
   Red_time_array(j) = Red_time_array(j) - 0.5D0*(s+1)
Enddo
!
! Forming A matrix
!
Do i = 1, s
   Do j = 1, (n+1)
      A(i,j) = Red_time_array(i)**(j-1)
  Enddo
Enddo
!
! Forming d vector
!
If ( s .EQ. (n+1) ) then
   Call FINDInv(A,AAA,(n+1),ErrorFlag)
Else
   AA = MATMUL(TRANSPOSE(A),A)
   Call FINDInv(AA,InvAA,(n+1),ErrorFlag)
   AAA = MATMUL( InvAA, TRANSPOSE(A) )
Endif
d = 0.D0
Do p = 1, s
   Do q = 1, (n+1)
      d(p) = d(p) + AAA(q,p)*Red_time_array(s+1)**(q-1)
   Enddo
Enddo
!Write(*,*)'d vector',d
!
! Extrapolating
!
Do i = 1,s
   table(s+1) = table(s+1) + d(i)*table(i)
Enddo
!
End subroutine least_squares_extr
!
!===============================================!
!  Subroutine FINdInv                           !
!===============================================!
! Finds inverse matrices
!
Subroutine FINDInv(matrix, inverse, n, ErrorFlag)
Implicit none
Integer :: n
Integer :: ErrorFlag  !Return error status. -1 for error, 0 for normal
Double precision, dimension(n,n) :: matrix  !Input matrix
Double precision, dimension(n,n) :: inverse !Inverted matrix
Logical :: FLAG = .TRUE.
Integer :: i, j, k, l
Double precision :: m
Double precision, dimension(n,2*n) :: augmatrix !augmented matrix
!
!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
	   	        END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END Subroutine FINDinv
!
End module dyn_util
