!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_xcfun_interface

   use xcfun
   use fde_max_block_length

   implicit none

   public parse_functional
   public get_functional_derv
   public fun_is_lda
   public fun_is_gga
   public fun_is_tau_mgga
   public xcfun_set_functional
   public xc_get_type

   private

   public functional
   type functional
      real(8)              :: w(XC_NR_PARAMS)
      integer              :: id
      logical              :: xcfun_id_is_set
      real(8), allocatable :: input(:, :)
      real(8), allocatable :: output(:, :)
   end type

   integer, parameter, public :: d0000000 =  0
   integer, parameter, public :: d1000000 =  1
   integer, parameter, public :: d0010000 =  2
   integer, parameter, public :: d0000100 =  3
   integer, parameter, public :: d2000000 =  4
   integer, parameter, public :: d0200000 =  5
   integer, parameter, public :: d0020000 =  6
   integer, parameter, public :: d0002000 =  7
   integer, parameter, public :: d1010000 =  8
   integer, parameter, public :: d0101000 =  9
   integer, parameter, public :: d1000100 = 10
   integer, parameter, public :: d0010100 = 11
   integer, parameter, public :: d3000000 = 12
   integer, parameter, public :: d0030000 = 13
   integer, parameter, public :: d1200000 = 14
   integer, parameter, public :: d1020000 = 15
   integer, parameter, public :: d1002000 = 16
   integer, parameter, public :: d2010000 = 17
   integer, parameter, public :: d0210000 = 18
   integer, parameter, public :: d1101000 = 19
   integer, parameter, public :: d0111000 = 20
   integer, parameter, public :: d0012000 = 21

   integer, parameter, public :: d0000010 = 22
   integer, parameter, public :: d0000020 = 23
   integer, parameter, public :: d1000010 = 24
   integer, parameter, public :: d0010010 = 25
   integer, parameter, public :: d0100001 = 26
   integer, parameter, public :: d0001001 = 27
   integer, parameter, public :: d0000002 = 28

   integer, parameter, public :: nr_nonzero_derv = 28

contains

  function get_weight(word, functional)

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: word
    character(*), intent(in)  :: functional
    real(8)                   :: get_weight
!   ----------------------------------------------------------------------------
    integer                   :: i
!   ----------------------------------------------------------------------------

    get_weight = 0.0d0

    if (word_contains(word, '=')) then

      i = index(word, '=')

      if (word(1:i-1) == functional) then
        read(word(i+1:len(word)), *) get_weight
      end if

    else

      if (word == functional) then
        get_weight = 1.0d0
      end if

    end if

  end function

#ifndef PRG_DIRAC
   subroutine parse_functional(line,     &
                               f,        &
                               hfx_out,  &
                               mu_out,   &
                               beta_out, &
                               set_hf_exchange_factor)
!     --------------------------------------------------------------------------
      character(80)              :: line
      type(functional)           :: f
      real(8), intent(out)       :: hfx_out
      real(8), intent(out)       :: mu_out
      real(8), intent(out)       :: beta_out
      logical, intent(in)        :: set_hf_exchange_factor
!     --------------------------------------------------------------------------
      character(80), allocatable :: word_array(:)
      integer                    :: i, nr_words
      real(8)                    :: w, h
      integer                    :: u
!     --------------------------------------------------------------------------

      write(*, '(a)') '* Using the automatic differentiation xc functional:'
      write(*, '(a)') '   (weight: functional)'

   end subroutine parse_functional
#else
   subroutine parse_functional(line,     &
                               f,        &
                               hfx_out,  &
                               mu_out,   &
                               beta_out, &
                               set_hf_exchange_factor)

!     radovan: line is made lowercase, then cut into words
!              if word contains "=" the weight is read
!              otherwise it is set to 1.0
!
!              equivalent functionals:
!              b3lyp
!              slaterx=0.8 bcorrx=0.72 hfx=0.2 vwnc=0.19 lypc=0.81

!     --------------------------------------------------------------------------
      character(80)              :: line
      type(functional)           :: f
      real(8), intent(out)       :: hfx_out
      real(8), intent(out)       :: mu_out
      real(8), intent(out)       :: beta_out
      logical, intent(in)        :: set_hf_exchange_factor
!     --------------------------------------------------------------------------
      character(80), allocatable :: word_array(:)
      integer                    :: i, nr_words
      real(8)                    :: w, h
      integer                    :: u
!     --------------------------------------------------------------------------

      write(*, '(a)') '* Using the automatic differentiation xc functional:'
      write(*, '(a)') '   (weight: functional)'

      hfx_out  = 0.0d0
      mu_out   = 0.0d0
      beta_out = 0.0d0

      f%w = 0.0d0

      line = lowercase(line)

      nr_words = word_count(line)

      allocate(word_array(nr_words))
      read(line, *) (word_array(i), i = 1, nr_words)

      do i = 1, nr_words


!       composite functionals
!       =====================

        w = get_weight(word_array(i), 'lda')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) + w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) + w
        end if

        w = get_weight(word_array(i), 'blyp')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKEX) = f%w(XC_BECKEX) + w
           f%w(XC_LYPC  ) = f%w(XC_LYPC  ) + w
        end if

        w = get_weight(word_array(i), 'b3lyp')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX   ) = f%w(XC_SLATERX   ) + 0.80d0*w
           f%w(XC_BECKECORRX) = f%w(XC_BECKECORRX) + 0.72d0*w
           f%w(XC_VWN5C     ) = f%w(XC_VWN5C     ) + 0.19d0*w
           f%w(XC_LYPC      ) = f%w(XC_LYPC      ) + 0.81d0*w
           hfx_out = hfx_out + 0.2d0*w
        end if

        w = get_weight(word_array(i), 'bhandh')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX   ) = f%w(XC_SLATERX   ) + 0.5d0*w
           f%w(XC_LYPC      ) = f%w(XC_LYPC      ) +       w
           hfx_out = hfx_out + 0.5d0*w
        end if

        w = get_weight(word_array(i), 'pp86')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PW86X) = f%w(XC_PW86X) + w
           f%w(XC_P86C)  = f%w(XC_P86C)  + w
        end if

        w = get_weight(word_array(i), 'm05')
        if (dabs(w) > tiny(0.0d0)) then
           h = 0.28d0
!radovan:  notice the difference of the x weight in m05 and m06
!          it seems that the hfx modification is already absorbed
!          in the parameters of m06x
!          this is not the case for m05x
           f%w(XC_M05X) = f%w(XC_M05X) + w*(1.0d0 - h)
           f%w(XC_M05C) = f%w(XC_M05C) + w
           hfx_out = hfx_out + h*w
        end if

        w = get_weight(word_array(i), 'm06')
        if (dabs(w) > tiny(0.0d0)) then
           h = 0.27d0
           f%w(XC_M06X) = f%w(XC_M06X) + w
           f%w(XC_M06C) = f%w(XC_M06C) + w
           hfx_out = hfx_out + h*w
        end if

        w = get_weight(word_array(i), 'm06l')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_M06LX) = f%w(XC_M06LX) + w
           f%w(XC_M06LC) = f%w(XC_M06LC) + w
        end if

        w = get_weight(word_array(i), 'm06-l')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_M06LX) = f%w(XC_M06LX) + w
           f%w(XC_M06LC) = f%w(XC_M06LC) + w
        end if

        w = get_weight(word_array(i), 'pbe')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbe'
           f%w(XC_PBEX) = f%w(XC_PBEX) + w
           f%w(XC_PBEC) = f%w(XC_PBEC) + w
        end if

        w = get_weight(word_array(i), 'pbe0')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PBEX) = f%w(XC_PBEX) + 0.75d0*w
           f%w(XC_PBEC) = f%w(XC_PBEC) +        w
           hfx_out = hfx_out + 0.25d0*w
        end if

        w = get_weight(word_array(i), 'kt1')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) +            w
           f%w(XC_KTX    ) = f%w(XC_KTX    ) -    0.006d0*w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) +            w
        end if

        w = get_weight(word_array(i), 'kt2')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) +  1.07173d0*w
           f%w(XC_KTX    ) = f%w(XC_KTX    ) -    0.006d0*w
           f%w(XC_VWN5C  ) = f%w(XC_VWN5C  ) + 0.576727d0*w
        end if


!       "atomic" functionals
!       ====================

        w = get_weight(word_array(i), 'pw86x')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_PW86X) = f%w(XC_PW86X) + w
        end if

        w = get_weight(word_array(i), 'p86c')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_P86C) = f%w(XC_P86C) + w
        end if

        w = get_weight(word_array(i), 'kin_pw91')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': kin_pw91'
          f%w(XC_PW91K) = f%w(XC_PW91K) + w
        end if

        w = get_weight(word_array(i), 'kin_tf')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': kin_tf'
          f%w(XC_TFK) = f%w(XC_TFK) + w
        end if

        w = get_weight(word_array(i), 'ktx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_KTX) = f%w(XC_KTX) + w
        end if

        w = get_weight(word_array(i), 'pbex')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbex'
           f%w(XC_PBEX) = f%w(XC_PBEX) + w
        end if

        w = get_weight(word_array(i), 'pbec')
        if (dabs(w) > tiny(0.0d0)) then
          write(*, '(f10.6, a)') w, ': pbec'
           f%w(XC_PBEC) = f%w(XC_PBEC) + w
        end if

        w = get_weight(word_array(i), 'lypc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LYPC) = f%w(XC_LYPC) + w
        end if

        w = get_weight(word_array(i), 'hfx')
        if (dabs(w) > tiny(0.0d0)) then
           hfx_out = hfx_out + w
        end if

        w = get_weight(word_array(i), 'ldaerfx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFX) = f%w(XC_LDAERFX) + w
        end if

        w = get_weight(word_array(i), 'ldaerfc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFC) = f%w(XC_LDAERFC) + w
        end if

        w = get_weight(word_array(i), 'ldaerfc_jt')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_LDAERFC_JT) = f%w(XC_LDAERFC_JT) + w
        end if

        w = get_weight(word_array(i), 'mu')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_RANGESEP_MU) = w
           mu_out = mu_out + w
        end if

        w = get_weight(word_array(i), 'beta')
        if (dabs(w) > tiny(0.0d0)) then
           beta_out = beta_out + w
        end if

        w = get_weight(word_array(i), 'slaterx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_SLATERX) = f%w(XC_SLATERX) + w
        end if

        w = get_weight(word_array(i), 'vwnc')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_VWN5C) = f%w(XC_VWN5C) + w
        end if

        w = get_weight(word_array(i), 'beckex')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKEX) = f%w(XC_BECKEX) + w
        end if

        w = get_weight(word_array(i), 'bcorrx')
        if (dabs(w) > tiny(0.0d0)) then
           f%w(XC_BECKECORRX) = f%w(XC_BECKECORRX) + w
        end if

      end do

      deallocate(word_array)

      write(*, '(a)')

!     check that not all weights are zero
      h = 0.0d0
      do i = 1, size(f%w)
         h = h + dabs(f%w(i))
      end do
      if (h < tiny(0.0d0)) then
             print *, 'XCFun functional not recognized from input'
             stop
      end if

   end subroutine
#endif

  function word_count(s)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: s
    integer                  :: word_count
!   ----------------------------------------------------------------------------
    integer                  :: i
    logical                  :: is_blank
!   ----------------------------------------------------------------------------

    word_count = 0

    if (len(s) <= 0) return

    is_blank = .true.

    do i = 1, len(s)
      if (s(i:i) == ' ') then
        is_blank = .true.
      else if (is_blank) then
        word_count = word_count + 1
        is_blank = .false.
      end if
    end do

  end function

  function word_contains(word, substring)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: word
    character(*), intent(in) :: substring
!   ----------------------------------------------------------------------------
    logical                  :: word_contains
!   ----------------------------------------------------------------------------

    word_contains = .false.
    if (index(word, substring) > 0) then
      word_contains = .true.
    end if

  end function

  function lowercase(s)

!   ----------------------------------------------------------------------------
    character(*), intent(in) :: s
    character(len(s))        :: lowercase
!   ----------------------------------------------------------------------------
    integer                  :: off, i, ia
!   ----------------------------------------------------------------------------

    lowercase = s

    off = iachar('a') - iachar('A')

    do i = 1, len(s)
      ia = iachar(s(i:i))
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
        lowercase(i:i) = achar(ia + off)
      end if
    enddo

  end function

  subroutine get_functional_derv(f,            &
                                 f_alda,       &
                                 f_xalda,      &
                                 order,        &
                                 block_length, &
                                 w,            &
                                 r_0,          &
                                 z_0,          &
                                 n,            &
                                 derv,         &
                                 alda_real,    &
                                 alda_imag,    &
                                 xalda,        &
                                 tau)

!   ----------------------------------------------------------------------------
    type(functional)              :: f
    type(functional)              :: f_alda, f_xalda
    integer,          intent(in)  :: order
    integer,          intent(in)  :: block_length
    real(8),          intent(in)  :: w(*)
    real(8),          intent(in)  :: r_0(*)
    real(8),          intent(in)  :: z_0(*)
    integer,          intent(in)  :: n
    real(8),          intent(out) :: derv(max_block_length, 0:n)
    logical, optional, intent(in) :: alda_real
    logical, optional, intent(in) :: alda_imag
    logical, optional, intent(in) :: xalda
    real(8), optional, intent(in) :: tau(*)
!   ----------------------------------------------------------------------------
    integer                       :: k, l, id
#ifdef PRG_DIRAC
    real(8), external             :: dftenergy
#endif
    logical                       :: do_alda_real
    logical                       :: do_alda_imag
    logical                       :: do_xalda
!   ----------------------------------------------------------------------------

    do_alda_real = .false.
    do_alda_imag = .false.
    do_xalda     = .false.

    if (present(alda_real)) do_alda_real = alda_real
    if (present(alda_imag)) do_alda_imag = alda_imag
    if (present(xalda))     do_xalda     = xalda

!   if less then linear response ignore alda
    if (order < 2) then
       do_alda_real = .false.
       do_alda_imag = .false.
    end if

    derv = 0.0d0

      f%input = 0.0d0
      do k = 1, block_length
         f%input(1, k) = r_0(k)
      end do
      if (fun_is_gga(f)) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
      end if
      if (fun_is_tau_mgga(f)) then
         do k = 1, block_length
            f%input(3, k) = z_0(k)
         end do
         do k = 1, block_length
            f%input(6, k) = tau(k)
         end do
      end if

      call xc_eval(f%id,         &
                   order,        &
                   block_length, &
                   f%input,      &
                   f%output)

      if (do_alda_real .or. do_alda_imag) then
         f_alda%input = 0.0d0
         do k = 1, block_length
            f_alda%input(1, k) = r_0(k)
         end do
         id = f_alda%id
         if (do_xalda) then
            id = f_xalda%id
         end if
         call xc_eval(id,           &
                      order,        &
                      block_length, &
                      f_alda%input, &
                      f_alda%output)
      end if

      do k = 1, block_length

         derv(k, d0000000) = f%output(1, k)

         if (order == 0) then
               call fde_quit('get_functional_derv: order 0; probably programming error')
         end if

         if (order > 0) then
               if (fun_is_lda(f)) then
                  derv(k, d1000000) = f%output(XC_D10, k)
               else if (fun_is_gga(f)) then
                  derv(k, d1000000) = f%output(XC_D10000, k)
                  derv(k, d0010000) = f%output(XC_D00100, k)
               else if (fun_is_tau_mgga(f)) then
                  derv(k, d1000000) = f%output(XC_D1000000, k)
                  derv(k, d0010000) = f%output(XC_D0010000, k)
                  derv(k, d0000010) = f%output(XC_D0000010, k)
               end if
         end if

         if (order > 1) then
!              order 2 - spin unpolarized
               if (do_alda_real) then
                  derv(k, d2000000) = f_alda%output(XC_D20, k)
               else
                  if (fun_is_lda(f)) then
                     derv(k, d2000000) = f%output(XC_D20, k)
                  else if (fun_is_gga(f)) then
                     derv(k, d2000000) = f%output(XC_D20000, k)
                     derv(k, d1010000) = f%output(XC_D10100, k)
                     derv(k, d0020000) = f%output(XC_D00200, k)
                     derv(k, d0010000) = f%output(XC_D00100, k)
                  else if (fun_is_tau_mgga(f)) then
                     derv(k, d0010000) = f%output(XC_D0010000, k)
                     derv(k, d2000000) = f%output(XC_D2000000, k)
                     derv(k, d1010000) = f%output(XC_D1010000, k)
                     derv(k, d0020000) = f%output(XC_D0020000, k)

                     derv(k, d0000020) = f%output(XC_D0000020, k)
                     derv(k, d0010010) = f%output(XC_D0010010, k)
                     derv(k, d1000010) = f%output(XC_D1000010, k)
                  end if
               end if
!              order 2 - spin polarized
               if (do_alda_imag) then
                  derv(k, d0200000) = f_alda%output(XC_D02, k)
               else
                  if (fun_is_lda(f)) then
                     derv(k, d0200000) = f%output(XC_D02, k)
                  else if (fun_is_gga(f)) then
                     derv(k, d0200000) = f%output(XC_D02000, k)
                     derv(k, d0101000) = f%output(XC_D01010, k)
                     derv(k, d0002000) = f%output(XC_D00020, k)
                     derv(k, d0000100) = f%output(XC_D00001, k)
                  else if (fun_is_tau_mgga(f)) then
                     derv(k, d0000100) = f%output(XC_D0000100, k)
                     derv(k, d0200000) = f%output(XC_D0200000, k)
                     derv(k, d0101000) = f%output(XC_D0101000, k)
                     derv(k, d0002000) = f%output(XC_D0002000, k)

                     derv(k, d0000002) = f%output(XC_D0000002, k)
                     derv(k, d0001001) = f%output(XC_D0001001, k)
                     derv(k, d0100001) = f%output(XC_D0100001, k)
                  end if
               end if
         end if

         if (order > 2) then
               if (do_alda_real .or. do_alda_imag) then
                  call fde_quit('get_functional_derv: xcfun alda qr not implemented')
               end if
               if (fun_is_lda(f)) then
                  derv(k, d2000000) = f%output(XC_D20, k)
                  derv(k, d0200000) = f%output(XC_D02, k)
                  derv(k, d1200000) = f%output(XC_D12, k)
                  derv(k, d3000000) = f%output(XC_D30, k)
               else if (fun_is_gga(f)) then
                  derv(k, d0010000) = f%output(XC_D00100, k)
                  derv(k, d0000100) = f%output(XC_D00001, k)
                  derv(k, d2000000) = f%output(XC_D20000, k)
                  derv(k, d1010000) = f%output(XC_D10100, k)
                  derv(k, d1000100) = f%output(XC_D10001, k)
                  derv(k, d0200000) = f%output(XC_D02000, k)
                  derv(k, d0101000) = f%output(XC_D01010, k)
                  derv(k, d0020000) = f%output(XC_D00200, k)
                  derv(k, d0010100) = f%output(XC_D00101, k)
                  derv(k, d0002000) = f%output(XC_D00020, k)
                  derv(k, d3000000) = f%output(XC_D30000, k)
                  derv(k, d2010000) = f%output(XC_D20100, k)
                  derv(k, d1200000) = f%output(XC_D12000, k)
                  derv(k, d1101000) = f%output(XC_D11010, k)
                  derv(k, d1020000) = f%output(XC_D10200, k)
                  derv(k, d1002000) = f%output(XC_D10020, k)
                  derv(k, d0210000) = f%output(XC_D02100, k)
                  derv(k, d0111000) = f%output(XC_D01110, k)
                  derv(k, d0030000) = f%output(XC_D00300, k)
                  derv(k, d0012000) = f%output(XC_D00120, k)
               else if (fun_is_tau_mgga(f)) then
                  call fde_quit('get_functional_derv: tau mgga qr not implemented')
               end if
         end if

         if (order > 3) then
               call fde_quit('get_functional_derv: order too high')
         end if

      end do

!   multiply by weight
    do l = 0, n
       do k = 1, block_length
          derv(k, l) = w(k)*derv(k, l)
       end do
    end do

  end subroutine

   function fun_is_lda(f)
      type(functional) :: f
      logical :: fun_is_lda
      fun_is_lda = (xc_get_type(f%id) == 0)
   end function

   function fun_is_gga(f)
      type(functional) :: f
      logical :: fun_is_gga
      fun_is_gga = (xc_get_type(f%id) == 1)
   end function

   function fun_is_tau_mgga(f)
      type(functional) :: f
      logical :: fun_is_tau_mgga
      fun_is_tau_mgga = (xc_get_type(f%id) == 2)
   end function

   subroutine xcfun_set_functional(f,              &
                                   response_order, &
                                   parallel_xc)

!     --------------------------------------------------------------------------
      type(functional), intent(inout) :: f
      integer,          intent(in)    :: response_order
      logical,          intent(in)    :: parallel_xc
!     --------------------------------------------------------------------------
      integer                         :: i, l
!     --------------------------------------------------------------------------

      if (.not. f%xcfun_id_is_set) then
         f%id = xc_new_functional()
         call xc_set_mode(f%id, XC_VARS_NS)

         do i = 1, XC_NR_PARAMS
            if (dabs(f%w(i)) > tiny(0.0d0)) then
               call xc_set_param(f%id, i, f%w(i))
            end if
         end do
         f%xcfun_id_is_set = .true.
      end if

      if (allocated(f%input)) then
         l = size(f%input)/max_block_length
      else
         l = 0
      end if
      if (xc_input_length(f%id) > l) then
         l = xc_input_length(f%id)
         if (allocated(f%input)) then
            deallocate(f%input)
         end if
         allocate(f%input(l, max_block_length))
      end if

      if (allocated(f%output)) then
         l = size(f%output)/max_block_length
      else
         l = 0
      end if
      if (xc_output_length(f%id, 1 + response_order) > l) then
         l = xc_output_length(f%id, 1 + response_order)
         if (allocated(f%output)) then
            deallocate(f%output)
         end if
         allocate(f%output(l, max_block_length))
      end if

   end subroutine

end module
