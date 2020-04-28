
/*
!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2020 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-rpbex.c:
   implementation of RPBEx functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001
   Z. Rinkevicius adapted for open shell systems: energy, first derivatives.
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.

   Derivatives in this file generated with SymPy using xcdiff by Olav Vahtras

    ~~~
    from math import pi
    from textwrap import indent

    import sympy
    import xcdiff


    THIS_FILE = f"~~~\n{open(__file__).read()}~~~"


    def rpbex():
        ra, rb, ga, gb = sympy.symbols('dp->rhoa, dp->rhob, dp->grada, dp->gradb')

        xa = ga/pow(ra, 4/3)
        xb = gb/pow(rb, 4/3)

        R = 0.804
        d = 0.066725
        mu = d*pi**2/3
        Sa = xa/(2*(6*pi**2)**(1/3))
        Sb = xb/(2*(6*pi**2)**(1/3))

        def F(S):
            return 1 + R*(1-sympy.exp(-mu*S**2/R))

        def Ea(n):
            return -3/(4*pi)*(3*pi**2)**(1/3)*n**(4/3)*F(Sa)

        def Eb(n):
            return -3/(4*pi)*(3*pi**2)**(1/3)*n**(4/3)*F(Sb)

        func = xcdiff.GGAFunctional(
            'RPBEx',
            ra,
            rb,
            ga,
            gb,
            0.5*Ea(2*ra),
            0.5*Eb(2*rb),
            info=indent(THIS_FILE, '    ')
        )
        return func


    func = rpbex()

    with open('dft/fun-rpbex.c', 'w') as f:
        print(func, file=f)
    ~~~
*/

#include <math.h>
#include <stdio.h>
#include "general.h"

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static integer rpbex_isgga(void) { return 1; }
static integer rpbex_read(const char* conf_line);
static real rpbex_energy(const FunDensProp* dp);
static void rpbex_first(FunFirstFuncDrv *ds,   real fac, const FunDensProp*);
static void rpbex_second(FunSecondFuncDrv *ds, real fac, const FunDensProp*);
static void rpbex_third(FunThirdFuncDrv *ds,   real fac, const FunDensProp*);
static void rpbex_fourth(FunFourthFuncDrv *ds, real fac, const FunDensProp*);

Functional RPBExFunctional = {
  "RPBEx",       /* name */
  rpbex_isgga,   /* gga-corrected */
   3,
  rpbex_read,
  NULL,
  rpbex_energy,
  rpbex_first,
  rpbex_second,
  rpbex_third,
  rpbex_fourth
};

/* IMPLEMENTATION PART */
static integer
rpbex_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
rpbex_energy(const FunDensProp* dp)
{
  return -0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))+-0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)));
}

static void
rpbex_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += (0.0089633469276457472*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df0010 += (-0.00672251019573431*dp->grada*pow(dp->rhoa, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df0100 += (0.0089633469276457472*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0001 += (-0.00672251019573431*dp->gradb*pow(dp->rhob, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
}

static void
rpbex_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += (0.0089633469276457472*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df0010 += (-0.00672251019573431*dp->grada*pow(dp->rhoa, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df2000 += (0.023902258473721991*pow(dp->grada, 2)*pow(dp->rhoa, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhoa, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df1010 += (dp->grada*(-8.0541040835315503e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.0) + 0.0089633469276457454*pow(dp->rhoa, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0020 += (0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df0100 += (0.0089633469276457472*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0001 += (-0.00672251019573431*dp->gradb*pow(dp->rhob, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0200 += (0.023902258473721991*pow(dp->gradb, 2)*pow(dp->rhob, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhob, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0101 += (dp->gradb*(-8.0541040835315503e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.0) + 0.0089633469276457454*pow(dp->rhob, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0002 += (0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
}

static void
rpbex_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += (0.0089633469276457472*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df0010 += (-0.00672251019573431*dp->grada*pow(dp->rhoa, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df2000 += (0.023902258473721991*pow(dp->grada, 2)*pow(dp->rhoa, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhoa, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df1010 += (dp->grada*(-8.0541040835315503e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.0) + 0.0089633469276457454*pow(dp->rhoa, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0020 += (0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df3000 += (0.01195112923686099*pow(dp->grada, 2)*pow(dp->rhoa, -4.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 3.7221029453963999*pow(dp->grada, 2)*pow(dp->rhoa, 0.33333333333333326)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(1.3826534867507642e-6*pow(dp->grada, 4)*pow(dp->rhoa, -11.0) - 0.0012694636513253738*pow(dp->grada, 2)*pow(dp->rhoa, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhoa, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.27571132928862213*pow(dp->rhoa, -1.6666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df2010 += (dp->grada*(0.00021477610889417465*pow(dp->grada, 2)*pow(dp->rhoa, -6.0) - 0.00672251019573431*pow(dp->grada, 2)*pow(dp->rhoa, -1.3333333333333333)*(0.00014353953542801603*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.043929594917815104*pow(dp->rhoa, -4.6666666666666661)) - 0.020914476164506739*pow(dp->rhoa, -3.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df1020 += ((0.011148441452295705*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) + 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) - 0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(0.00034621735945237478*pow(dp->grada, 2)*pow(dp->rhoa, -6.333333333333333) - 0.019265124171230916*pow(dp->rhoa, -3.6666666666666665)))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0030 += (-0.93052573634909996*dp->grada*pow(dp->rhoa, 1.3333333333333333)*(5.8330693972297879e-7*pow(dp->grada, 2)*pow(dp->rhoa, -8.0) - 0.00019474726469196081*pow(dp->rhoa, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df0100 += (0.0089633469276457472*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0001 += (-0.00672251019573431*dp->gradb*pow(dp->rhob, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0200 += (0.023902258473721991*pow(dp->gradb, 2)*pow(dp->rhob, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhob, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0101 += (dp->gradb*(-8.0541040835315503e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.0) + 0.0089633469276457454*pow(dp->rhob, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0002 += (0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0300 += (0.01195112923686099*pow(dp->gradb, 2)*pow(dp->rhob, -4.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 3.7221029453963999*pow(dp->gradb, 2)*pow(dp->rhob, 0.33333333333333326)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(1.3826534867507642e-6*pow(dp->gradb, 4)*pow(dp->rhob, -11.0) - 0.0012694636513253738*pow(dp->gradb, 2)*pow(dp->rhob, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhob, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.27571132928862213*pow(dp->rhob, -1.6666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0201 += (dp->gradb*(0.00021477610889417465*pow(dp->gradb, 2)*pow(dp->rhob, -6.0) - 0.00672251019573431*pow(dp->gradb, 2)*pow(dp->rhob, -1.3333333333333333)*(0.00014353953542801603*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.043929594917815104*pow(dp->rhob, -4.6666666666666661)) - 0.020914476164506739*pow(dp->rhob, -3.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0102 += ((0.011148441452295705*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) + 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) - 0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(0.00034621735945237478*pow(dp->gradb, 2)*pow(dp->rhob, -6.333333333333333) - 0.019265124171230916*pow(dp->rhob, -3.6666666666666665)))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0003 += (-0.93052573634909996*dp->gradb*pow(dp->rhob, 1.3333333333333333)*(5.8330693972297879e-7*pow(dp->gradb, 2)*pow(dp->rhob, -8.0) - 0.00019474726469196081*pow(dp->rhob, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
}

static void
rpbex_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
  ds->df1000 += (0.0089633469276457472*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df0010 += (-0.00672251019573431*dp->grada*pow(dp->rhoa, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df2000 += (0.023902258473721991*pow(dp->grada, 2)*pow(dp->rhoa, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhoa, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df1010 += (dp->grada*(-8.0541040835315503e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.0) + 0.0089633469276457454*pow(dp->rhoa, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0020 += (0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df3000 += (0.01195112923686099*pow(dp->grada, 2)*pow(dp->rhoa, -4.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 3.7221029453963999*pow(dp->grada, 2)*pow(dp->rhoa, 0.33333333333333326)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(1.3826534867507642e-6*pow(dp->grada, 4)*pow(dp->rhoa, -11.0) - 0.0012694636513253738*pow(dp->grada, 2)*pow(dp->rhoa, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhoa, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 0.27571132928862213*pow(dp->rhoa, -1.6666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df2010 += (dp->grada*(0.00021477610889417465*pow(dp->grada, 2)*pow(dp->rhoa, -6.0) - 0.00672251019573431*pow(dp->grada, 2)*pow(dp->rhoa, -1.3333333333333333)*(0.00014353953542801603*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.043929594917815104*pow(dp->rhoa, -4.6666666666666661)) - 0.020914476164506739*pow(dp->rhoa, -3.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df1020 += ((0.011148441452295705*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) + 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) - 0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(0.00034621735945237478*pow(dp->grada, 2)*pow(dp->rhoa, -6.333333333333333) - 0.019265124171230916*pow(dp->rhoa, -3.6666666666666665)))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0030 += (-0.93052573634909996*dp->grada*pow(dp->rhoa, 1.3333333333333333)*(5.8330693972297879e-7*pow(dp->grada, 2)*pow(dp->rhoa, -8.0) - 0.00019474726469196081*pow(dp->rhoa, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df4000 += (-0.010623225988320882*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333)*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 2.481401963597599*pow(dp->grada, 2)*pow(dp->rhoa, -0.66666666666666674)*(0.00011540578648412491*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.035319394313923343*pow(dp->rhoa, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) + 4.9628039271951998*pow(dp->grada, 2)*pow(dp->rhoa, 0.33333333333333326)*(1.3826534867507642e-6*pow(dp->grada, 4)*pow(dp->rhoa, -11.0) - 0.0012694636513253738*pow(dp->grada, 2)*pow(dp->rhoa, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhoa, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(-1.6565292977636101e-8*pow(dp->grada, 6)*pow(dp->rhoa, -14.666666666666666) + 1.5209188354258407e-5*pow(dp->grada, 4)*pow(dp->rhoa, -12.0) + 1.5209188354258405e-5*pow(dp->grada, 4)*pow(dp->rhoa, -11.999999999999998) - 0.012553584996439805*pow(dp->grada, 2)*pow(dp->rhoa, -9.3333333333333321) + 0.93400176074597263*pow(dp->rhoa, -6.6666666666666661))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)) - 0.45951888214770359*pow(dp->rhoa, -2.666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665))))*factor;
  ds->df3010 += (dp->grada*(-0.00075171638112961131*pow(dp->grada, 2)*pow(dp->rhoa, -7.0) + 0.026890040782937236*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*(0.00014353953542801603*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.043929594917815104*pow(dp->rhoa, -4.6666666666666661)) - 0.00672251019573431*pow(dp->grada, 2)*pow(dp->rhoa, -1.3333333333333333)*(1.7197182671029402e-6*pow(dp->grada, 4)*pow(dp->rhoa, -11.0) - 0.0015789348897081764*pow(dp->grada, 2)*pow(dp->rhoa, -8.3333333333333321) + 0.20500477628313712*pow(dp->rhoa, -5.6666666666666661)) + 0.069714920548355791*pow(dp->rhoa, -4.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df2020 += ((0.029729177206121876*pow(dp->grada, 2)*pow(dp->rhoa, -3.333333333333333)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) - 0.022296882904591409*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*(0.00034621735945237478*pow(dp->grada, 2)*pow(dp->rhoa, -6.333333333333333) - 0.019265124171230916*pow(dp->rhoa, -3.6666666666666665)) + 0.93052573634909996*pow(dp->grada, 2)*pow(dp->rhoa, 1.3333333333333333)*(0.00014353953542801603*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.043929594917815104*pow(dp->rhoa, -4.6666666666666661))*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) + 0.41356699393293317*pow(dp->rhoa, -0.66666666666666674)*(6.491575489732027e-5*pow(dp->grada, 2)*pow(dp->rhoa, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhoa, -2.6666666666666665)) - 2.4814019635975999*pow(dp->rhoa, 0.33333333333333326)*(0.00034621735945237478*pow(dp->grada, 2)*pow(dp->rhoa, -6.333333333333333) - 0.019265124171230916*pow(dp->rhoa, -3.6666666666666665)) + 0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(0.0021927099431983733*pow(dp->grada, 2)*pow(dp->rhoa, -7.333333333333333) - 0.070638788627846685*pow(dp->rhoa, -4.6666666666666661)))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df1030 += (dp->grada*(-0.011148441452295705*pow(dp->grada, 2)*pow(dp->rhoa, -2.333333333333333)*(5.8330693972297879e-7*pow(dp->grada, 2)*pow(dp->rhoa, -8.0) - 0.00019474726469196081*pow(dp->rhoa, -5.333333333333333)) - 1.2407009817988*pow(dp->rhoa, 0.33333333333333326)*(5.8330693972297879e-7*pow(dp->grada, 2)*pow(dp->rhoa, -8.0) - 0.00019474726469196081*pow(dp->rhoa, -5.333333333333333)) + 0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(4.6664555177838303e-6*pow(dp->grada, 2)*pow(dp->rhoa, -9.0) - 0.0010386520783571243*pow(dp->rhoa, -6.333333333333333)))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;
  ds->df0040 += (0.93052573634909996*pow(dp->rhoa, 1.3333333333333333)*(5.2413622312051744e-9*pow(dp->grada, 4)*pow(dp->rhoa, -10.666666666666666) - 3.4998416383378725e-6*pow(dp->grada, 2)*pow(dp->rhoa, -8.0) + 0.00019474726469196081*pow(dp->rhoa, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->grada, 2)*pow(dp->rhoa, -2.6666666666666665)))*factor;

  ds->df0100 += (0.0089633469276457472*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0001 += (-0.00672251019573431*dp->gradb*pow(dp->rhob, -1.3333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0200 += (0.023902258473721991*pow(dp->gradb, 2)*pow(dp->rhob, -3.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 0.41356699393293317*pow(dp->rhob, -0.66666666666666674)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0101 += (dp->gradb*(-8.0541040835315503e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.0) + 0.0089633469276457454*pow(dp->rhob, -2.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0002 += (0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0300 += (0.01195112923686099*pow(dp->gradb, 2)*pow(dp->rhob, -4.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 3.7221029453963999*pow(dp->gradb, 2)*pow(dp->rhob, 0.33333333333333326)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(1.3826534867507642e-6*pow(dp->gradb, 4)*pow(dp->rhob, -11.0) - 0.0012694636513253738*pow(dp->gradb, 2)*pow(dp->rhob, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhob, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 0.27571132928862213*pow(dp->rhob, -1.6666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0201 += (dp->gradb*(0.00021477610889417465*pow(dp->gradb, 2)*pow(dp->rhob, -6.0) - 0.00672251019573431*pow(dp->gradb, 2)*pow(dp->rhob, -1.3333333333333333)*(0.00014353953542801603*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.043929594917815104*pow(dp->rhob, -4.6666666666666661)) - 0.020914476164506739*pow(dp->rhob, -3.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0102 += ((0.011148441452295705*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) + 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) - 0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(0.00034621735945237478*pow(dp->gradb, 2)*pow(dp->rhob, -6.333333333333333) - 0.019265124171230916*pow(dp->rhob, -3.6666666666666665)))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0003 += (-0.93052573634909996*dp->gradb*pow(dp->rhob, 1.3333333333333333)*(5.8330693972297879e-7*pow(dp->gradb, 2)*pow(dp->rhob, -8.0) - 0.00019474726469196081*pow(dp->rhob, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;

  ds->df0400 += (-0.010623225988320882*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333)*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 2.481401963597599*pow(dp->gradb, 2)*pow(dp->rhob, -0.66666666666666674)*(0.00011540578648412491*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.035319394313923343*pow(dp->rhob, -4.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) + 4.9628039271951998*pow(dp->gradb, 2)*pow(dp->rhob, 0.33333333333333326)*(1.3826534867507642e-6*pow(dp->gradb, 4)*pow(dp->rhob, -11.0) - 0.0012694636513253738*pow(dp->gradb, 2)*pow(dp->rhob, -8.3333333333333321) + 0.16482384013164225*pow(dp->rhob, -5.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(-1.6565292977636101e-8*pow(dp->gradb, 6)*pow(dp->rhob, -14.666666666666666) + 1.5209188354258407e-5*pow(dp->gradb, 4)*pow(dp->rhob, -12.0) + 1.5209188354258405e-5*pow(dp->gradb, 4)*pow(dp->rhob, -11.999999999999998) - 0.012553584996439805*pow(dp->gradb, 2)*pow(dp->rhob, -9.3333333333333321) + 0.93400176074597263*pow(dp->rhob, -6.6666666666666661))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)) - 0.45951888214770359*pow(dp->rhob, -2.666666666666667)*(1.804 - 0.80400000000000005*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665))))*factor;
  ds->df0301 += (dp->gradb*(-0.00075171638112961131*pow(dp->gradb, 2)*pow(dp->rhob, -7.0) + 0.026890040782937236*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*(0.00014353953542801603*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.043929594917815104*pow(dp->rhob, -4.6666666666666661)) - 0.00672251019573431*pow(dp->gradb, 2)*pow(dp->rhob, -1.3333333333333333)*(1.7197182671029402e-6*pow(dp->gradb, 4)*pow(dp->rhob, -11.0) - 0.0015789348897081764*pow(dp->gradb, 2)*pow(dp->rhob, -8.3333333333333321) + 0.20500477628313712*pow(dp->rhob, -5.6666666666666661)) + 0.069714920548355791*pow(dp->rhob, -4.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0202 += ((0.029729177206121876*pow(dp->gradb, 2)*pow(dp->rhob, -3.333333333333333)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) - 0.022296882904591409*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*(0.00034621735945237478*pow(dp->gradb, 2)*pow(dp->rhob, -6.333333333333333) - 0.019265124171230916*pow(dp->rhob, -3.6666666666666665)) + 0.93052573634909996*pow(dp->gradb, 2)*pow(dp->rhob, 1.3333333333333333)*(0.00014353953542801603*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.043929594917815104*pow(dp->rhob, -4.6666666666666661))*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) + 0.41356699393293317*pow(dp->rhob, -0.66666666666666674)*(6.491575489732027e-5*pow(dp->gradb, 2)*pow(dp->rhob, -5.333333333333333) - 0.0072244215642115941*pow(dp->rhob, -2.6666666666666665)) - 2.4814019635975999*pow(dp->rhob, 0.33333333333333326)*(0.00034621735945237478*pow(dp->gradb, 2)*pow(dp->rhob, -6.333333333333333) - 0.019265124171230916*pow(dp->rhob, -3.6666666666666665)) + 0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(0.0021927099431983733*pow(dp->gradb, 2)*pow(dp->rhob, -7.333333333333333) - 0.070638788627846685*pow(dp->rhob, -4.6666666666666661)))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0103 += (dp->gradb*(-0.011148441452295705*pow(dp->gradb, 2)*pow(dp->rhob, -2.333333333333333)*(5.8330693972297879e-7*pow(dp->gradb, 2)*pow(dp->rhob, -8.0) - 0.00019474726469196081*pow(dp->rhob, -5.333333333333333)) - 1.2407009817988*pow(dp->rhob, 0.33333333333333326)*(5.8330693972297879e-7*pow(dp->gradb, 2)*pow(dp->rhob, -8.0) - 0.00019474726469196081*pow(dp->rhob, -5.333333333333333)) + 0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(4.6664555177838303e-6*pow(dp->gradb, 2)*pow(dp->rhob, -9.0) - 0.0010386520783571243*pow(dp->rhob, -6.333333333333333)))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
  ds->df0004 += (0.93052573634909996*pow(dp->rhob, 1.3333333333333333)*(5.2413622312051744e-9*pow(dp->gradb, 4)*pow(dp->rhob, -10.666666666666666) - 3.4998416383378725e-6*pow(dp->gradb, 2)*pow(dp->rhob, -8.0) + 0.00019474726469196081*pow(dp->rhob, -5.333333333333333))*exp(-0.0044927994802310906*pow(dp->gradb, 2)*pow(dp->rhob, -2.6666666666666665)))*factor;
}

