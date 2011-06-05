!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module of useful astronomical constants.

module Numerical_Constants_Astronomical
  !% Contains various useful astronomical constants.
  use FGSL
  use Numerical_Constants_Prefixes
  use Numerical_Constants_Atomic
  public
  
  ! Solar mass (in kg).
  double precision, parameter :: massSolar=FGSL_CONST_MKSA_SOLAR_MASS

  ! Solar radius (in m; Allen's Astrophysical Quantities, page 340).
  double precision, parameter :: radiusSolar=6.95508d8

  ! Solar luminosity (in W; Allen's Astrophysical Quantities, page 340).
  double precision, parameter :: luminositySolar=3.845d26

  ! Solar composition (Allen's Atrophysical Quantities, page 28).
  double precision, parameter :: hydrogenByMassSolar=0.707d0
  double precision, parameter :: heliumByMassSolar  =0.274d0
  double precision, parameter :: metallicitySolar   =0.0188d0

  ! Primordial composition.
  double precision, parameter :: hydrogenByMassPrimordial=0.778d0
  double precision, parameter :: heliumByMassPrimordial  =0.222d0
  double precision, parameter :: metallicityPrimordial   =5.36d-10
  double precision, parameter :: meanAtomicMassPrimordial=1.0d0/(2.0d0*hydrogenByMassPrimordial/atomicMassHydrogen+3.0d0&
       &*heliumByMassPrimordial/atomicMassHelium)

  ! Megaparsec (in m).
  double precision, parameter :: parsec    =FGSL_CONST_MKSA_PARSEC
  double precision, parameter :: megaParsec=mega*parsec

  ! Years and related quantities (in s).
  double precision, parameter :: year=31558149.8d0 ! Sidereal year.
  double precision, parameter :: gigaYear=giga*year

end module Numerical_Constants_Astronomical
