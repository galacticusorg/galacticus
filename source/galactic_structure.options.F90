!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which provides various internal option codes for the galactic structure functions.

module Galactic_Structure_Options
  !% Provides various internal option codes for the galactic structure functions.
  implicit none
  public

  ! Values used to represent different mass types.
  integer,          parameter :: massTypeAll      =0
  integer,          parameter :: massTypeDark     =1
  integer,          parameter :: massTypeBaryonic =2
  integer,          parameter :: massTypeGalactic =3
  integer,          parameter :: massTypeGaseous  =4
  integer,          parameter :: massTypeStellar  =5
  integer,          parameter :: massTypeBlackHole=6

  ! Values used to represent different component types.
  integer,          parameter :: componentTypeAll      =0
  integer,          parameter :: componentTypeDisk     =1
  integer,          parameter :: componentTypeSpheroid =2
  integer,          parameter :: componentTypeHotHalo  =3
  integer,          parameter :: componentTypeDarkHalo =4
  integer,          parameter :: componentTypeBlackHole=5

  ! Coordinate system options.
  integer,          parameter :: coordinateSystemSpherical  =1
  integer,          parameter :: coordinateSystemCylindrical=2
  integer,          parameter :: coordinateSystemCartesian  =3

  ! Weighting options.
  integer,          parameter :: weightByMass      =0
  integer,          parameter :: weightByLuminosity=1
  integer,          parameter :: weightIndexNull   =0

  ! Suitably large value to represent infinite radius.
  double precision, parameter :: radiusLarge=1.0d10

end module Galactic_Structure_Options
