!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements an N-body dark matter halo mass error class using the model of \cite{trenti_how_2010}.

  !# <nbodyHaloMassError name="nbodyHaloMassErrorTrenti2010">
  !#  <description>An N-body dark matter halo mass error class using the model of \cite{trenti_how_2010}.</description>
  !# </nbodyHaloMassError>
  type, extends(nbodyHaloMassErrorPowerLaw) :: nbodyHaloMassErrorTrenti2010
     !% An N-body halo mass error class using the model of \cite{trenti_how_2010}.
     private
  end type nbodyHaloMassErrorTrenti2010

  interface nbodyHaloMassErrorTrenti2010
     !% Constructors for the {\normalfont \ttfamily trenti2010} N-body halo mass error class.
     module procedure nbodyHaloMassErrorTrenti2010Parameters
     module procedure nbodyHaloMassErrorTrenti2010Internal
  end interface nbodyHaloMassErrorTrenti2010

contains

  function nbodyHaloMassErrorTrenti2010Parameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily trenti2010} N-body halo mass error class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (nbodyHaloMassErrorTrenti2010)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: massParticle
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>massParticle</variable>
    !#   <description>The mass of the particle in the N-body simulation in which halos were found.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyHaloMassErrorTrenti2010(massParticle)
    return
  end function nbodyHaloMassErrorTrenti2010Parameters

  function nbodyHaloMassErrorTrenti2010Internal(massParticle) result(self)
    !% Internal constructor for the {\normalfont \ttfamily trenti2010} N-body halo mass error class. \cite{trenti_how_2010} report
    !% a normalization of the fractional error in particle number of 0.15 at $N=1000$ particles. Since this is based on
    !% comparisons of halos in simulations differing in number of particles by a factor $8$ this actually overestimates the
    !% normalization by a factor $\sqrt{5/4}$. Therefore, we use a normalization of $0.135$ here.
    implicit none
    type            (nbodyHaloMassErrorTrenti2010)                :: self
    double precision                              , intent(in   ) :: massParticle
    double precision                              , parameter     :: exponent               =-1.000d0/3.0d0
    double precision                              , parameter     :: normalization          =+0.135d0
    double precision                              , parameter     :: particleNumberReference=+1.000d3
    
    self%normalizationSquared          =(normalization*(powerLawMassReference/particleNumberReference/massParticle)**exponent)**2
    self%exponent                      =                                                                             exponent
    self%fractionalErrorHighMassSquared=+0.0d0
    return
  end function nbodyHaloMassErrorTrenti2010Internal
