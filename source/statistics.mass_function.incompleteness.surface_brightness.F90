!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implements calculations of incompleteness assuming a normal distribution of surface brightnesses.
  
  !# <massFunctionIncompleteness name="massFunctionIncompletenessSurfaceBrightness">
  !#  <description>Computes incompleteness assuming a normal distribution of surface brightnesses.</description>
  !# </massFunctionIncompleteness>

  type, extends(massFunctionIncompletenessClass) :: massFunctionIncompletenessSurfaceBrightness
     !% A class implementing incompleteness calculations assuming a normal distribution of surface brightnesses.
     private
     double precision :: surfaceBrightnessLimit       , surfaceBrightnessZeroPoint  , &
          &              surfaceBrightnessModelSlope  , surfaceBrightnessModelOffset, &
          &              surfaceBrightnessModelScatter
   contains
     procedure :: completeness => surfaceBrightnessSurfaceBrightnessness
  end type massFunctionIncompletenessSurfaceBrightness

  interface massFunctionIncompletenessSurfaceBrightness
     !% Constructors for the ``surface brightness'' incompleteness class.
     module procedure surfaceBrightnessDefaultConstructor
     module procedure surfaceBrightnessConstructor
  end interface massFunctionIncompletenessSurfaceBrightness

  ! Initialization state.
  logical          :: surfaceBrightnessInitialized=.false.

  ! Default values.
  double precision :: surfaceBrightnessDefaultLimit      , surfaceBrightnessDefaultZeroPoint  , &
       &              surfaceBrightnessDefaultModelSlope , surfaceBrightnessDefaultModelOffset, &
       &              surfaceBrightnessDefaultModelScatter

contains

  function surfaceBrightnessDefaultConstructor()
    !% Default constructor for the ``surface brightness'' incompleteness class.
    use Input_Parameters
    implicit none
    type(massFunctionIncompletenessSurfaceBrightness) :: surfaceBrightnessDefaultConstructor
    
    if (.not.surfaceBrightnessInitialized) then
       !$omp critical (massFunctionIncompletenessSurfaceBrightnessInitialize)
       if (.not.surfaceBrightnessInitialized) then
          !# <inputParameter>
          !#   <name>massFunctionIncompletenessSurfaceBrightnessLimit</name>
          !#   <cardinality>1</cardinality>
          !#   <description>Limiting surface brightness for mass function incompleteness calculations.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>surfaceBrightnessDefaultLimit</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>massFunctionIncompletenessSurfaceBrightnessZeroPoint</name>
          !#   <cardinality>1</cardinality>
          !#   <description>Mass zero point for the mass function incompleteness surface brightness model, i.e. $M_0$ in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>surfaceBrightnessDefaultZeroPoint</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>massFunctionIncompletenessSurfaceBrightnessModelSlope</name>
          !#   <cardinality>1</cardinality>
          !#   <description>Slope of mass function incompleteness surface brightness model, i.e. $\alpha$ in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>surfaceBrightnessDefaultModelSlope</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>massFunctionIncompletenessSurfaceBrightnessModelOffset</name>
          !#   <cardinality>1</cardinality>
          !#   <description>Offset in the mass function incompleteness surface brightness model, i.e. $beta$in $\mu(M) = \alpha \log_{10}(M/M_0)+\beta$.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>surfaceBrightnessDefaultModelOffset</variable>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>massFunctionIncompletenessSurfaceBrightnessModelScatter</name>
          !#   <cardinality>1</cardinality>
          !#   <description>Scatter in the mass function incompleteness surface brightness model.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !#   <variable>surfaceBrightnessDefaultModelScatter</variable>
          !# </inputParameter>
          ! Record that we are now initialized.
          surfaceBrightnessInitialized=.true.
       end if
       !$omp end critical (massFunctionIncompletenessSurfaceBrightnessInitialize)
    end if
    surfaceBrightnessDefaultConstructor=surfaceBrightnessConstructor(                                      &
         &                                                           surfaceBrightnessDefaultLimit       , &
         &                                                           surfaceBrightnessDefaultZeroPoint   , &
         &                                                           surfaceBrightnessDefaultModelSlope  , &
         &                                                           surfaceBrightnessDefaultModelOffset , &
         &                                                           surfaceBrightnessDefaultModelScatter  &
         &                                                          )
    return
  end function surfaceBrightnessDefaultConstructor

  function surfaceBrightnessConstructor(surfaceBrightnessLimit,surfaceBrightnessZeroPoint,surfaceBrightnessModelSlope,surfaceBrightnessModelOffset,surfaceBrightnessModelScatter)
    !% Generic constructor for the ``surface brightness'' incompleteness class.
    implicit none
    type            (massFunctionIncompletenessSurfaceBrightness) :: surfaceBrightnessConstructor
    double precision                                              :: surfaceBrightnessLimit       , surfaceBrightnessZeroPoint  , &
         &                                                           surfaceBrightnessModelSlope  , surfaceBrightnessModelOffset, &
         &                                                           surfaceBrightnessModelScatter
    
    surfaceBrightnessConstructor%surfaceBrightnessLimit       =surfaceBrightnessLimit
    surfaceBrightnessConstructor%surfaceBrightnessZeroPoint   =surfaceBrightnessZeroPoint
    surfaceBrightnessConstructor%surfaceBrightnessModelSlope  =surfaceBrightnessModelSlope
    surfaceBrightnessConstructor%surfaceBrightnessModelOffset =surfaceBrightnessModelOffset
    surfaceBrightnessConstructor%surfaceBrightnessModelScatter=surfaceBrightnessModelScatter
   return
  end function surfaceBrightnessConstructor

  double precision function surfaceBrightnessSurfaceBrightnessness(self,mass)
    !% Return the completeness.
    use Error_Functions
    implicit none
    class           (massFunctionIncompletenessSurfaceBrightness), intent(inout) :: self
    double precision                                             , intent(in   ) :: mass
    double precision                                                             :: surfaceBrightnessMean, surfaceBrightnessLimitNormalized

    surfaceBrightnessMean                 =self%surfaceBrightnessModelSlope*log10(mass/self%surfaceBrightnessZeroPoint)+self%surfaceBrightnessModelOffset
    surfaceBrightnessLimitNormalized      =(self%surfaceBrightnessLimit-surfaceBrightnessMean)/self%surfaceBrightnessModelScatter    
    surfaceBrightnessSurfaceBrightnessness=0.5d0*(1.0d0+Error_Function(surfaceBrightnessLimitNormalized))
    return
  end function surfaceBrightnessSurfaceBrightnessness
