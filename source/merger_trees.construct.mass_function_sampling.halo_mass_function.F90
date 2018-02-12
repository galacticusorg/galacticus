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

  !% Implementation of a merger tree halo mass function sampling class in which the sampling rate is proportional to the halo mass function.

  use Halo_Mass_Functions

  !# <mergerTreeHaloMassFunctionSampling name="mergerTreeHaloMassFunctionSamplingHaloMassFunction" defaultThreadPrivate="yes">
  !#  <description>A merger tree halo mass function sampling class in which the sampling rate is proportional to the halo mass function.</description>
  !# </mergerTreeHaloMassFunctionSampling>
  type, extends(mergerTreeHaloMassFunctionSamplingClass) :: mergerTreeHaloMassFunctionSamplingHaloMassFunction
     !% Implementation of merger tree halo mass function sampling class in which the sampling rate is proportional to the halo mass function.
     private
     class           (haloMassFunctionClass), pointer :: haloMassFunction_
     double precision                                 :: abundanceMinimum, abundanceMaximum, &
          &                                              modifier1       , modifier2
   contains
     final     ::           haloMassFunctionDestructor
     procedure :: sample => haloMassFunctionSample
  end type mergerTreeHaloMassFunctionSamplingHaloMassFunction

  interface mergerTreeHaloMassFunctionSamplingHaloMassFunction
     !% Constructors for the {\normalfont \ttfamily haloMassFunction} merger tree halo mass function sampling class.
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface mergerTreeHaloMassFunctionSamplingHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the haloMassFunction merger tree halo mass function sampling class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (mergerTreeHaloMassFunctionSamplingHaloMassFunction)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (haloMassFunctionClass                             ), pointer       :: haloMassFunction_
    double precision                                                                    :: abundanceMinimum , abundanceMaximum, &
         &                                                                                 modifier1        , modifier2

    !# <inputParameter>
    !#   <name>abundanceMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-1.0d0</defaultValue>
    !#   <description>The abundance (in units of Mpc$^{-3}$) below which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>abundanceMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>-1.0d0</defaultValue>
    !#   <description>The abundance (in units of Mpc$^{-3}$) above which to truncate the halo mass function when sampling halo masses for tree construction. A negative value indicates no truncation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>modifier1</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Coefficient of the polynomial modifier applied to the halo mass function when sampling halo masses for tree construction.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>modifier2</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Coefficient of the polynomial modifier applied to the halo mass function when sampling halo masses for tree construction.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="parameters"/>
    self=mergerTreeHaloMassFunctionSamplingHaloMassFunction(abundanceMinimum,abundanceMaximum,modifier1,modifier2,haloMassFunction_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(abundanceMinimum,abundanceMaximum,modifier1,modifier2,haloMassFunction_) result(self)
    !% Internal constructor for the haloMassFunction merger tree halo mass function sampling class.
    implicit none
    type            (mergerTreeHaloMassFunctionSamplingHaloMassFunction)                        :: self
    double precision                                                    , intent(in   )         :: abundanceMinimum, abundanceMaximum, &
         &                                                                                         modifier1       , modifier2
    class           (haloMassFunctionClass                             ), intent(in   ), target :: haloMassFunction_
    !# <constructorAssign variables="abundanceMinimum, abundanceMaximum, modifier1, modifier2, *haloMassFunction_"/>
    
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloMassFunction} merger tree halo mass sampling class.
    implicit none
    type(mergerTreeHaloMassFunctionSamplingHaloMassFunction), intent(inout) :: self

    !# <objectDestructor name="self%haloMassFunction_"/>
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionSample(self,mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a volume-limited sampling.
    implicit none
    class           (mergerTreeHaloMassFunctionSamplingHaloMassFunction), intent(inout) :: self
    double precision                                                    , intent(in   ) :: mass                , massMaximum, &
         &                                                                                 massMinimum         , time
    double precision                                                    , parameter     :: massZeroPoint=1.0d13
    !GCC$ attributes unused :: massMinimum, massMaximum

    ! Construct sampling rate.
    haloMassFunctionSample=+                                         mass         &
         &                 *self%haloMassFunction_%differential(time,mass)        &
         &                 *10.0d0**(                                             &
         &                           +self%modifier1*log10(mass/massZeroPoint)    &
         &                           +self%modifier2*log10(mass/massZeroPoint)**2 &
         &                          )
    ! Limit sampling rate.
    if (self%abundanceMinimum > 0.0d0)  &
         & haloMassFunctionSample       &
         & =max(                        &
         &      haloMassFunctionSample, &
         &      self%abundanceMinimum   &
         &     )
    if (self%abundanceMaximum > 0.0d0)  &
         & haloMassFunctionSample       &
         & =min(                        &
         &      haloMassFunctionSample, &
         &      self%abundanceMaximum   &
         &     )
    return
  end function haloMassFunctionSample
