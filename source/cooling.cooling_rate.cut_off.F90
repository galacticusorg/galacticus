!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a cooling rate class which modifies another cooling rate by cutting off cooling above some virial velocity.

  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !# <coolingRate name="coolingRateCutOff">
  !#  <description>A cooling rate class which modifies another cooling rate by cutting off cooling above some virial velocity.</description>
  !# </coolingRate>
  type, extends(coolingRateClass) :: coolingRateCutOff
     !% Implementation of cooling rate class which modifies another cooling rate by cutting off cooling above some virial velocity.
     private
     class           (coolingRateClass        ), pointer :: coolingRate_         => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     ! Parameters controlling the cut off.
     double precision                            :: velocityCutOff  , timeCutOff
     integer                                     :: whenCutOff
     logical                                     :: useFormationNode
   contains
     final     ::         cutOffDestructor
     procedure :: rate => cutOffRate
  end type coolingRateCutOff

  interface coolingRateCutOff
     !% Constructors for the cut off cooling rate class.
     module procedure cutOffConstructorParameters
     module procedure cutOffConstructorInternal
  end interface coolingRateCutOff

  ! Enumeration for whether cut off is before or after the given epoch.
  !# <enumeration>
  !#  <name>cutOffWhen</name>
  !#  <description>Specifies whether cooling is cut off before or after the given epoch.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <entry label="before" />
  !#  <entry label="after"  />
  !# </enumeration>

contains

  function cutOffConstructorParameters(parameters) result(self)
    !% Constructor for the cut off cooling rate class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateCutOff       )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (coolingRateClass        ), pointer       :: coolingRate_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    double precision                                          :: velocityCutOff      , redshiftCutOff
    logical                                                   :: useFormationNode
    type            (varying_string          )                :: whenCutOff

    !# <inputParameter>
    !#   <name>useFormationNode</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Specifies whether to use the virial velocity of the formation node or current node in the cooling rate ``cut-off'' modifier.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>velocityCutOff</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The velocity below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftCutOff</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The redshift below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>whenCutOff</name>
    !#   <defaultValue>var_str('after')</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Specifies whether cooling is cut off before or after {\normalfont \ttfamily [redshiftCutOff]}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="coolingRate"         name="coolingRate_"         source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=coolingRateCutOff(                                                                        &
         &                                                                   velocityCutOff      , &
         &                 cosmologyFunctions_ %cosmicTime                 (                       &
         &                  cosmologyFunctions_%expansionFactorFromRedshift (                      &
         &                                                                   redshiftCutOff        &
         &                                                                  )                      &
         &                                                                 )                     , &
         &                 enumerationCutOffWhenEncode                     (                       &
         &                                                                   char(whenCutOff)      &
         &                                                                 )                     , &
         &                                                                   useFormationNode    , &
         &                                                                   darkMatterHaloScale_, &
         &                                                                   coolingRate_          &
         &                )
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    !# <objectDestructor name="coolingRate_"        />
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function cutOffConstructorParameters

  function cutOffConstructorInternal(velocityCutOff,timeCutOff,whenCutOff,useFormationNode,darkMatterHaloScale_,coolingRate_) result(self)
    !% Internal constructor for the cut off cooling rate class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (coolingRateCutOff       )                        :: self
    double precision                          , intent(in   )         :: velocityCutOff      , timeCutOff
    integer                                   , intent(in   )         :: whenCutOff
    logical                                   , intent(in   )         :: useFormationNode
    class           (darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class           (coolingRateClass        ), intent(in   ), target :: coolingRate_

    !# <constructorAssign variables="velocityCutOff, timeCutOff, whenCutOff, useFormationNode, *coolingRate_, *darkMatterHaloScale_"/>
    ! Validate "whenCutOff" argument.
    if (.not.enumerationCutOffWhenIsValid(whenCutOff)) call Galacticus_Error_Report('[whenCutOff] is invalid'//{introspection:location})
    return
  end function cutOffConstructorInternal

  subroutine cutOffDestructor(self)
    !% Destructor for the cut off cooling rate class.
    implicit none
    type(coolingRateCutOff), intent(inout) :: self

    !# <objectDestructor name="self%coolingRate_"        />
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine cutOffDestructor

  double precision function cutOffRate(self,node)
    !% Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is cut off
    !% before/after a given epoch and below a given virial velocity.
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (coolingRateCutOff ), intent(inout) :: self
    type            (treeNode          ), intent(inout) :: node
    class           (nodeComponentBasic), pointer       :: basic
    double precision                                    :: velocityVirial

    ! Test for halos where cooling should be cut off.
    select case (self%useFormationNode)
    case (.false.)
       velocityVirial=self%darkMatterHaloScale_%virialVelocity(node              )
    case (.true. )
       velocityVirial=self%darkMatterHaloScale_%virialVelocity(node%formationNode)
    end select
    basic => node%basic()
    if     (                                                                             &
         &  (                                                                            &
         &   (basic%time() >= self%timeCutOff .and. self%whenCutOff == cutOffWhenAfter ) &
         &    .or.                                                                       &
         &   (basic%time() <= self%timeCutOff .and. self%whenCutOff == cutOffWhenBefore) &
         &  )                                                                            &
         &   .and.                                                                       &
         &  velocityVirial            <= self%velocityCutOff                             &
         & ) then
       cutOffRate=0.0d0
    else
       cutOffRate=self%coolingRate_%rate(node)
    end if
    return
  end function cutOffRate

