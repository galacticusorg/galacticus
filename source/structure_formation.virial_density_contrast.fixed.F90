!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of fixed dark matter halo virial density contrasts.

  !# <virialDensityContrast name="virialDensityContrastFixed">
  !#  <description>Fixed dark matter halo virial density contrasts.</description>
  !# </virialDensityContrast>

  type, extends(virialDensityContrastClass) :: virialDensityContrastFixed
     !% A dark matter halo virial density contrast class assuming fixed contrast.
     private
     double precision :: densityContrastValue
     integer          :: densityType
   contains
     procedure :: densityContrast             => fixedDensityContrast
     procedure :: densityContrastRateOfChange => fixedDensityContrastRateOfChange
  end type virialDensityContrastFixed

  interface virialDensityContrastFixed
     !% Constructors for the {\tt fixed} dark matter halo virial density contrast class.
     module procedure fixedDefaultConstructor
     module procedure fixedConstructor
  end interface virialDensityContrastFixed

  ! Initialization status.
  logical                                             :: fixedInitialized                             =.false.

  ! The type of reference density to use.
  integer                         , parameter, public :: virialDensityContrastFixedDensityTypeCritical=0
  integer                         , parameter, public :: virialDensityContrastFixedDensityTypeMean    =1

  ! Parameters for the default implementation.
  double precision                                    :: virialDensityContrastFixedValue
  type            (varying_string)                    :: virialDensityContrastFixedTypeText
  integer                                             :: virialDensityContrastFixedType
  
contains

  function fixedDefaultConstructor()
    !% Default constructor for the {\tt fixed} dark matter halo virial density contrast class.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(virialDensityContrastFixed), target :: fixedDefaultConstructor

    if (.not.fixedInitialized) then
       !$omp critical(fixedDefaultInitialize)
       if (.not.fixedInitialized) then
          !@ <inputParameter>
          !@   <name>virialDensityContrastFixed</name>
          !@   <defaultValue>200</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The virial density contrast to use in the fixed value model.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("virialDensityContrastFixed"    ,virialDensityContrastFixedValue   ,defaultValue=200.0d0           )
          !@ <inputParameter>
          !@   <name>virialDensityContrastFixedType</name>
          !@   <defaultValue>criticalDensity</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The reference density to use in the fixed value virial density contrast model. Either of {\tt critical density} and {\tt mean density} are allowed.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("virialDensityContrastFixedType",virialDensityContrastFixedTypeText,defaultValue='criticalDensity')
          select case (char(virialDensityContrastFixedTypeText))
          case ("criticalDensity")
             virialDensityContrastFixedType=virialDensityContrastFixedDensityTypeCritical
          case ("meanDensity"    )
             virialDensityContrastFixedType=virialDensityContrastFixedDensityTypeMean
          case default
             call Galacticus_Error_Report('fixedDefaultConstructor','[virialDensityContrastFixedType] must be either "criticalDensity" or "meanDensity"')
          end select
          fixedInitialized=.true.
       end if
       !$omp end critical(fixedDefaultInitialize)
    end if
    ! Construct the object.
    fixedDefaultConstructor=fixedConstructor(virialDensityContrastFixedValue,virialDensityContrastFixedType)
    return
  end function fixedDefaultConstructor

  function fixedConstructor(densityContrast,densityType)
    !% Constructor for the {\tt fixed} dark matter halo virial density contrast class.
    use Galacticus_Error
    implicit none
    type            (virialDensityContrastFixed)                :: fixedConstructor
    double precision                            , intent(in   ) :: densityContrast
    integer                                     , intent(in   ) :: densityType

    if     (                                                              &
         &   densityType /= virialDensityContrastFixedDensityTypeCritical &
         &  .and.                                                         &
         &   densityType /= virialDensityContrastFixedDensityTypeMean     &
         & ) call Galacticus_Error_Report('fixedConstructor','invalid densityType')
    fixedConstructor%densityContrastValue=densityContrast
    fixedConstructor%densityType         =densityType
    return
  end function fixedConstructor
  
  double precision function fixedDensityContrast(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming a fixed contrast.
    use Cosmology_Functions
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   ), optional :: time               , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass   ), pointer                 :: cosmologyFunctions_

    ! Set the density contrast.
    fixedDensityContrast=self%densityContrastValue
    ! If density contrast is specified relative to critical density, convert to mean density.
    if (self%densityType == virialDensityContrastFixedDensityTypeCritical) then
       cosmologyFunctions_ => cosmologyFunctions()
       fixedDensityContrast                                              &
            & =fixedDensityContrast                                      &
            & /cosmologyFunctions_ %omegaMatterEpochal(                  &
            &   cosmologyFunctions_%epochTime          (                 &
            &                                           time           , &
            &                                           expansionFactor, &
            &                                           collapsing       &
            &                                          )                 &
            &                                         )
    end if
    return
  end function fixedDensityContrast

  double precision function fixedDensityContrastRateOfChange(self,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, assuming a fixed contrast.
    use Cosmology_Functions
    implicit none
    class           (virialDensityContrastFixed), intent(inout)           :: self
    double precision                            , intent(in   ), optional :: time      , expansionFactor
    logical                                     , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass   ), pointer                 :: cosmologyFunctions_
    double precision                                                      :: epochTime

    ! Zero rate of change for fixed density contrast.
    fixedDensityContrastRateOfChange=0.0d0
    ! If density contrast is defined relative to critical density, include the rate of change of
    ! critical mean density.
    if (self%densityType == virialDensityContrastFixedDensityTypeCritical) then
       cosmologyFunctions_ => cosmologyFunctions()
       epochTime=cosmologyFunctions_%epochTime(                 &
            &                                  time           , &
            &                                  expansionFactor, &
            &                                  collapsing       &
            &                                 )
       fixedDensityContrastRateOfChange                                    &
            & =-self%densityContrastValue                                  &
            &  *cosmologyFunctions_ %omegaMatterRateOfChange(epochTime)    &
            &  /cosmologyFunctions_ %omegaMatterEpochal     (epochTime)**2
    end if
    return
  end function fixedDensityContrastRateOfChange
