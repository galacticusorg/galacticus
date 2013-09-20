!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements geometries of galaxy surveys.

module Geometry_Surveys
  !% Implements geometries of galaxy surveys.
  use ISO_Varying_String
  implicit none
  private
  public :: Geometry_Survey_Distance_Maximum, Geometry_Survey_Solid_Angle, Geometry_Survey_Volume_Maximum,&
       & Geometry_Survey_Window_Functions
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: moduleInitialized=.false.

  ! Name of survey geometry method used.
  type(varying_string) :: surveyGeometryMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Geometry_Survey_Distance_Maximum), pointer :: Geometry_Survey_Distance_Maximum_Get => null()
  procedure(Geometry_Survey_Solid_Angle     ), pointer :: Geometry_Survey_Solid_Angle_Get      => null()
  procedure(Geometry_Survey_Volume_Maximum  ), pointer :: Geometry_Survey_Volume_Maximum_Get   => null()
  procedure(Geometry_Survey_Window_Functions), pointer :: Geometry_Survey_Window_Functions_Get => null()

contains

  subroutine Geometry_Surveys_Initialize
    !% Initialize the spheroid star formation timecale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="surveyGeometryMethod" type="moduleUse">
    include 'geometry.surveys.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Geometry_Surveys_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the survey geometry method parameter.
          !@ <inputParameter>
          !@   <name>surveyGeometryMethod</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing survey geometries.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('surveyGeometryMethod',surveyGeometryMethod)
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="surveyGeometryMethod" type="functionCall" functionType="void">
          !#  <functionArgs>surveyGeometryMethod,Geometry_Survey_Distance_Maximum_Get,Geometry_Survey_Solid_Angle_Get,Geometry_Survey_Volume_Maximum_Get,Geometry_Survey_Window_Functions_Get</functionArgs>
          include 'geometry.surveys.inc'
          !# </include>
          if     (                                                                                        &
               &  .not.(                                                                                  &
               &             associated(Geometry_Survey_Distance_Maximum_Get)                             &
               &        .and.associated(Geometry_Survey_Volume_Maximum_Get  )                             &
               &        .and.associated(Geometry_Survey_Window_Functions_Get)                             &
               &       )                                                                                  &
               & ) call Galacticus_Error_Report(                                                          &
               &                                'Geometry_Surveys_Initialize',                            &
               &                                'method '//char(surveyGeometryMethod)//' is unrecognized' &
               &                               )
          moduleInitialized=.true.
       end if
       !$omp end critical(Geometry_Surveys_Initialization) 
    end if
    return
  end subroutine Geometry_Surveys_Initialize

  double precision function Geometry_Survey_Distance_Maximum(mass)
    !% Returns the maximum distance (in Mpc) at which a galaxy of the specified {\tt mass} (in $M_\odot$) could be detected.
    implicit none
    double precision, intent(in) :: mass

    ! Initialize the module.
    call Geometry_Surveys_Initialize()

    ! Get the energy using the selected method.
    Geometry_Survey_Distance_Maximum=Geometry_Survey_Distance_Maximum_Get(mass)
    
    return
  end function Geometry_Survey_Distance_Maximum
  
  double precision function Geometry_Survey_Solid_Angle()
    !% Returns the solid angle (in steradians) of the survey.
    implicit none

    ! Initialize the module.
    call Geometry_Surveys_Initialize()

    ! Get the energy using the selected method.
    Geometry_Survey_Solid_Angle=Geometry_Survey_Solid_Angle_Get()
    
    return
  end function Geometry_Survey_Solid_Angle
  
  double precision function Geometry_Survey_Volume_Maximum(mass)
    !% Returns the maximum volume (in Mpc$^3$) at which a galaxy of the specified {\tt mass} (in $M_\odot$) could be detected.
    implicit none
    double precision, intent(in) :: mass

    ! Initialize the module.
    call Geometry_Surveys_Initialize()

    ! Get the energy using the selected method.
    Geometry_Survey_Volume_Maximum=Geometry_Survey_Volume_Maximum_Get(mass)
    
    return
  end function Geometry_Survey_Volume_Maximum
  
  subroutine Geometry_Survey_Window_Functions(mass1,mass2,boxLength,gridCount,windowFunction1,windowFunction2)
    !% Returns the window functions on a grid of the specified size ({\tt gridCount} cells in each dimension) for galaxies of the
    !% specified {\tt mass1} and {\tt mass2} (in $M_\odot$). The {\tt boxLength} should be set to an appropriate value to fully
    !% enclose (with sufficient buffering to allow for Fourier transformation) the two window functions.
    use, intrinsic :: ISO_C_Binding
    implicit none
    double precision         , intent(in   )                                           :: mass1,mass2
    integer                  , intent(in   )                                           :: gridCount
    double precision         , intent(  out)                                           :: boxLength
    complex(c_double_complex), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2


    ! Initialize the module.
    call Geometry_Surveys_Initialize()

    ! Call the function that does the work.
    call Geometry_Survey_Window_Functions_Get(mass1,mass2,boxLength,gridCount,windowFunction1,windowFunction2)
    
    return
  end subroutine Geometry_Survey_Window_Functions
  
end module Geometry_Surveys
