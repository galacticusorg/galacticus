!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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
  
  !% Contains a module which implements a \cite{tinker_towardhalo_2008} dark matter halo mass function class.
  use Cosmological_Mass_Variance
  use Linear_Growth
  use Cosmology_Functions
  use Virial_Density_Contrast
  use Tables

  !# <enumeration>
  !#  <name>tinker2008Parameter</name>
  !#  <description>Enumeration of parameters for the {\normalfont \ttfamily tinker2008} halo mass function class.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <entry label="a"            />
  !#  <entry label="b"            />
  !#  <entry label="c"            />
  !#  <entry label="normalization"/>
  !# </enumeration>

  !# <haloMassFunction name="haloMassFunctionTinker2008">
  !#  <description>The halo mass function is computed from the function given by \cite{tinker_towardhalo_2008}.</description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionTinker2008
     !% A halo mass function class using the fitting function of \cite{tinker_towardhalo_2008}.
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
     class           (linearGrowthClass            ), pointer :: linearGrowth_
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_
     type            (table1DGeneric               )          :: densityContrast
     double precision                                         :: time                     , mass, &
          &                                                      massFunction
   contains
     final     ::                 tinker2008Destructor
     procedure :: differential => tinker2008Differential
  end type haloMassFunctionTinker2008

  interface haloMassFunctionTinker2008
     !% Constructors for the {\normalfont \ttfamily tinker2008} halo mass function class.
     module procedure tinker2008ConstructorParameters
     module procedure tinker2008ConstructorInternal
  end interface haloMassFunctionTinker2008

contains

  function tinker2008ConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily tinker2008} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionTinker2008    )                :: tinker2008ConstructorParameters
    type(inputParameters               ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_    
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass            ), pointer       :: linearGrowth_
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    tinker2008ConstructorParameters=tinker2008ConstructorInternal(                            &
         &                                                         cosmologyParameters_     , &
         &                                                         cosmologicalMassVariance_, &
         &                                                         linearGrowth_            , &
         &                                                         cosmologyFunctions_        &
         &                                                        )
    return
  end function tinker2008ConstructorParameters

  function tinker2008ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_)
    !% Internal constructor for the {\normalfont \ttfamily tinker2008} halo mass function class.
    use FoX_DOM
    use IO_XML
    use Galacticus_Error
    use Galacticus_Input_Paths
    use Table_Labels
    implicit none
    type            (haloMassFunctionTinker2008   )                             :: tinker2008ConstructorInternal
    class           (cosmologyParametersClass     ), target     , intent(in   ) :: cosmologyParameters_    
    class           (cosmologicalMassVarianceClass), target     , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), target     , intent(in   ) :: linearGrowth_
    class           (cosmologyFunctionsClass      ), target     , intent(in   ) :: cosmologyFunctions_
    type            (node                         ), pointer                    :: doc                          , columnsElement, &
         &                                                                         columnElement
    double precision                               , allocatable, dimension(:)  :: dataTmp
    integer                                                                     :: i                            , ioStatus
    
    tinker2008ConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    tinker2008ConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    tinker2008ConstructorInternal%linearGrowth_             => linearGrowth_
    tinker2008ConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    tinker2008ConstructorInternal%time                      =  -1.0d0
    tinker2008ConstructorInternal%mass                      =  -1.0d0
    ! Read the data file which gives fitting parameters as a function of halo overdensity.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(Galacticus_Input_Path())//"data/darkMatter/Halo_Mass_Function_Parameters_Tinker_2008.xml",ioStat=ioStatus)
    if (ioStatus /= 0) call Galacticus_Error_Report('tinker2008ConstructorInternal','Unable to find data file')
    columnsElement => XML_Get_First_Element_By_Tag_Name(doc           ,"columns"        )
    columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"densityContrast")
    call XML_Array_Read(columnElement,"data",dataTmp)
    call tinker2008ConstructorInternal%densityContrast%create(dataTmp,tinker2008ParameterCount,[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
    deallocate(dataTmp)
    do i=tinker2008ParameterMin,tinker2008ParameterMax
       columnElement => XML_Get_First_Element_By_Tag_Name(columnsElement,char(enumerationTinker2008ParameterDecode(i,includePrefix=.false.)))
       call XML_Array_Read(columnElement,"data",dataTmp)
       call tinker2008ConstructorInternal%densityContrast%populate(dataTmp,table=i+1)
       deallocate(dataTmp)
    end do
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    return
  end function tinker2008ConstructorInternal
  
  subroutine tinker2008Destructor(self)
    !% Destructor for the {\normalfont \ttfamily tinker2008} halo mass function class.
    implicit none
    type(haloMassFunctionTinker2008), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"       />
    !# <objectDestructor name="self%cosmologicalMassVariance_" />
    !# <objectDestructor name="self%linearGrowth_"             />
    !# <objectDestructor name="self%cosmologyParameters_"      />
    return
  end subroutine tinker2008Destructor

  double precision function tinker2008Differential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    use FGSL
    use Numerical_Interpolation
    use Table_Labels
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout)                                            :: self
    double precision                            , intent(in   )                                            :: time                  , mass
    double precision                            , dimension(tinker2008ParameterMin:tinker2008ParameterMax) :: parameters
    class           (virialDensityContrastClass), pointer                                                  :: virialDensityContrast_
    integer                                                                                                :: i
    double precision                                                                                       :: alpha                 , sigma, &
         &                                                                                                    densityContrast       , a    , &
         &                                                                                                    alphaDensityContrast  , b    , &
         &                                                                                                    expansionFactor       , c    , &
         &                                                                                                    normalization 

    ! Update fitting function parameters if the time differs from that on the previous call.
    if (time /= self%time .or. mass /= self%mass) then
       ! Get halo virial density contrast, and expansion factor.
       virialDensityContrast_ =>      virialDensityContrast                 (         )
       expansionFactor        =  self%cosmologyFunctions_   %expansionFactor(     time)
       densityContrast        =       virialDensityContrast_%densityContrast(mass,time)
       ! Compute coefficients of fitting function.
       do i=tinker2008ParameterMin,tinker2008ParameterMax
          parameters(i)=self%densityContrast%interpolate(densityContrast,table=i+1)
       end do
       ! Extrapolate to higher redshift using redshift scalings given by Tinker et al. (2008; eqns. 5-8).
       alphaDensityContrast=+10.0d0**(                          &
            &                         -(                        &
            &                           +0.75d0                 &
            &                           /log10(                 &
            &                                  +densityContrast &
            &                                  /75.0d0          &
            &                                 )                 &
            &                          )**1.2d0                 &
            &                        )
       normalization       =parameters(tinker2008ParameterNormalization)*expansionFactor**0.14d0
       a                   =parameters(tinker2008ParameterA            )*expansionFactor**0.06d0
       b                   =parameters(tinker2008ParameterB            )*expansionFactor**alphaDensityContrast
       c                   =parameters(tinker2008ParameterC            )
       ! Compute and store the mass function.
       sigma            =+    self%cosmologicalMassVariance_%rootVariance                   (mass)  &
            &            *    self%linearGrowth_            %value                          (time)
       alpha            =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass))
       self%massFunction=+    self%cosmologyParameters_     %OmegaMatter                    (    )  &
            &            *    self%cosmologyParameters_     %densityCritical                (    )  &
            &            /mass**2                                                                   &
            &            *alpha                                                                     &
            &            *normalization                                                             &
            &            *exp(                                                                      &
            &                 -c                                                                    &
            &                  /sigma**2                                                            &
            &                )                                                                      &
            &            *(                                                                         &
            &              +1.0d0                                                                   &
            &              +(                                                                       &
            &                +b                                                                     &
            &                /sigma                                                                 &
            &               )**a                                                                    &
            &             )
       ! Store the time and mass.
       self%time=time
       self%mass=mass
    end if
    ! Return the stored mass function.
    tinker2008Differential=self%massFunction
    return
  end function tinker2008Differential
