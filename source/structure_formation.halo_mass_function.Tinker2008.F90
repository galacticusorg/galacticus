!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implements a \cite{tinker_towardhalo_2008} dark matter halo mass function class.
  !!}
  use :: Tables                 , only : table1DGeneric
  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <enumeration>
   <name>tinker2008Parameter</name>
   <description>Enumeration of parameters for the {\normalfont \ttfamily tinker2008} halo mass function class.</description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <entry label="a"            />
   <entry label="b"            />
   <entry label="c"            />
   <entry label="normalization"/>
  </enumeration>
  !!]

  !![
  <haloMassFunction name="haloMassFunctionTinker2008">
   <description>
    A dark matter halo mass function class using the function given by \cite{tinker_towardhalo_2008}, and using their fits for
    the parameter values at the appropriate virial density contrast.
   </description>
   <stateStorable>
    <exclude variables="densityContrast"/>
   </stateStorable>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionTinker2008Form) :: haloMassFunctionTinker2008
     !!{
     A halo mass function class using the fitting function of \cite{tinker_towardhalo_2008}, and using their fits for the parameter values.
     !!}
     private
     class           (virialDensityContrastClass), pointer                                                  :: virialDensityContrast_ => null()
     type            (table1DGeneric            )                                                           :: densityContrast
     double precision                            , dimension(tinker2008ParameterMin:tinker2008ParameterMax) :: parameters
     double precision                                                                                       :: alphaDensityContrast            , timeParameters, &
          &                                                                                                    massParameters
   contains
     !![
     <methods>
       <method description="Evaluate hyper-parameters needed for the fitting function." method="parametersEvaluate" />
     </methods>
     !!]
     final     ::                       tinker2008Destructor
     procedure :: normalization      => tinker2008Normalization
     procedure :: a                  => tinker2008A
     procedure :: b                  => tinker2008B
     procedure :: c                  => tinker2008C
     procedure :: parametersEvaluate => tinker2008ParametersEvaluate
  end type haloMassFunctionTinker2008

  interface haloMassFunctionTinker2008
     !!{
     Constructors for the \refClass{haloMassFunctionTinker2008} halo mass function class.
     !!}
     module procedure tinker2008ConstructorParameters
     module procedure tinker2008ConstructorInternal
  end interface haloMassFunctionTinker2008

contains

  function tinker2008ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionTinker2008} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionTinker2008   )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass            ), pointer       :: linearGrowth_
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class(virialDensityContrastClass   ), pointer       :: virialDensityContrast_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self=haloMassFunctionTinker2008(                           &
         &                          cosmologyParameters_     , &
         &                          cosmologicalMassVariance_, &
         &                          linearGrowth_            , &
         &                          cosmologyFunctions_      , &
         &                          virialDensityContrast_     &
         &                         )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function tinker2008ConstructorParameters

  function tinker2008ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionTinker2008} halo mass function class.
    !!}
    use :: File_Utilities    , only : File_Exists
    use :: FoX_DOM           , only : destroy                     , node
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath                   , pathTypeDataStatic
    use :: IO_XML            , only : XML_Array_Read              , XML_Get_First_Element_By_Tag_Name, XML_Parse
    use :: ISO_Varying_String, only : varying_string
    use :: Table_Labels      , only : extrapolationTypeExtrapolate
    implicit none
    type            (haloMassFunctionTinker2008   )                             :: self
    class           (cosmologyParametersClass     ), target     , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), target     , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), target     , intent(in   ) :: linearGrowth_
    class           (cosmologyFunctionsClass      ), target     , intent(in   ) :: cosmologyFunctions_
    class           (virialDensityContrastClass   ), target     , intent(in   ) :: virialDensityContrast_
    type            (node                         ), pointer                    :: doc                          , columnsElement, &
         &                                                                         columnElement
    double precision                               , allocatable, dimension(:)  :: dataTmp
    integer                                                                     :: i                            , ioStatus
    type            (varying_string               )                             :: parameterFileName
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *linearGrowth_, *cosmologyFunctions_, *virialDensityContrast_"/>
    !!]

    self%time          =-1.0d0
    self%mass          =-1.0d0
    self%timeParameters=-1.0d0
    self%massParameters=-1.0d0
    ! Read the data file which gives fitting parameters as a function of halo overdensity.
    parameterFileName=inputPath(pathTypeDataStatic)//"darkMatter/Halo_Mass_Function_Parameters_Tinker_2008.xml"
    if (.not.File_Exists(parameterFileName)) call Error_Report('Unable to find data file "'//parameterFileName//'"'//{introspection:location})
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(parameterFileName,ioStat=ioStatus)
    if (ioStatus /= 0) call Error_Report('Unable to parse data file "'//parameterFileName//'"'//{introspection:location})
    columnsElement => XML_Get_First_Element_By_Tag_Name(doc           ,"columns"        )
    columnElement  => XML_Get_First_Element_By_Tag_Name(columnsElement,"densityContrast")
    call XML_Array_Read(columnElement,"data",dataTmp)
    call self%densityContrast%create(dataTmp,tinker2008ParameterCount,[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
    deallocate(dataTmp)
    do i=tinker2008ParameterMin,tinker2008ParameterMax       
       columnElement => XML_Get_First_Element_By_Tag_Name(columnsElement,char(enumerationTinker2008ParameterDecode(i,includePrefix=.false.)))
       call XML_Array_Read(columnElement,"data",dataTmp)
       call self%densityContrast%populate(dataTmp,table=i+1)
       deallocate(dataTmp)
    end do
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    return
  end function tinker2008ConstructorInternal

  subroutine tinker2008Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionTinker2008} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionTinker2008), intent(inout) :: self

    call self%densityContrast%destroy()
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%virialDensityContrast_"   />
    !!]
    return
  end subroutine tinker2008Destructor

  subroutine tinker2008ParametersEvaluate(self,time,mass)
    !!{
    Evaluate interpolating parameters for the {\normalfont \ttfamily tinker2008} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout) :: self
    double precision                            , intent(in   ) :: time           , mass
    integer                                                     :: i
    double precision                                            :: expansionFactor, densityContrast

    if     (                             &
         &   time /= self%timeParameters &
         &  .or.                         &
         &   mass /= self%massParameters &
         & ) then
       ! Get halo virial density contrast, and expansion factor.
       expansionFactor=self%cosmologyFunctions_   %expansionFactor(     time)
       densityContrast=self%virialDensityContrast_%densityContrast(mass,time)
       ! Compute coefficients of fitting function.
       do i=tinker2008ParameterMin,tinker2008ParameterMax
          self%parameters(i)=self%densityContrast%interpolate(densityContrast,table=i+1)
       end do
       ! Extrapolate to higher redshift using redshift scalings given by Tinker et al. (2008; eqns. 5-8).
       self%alphaDensityContrast=+10.0d0**(                          &
            &                              -(                        &
            &                                +0.75d0                 &
            &                                /log10(                 &
            &                                       +densityContrast &
            &                                       /75.0d0          &
            &                                      )                 &
            &                               )**1.2d0                 &
            &                             )
    end if
    return
  end subroutine tinker2008ParametersEvaluate

  double precision function tinker2008Normalization(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily tinker2008} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout) :: self
    double precision                            , intent(in   ) :: time           , mass
    double precision                                            :: expansionFactor

    call self%parametersEvaluate(time,mass)
    expansionFactor        =self%cosmologyFunctions_%expansionFactor(                            time   )
    tinker2008Normalization=self%parameters                         (tinker2008ParameterNormalization%ID)*expansionFactor**0.14d0
    return
  end function tinker2008Normalization

  double precision function tinker2008A(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily tinker2008} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout) :: self
    double precision                            , intent(in   ) :: time           , mass
    double precision                                            :: expansionFactor

    call self%parametersEvaluate(time,mass)
    expansionFactor        =self%cosmologyFunctions_%expansionFactor(                            time   )
    tinker2008A            =self%parameters                         (tinker2008ParameterA            %ID)*expansionFactor**0.06d0
    return
  end function tinker2008A

  double precision function tinker2008B(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily tinker2008} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout) :: self
    double precision                            , intent(in   ) :: time           , mass
    double precision                                            :: expansionFactor

    call self%parametersEvaluate(time,mass)
    expansionFactor        =self%cosmologyFunctions_%expansionFactor(                            time   )
    tinker2008B            =self%parameters                         (tinker2008ParameterB            %ID)*expansionFactor**self%alphaDensityContrast
    return
  end function tinker2008B

  double precision function tinker2008C(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily tinker2008} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008), intent(inout) :: self
    double precision                            , intent(in   ) :: time           , mass
    double precision                                            :: expansionFactor

    call self%parametersEvaluate(time,mass)
    expansionFactor        =self%cosmologyFunctions_%expansionFactor(                            time   )
    tinker2008C            =self%parameters                         (tinker2008ParameterC            %ID)
    return
  end function tinker2008C

