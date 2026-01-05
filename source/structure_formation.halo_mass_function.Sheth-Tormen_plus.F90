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
  Contains a class which implements an extended version of the \cite{sheth_ellipsoidal_2001} dark matter halo mass function class.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass
  
  !![
  <haloMassFunction name="haloMassFunctionShethTormenPlus">
   <description>The halo mass function is computed using an extended version of the \cite{sheth_ellipsoidal_2001} fitting function.</description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionShethTormenPlus
     !!{
     A halo mass function class using an extended version of the \cite{sheth_ellipsoidal_2001} fitting function.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (linearGrowthClass            ), pointer :: linearGrowth_             => null()
     double precision                                         :: aValue                             , pValue, &
          &                                                      normalizationValue                 , qValue, &
          &                                                      bValue                             , cValue, &
          &                                                      dValue
   contains
     !![
     <methods>
       <method description="Return the parameter $\bar{a}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="a"            />
       <method description="Return the parameter $\bar{b}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="b"            />
       <method description="Return the parameter $\bar{c}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="c"            />
       <method description="Return the parameter $\bar{d}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="d"            />
       <method description="Return the parameter $\bar{p}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="p"            />
       <method description="Return the parameter $\bar{q}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="q"            />
       <method description="Return the parameter $\bar{A}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit." method="normalization"/>
     </methods>
     !!]
     final     ::                  shethTormenPlusDestructor
     procedure :: differential  => shethTormenPlusDifferential
     procedure :: a             => shethTormenPlusA
     procedure :: b             => shethTormenPlusB
     procedure :: c             => shethTormenPlusC
     procedure :: d             => shethTormenPlusD
     procedure :: p             => shethTormenPlusP
     procedure :: q             => shethTormenPlusQ
     procedure :: normalization => shethTormenPlusNormalization
     procedure :: descriptor    => shethTormenPlusDescriptor
  end type haloMassFunctionShethTormenPlus

  interface haloMassFunctionShethTormenPlus
     !!{
     Constructors for the \refClass{haloMassFunctionShethTormenPlus} halo mass function class.
     !!}
     module procedure shethTormenPlusConstructorParameters
     module procedure shethTormenPlusConstructorInternal
  end interface haloMassFunctionShethTormenPlus

contains

  function shethTormenPlusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionShethTormenPlus} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionShethTormenPlus)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class           (linearGrowthClass              ), pointer       :: linearGrowth_
    double precision                                                 :: a                        , p, &
         &                                                              normalization            , q, &
         &                                                              b                        , c, &
         &                                                              d

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.788d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{a}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>b</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\bar{b}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>c</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\bar{c}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>d</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\bar{d}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>p</name>
      <source>parameters</source>
      <defaultValue>0.807d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{p}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>q</name>
      <source>parameters</source>
      <defaultValue>1.795d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{q}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.333d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The normalization parameter $\bar{A}$ in the extended \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    !!]
    self=haloMassFunctionShethTormenPlus(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_,a,b,c,d,p,q,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function shethTormenPlusConstructorParameters

  function shethTormenPlusConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_,a,b,c,d,p,q,normalization) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionShethTormenPlus} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionShethTormenPlus)                        :: self
    class           (cosmologyParametersClass       ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass  ), target, intent(in   ) :: cosmologicalMassVariance_
    class           (criticalOverdensityClass       ), target, intent(in   ) :: criticalOverdensity_
    class           (linearGrowthClass              ), target, intent(in   ) :: linearGrowth_
    double precision                                         , intent(in   ) :: a                        , p, &
         &                                                                      normalization            , q, &
         &                                                                      b                        , c, &
         &                                                                      d
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *criticalOverdensity_, *linearGrowth_"/>
    !!]

    self%            aValue=a
    self%            bValue=b
    self%            cValue=c
    self%            dValue=d
    self%            pValue=p
    self%            qValue=q
    self%normalizationValue=normalization
    return
  end function shethTormenPlusConstructorInternal

  subroutine shethTormenPlusDestructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionShethTormenPlus} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionShethTormenPlus), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine shethTormenPlusDestructor

  double precision function shethTormenPlusDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: time      , mass
    type            (treeNode                       ), intent(inout), optional :: node
    double precision                                                           :: alpha     , nu          , &
         &                                                                        nuPrime   , massVariance, &
         &                                                                        growthRate

    ! Determine the mass variance. If zero, return zero mass function.
    massVariance=self%cosmologicalMassVariance_%rootVariance(mass,time)
    if (massVariance <= 0.0d0) then
       shethTormenPlusDifferential=0.0d0
       return
    end if
    ! Compute the mass function.
    nu                     =+(                                                                &
         &                    +self%criticalOverdensity_%value(time=time,mass=mass,node=node) &
         &                    /massVariance                                                   &
         &                   )**2
    if (nu <= 0.0d0) then
       shethTormenPlusDifferential=0.0d0
       return
    end if
    nuPrime                =+self%a(time,mass)                                                &
         &                  *nu
    alpha                  =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time))
    growthRate             =self%linearGrowth_%logarithmicDerivativeExpansionFactor(time=time)
    shethTormenPlusDifferential=+self%cosmologyParameters_%OmegaMatter    () &
         &                      *self%cosmologyParameters_%densityCritical() &
         &                      /mass**2                                     &
         &                      *alpha                                       &
         &                      *sqrt(                                       &
         &                            +2.0d0                                 &
         &                            *nuPrime**self%q(time,mass)            &
         &                            /Pi                                    &
         &                           )                                       &
         &                      *self%normalization(time,mass)               &
         &                      *(                                           &
         &                        +1.0d0                                     &
         &                        +1.0d0                                     &
         &                        /nuPrime**self%p(time,mass)                &
         &                       )                                           &
         &                      *exp(                                        &
         &                           -0.5d0                                  &
         &                           *nuPrime                                &
         &                           -self%b(time,mass)                      &
         &                           *growthRate**self%c(time,mass)          &
         &                           *nuPrime   **self%d(time,mass)          &
         &                          )
    return
  end function shethTormenPlusDifferential

  double precision function shethTormenPlusA(self,time,mass)
    !!{
    Return the parameter $\bar{a}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusA=self%aValue
    return
  end function shethTormenPlusA

  double precision function shethTormenPlusB(self,time,mass)
    !!{
    Return the parameter $\bar{b}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusB=self%bValue
    return
  end function shethTormenPlusB

  double precision function shethTormenPlusC(self,time,mass)
    !!{
    Return the parameter $\bar{c}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusC=self%cValue
    return
  end function shethTormenPlusC

  double precision function shethTormenPlusD(self,time,mass)
    !!{
    Return the parameter $\bar{d}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusD=self%dValue
    return
  end function shethTormenPlusD

  double precision function shethTormenPlusP(self,time,mass)
    !!{
    Return the parameter $\bar{p}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusP=self%pValue
    return
  end function shethTormenPlusP

  double precision function shethTormenPlusQ(self,time,mass)
    !!{
    Return the parameter $\bar{q}$ in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusQ=self%qValue
    return
  end function shethTormenPlusQ

  double precision function shethTormenPlusNormalization(self,time,mass)
    !!{
    Return the normalization, $\bar{A}$, in the {\normalfont \ttfamily shethTormenPlus} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormenPlus), intent(inout) :: self
    double precision                                 , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    shethTormenPlusNormalization=self%normalizationValue
    return
  end function shethTormenPlusNormalization

  subroutine shethTormenPlusDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (haloMassFunctionShethTormenPlus), intent(inout)           :: self
    type     (inputParameters                ), intent(inout)           :: descriptor
    logical                                   , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                         )                          :: parameterLabel
    type     (inputParameters                )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('haloMassFunction','shethTormenPlus')
    parameters=descriptor%subparameters('haloMassFunction')
    write (parameterLabel,'(e17.10)') self%aValue
    call parameters%addParameter('a'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%bValue
    call parameters%addParameter('b'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%cValue
    call parameters%addParameter('c'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%dValue
    call parameters%addParameter('d'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%pValue
    call parameters%addParameter('p'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%qValue
    call parameters%addParameter('q'            ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%normalizationValue
    call parameters%addParameter('normalization',trim(adjustl(parameterLabel)))
    call self%cosmologyParameters_     %descriptor(parameters,includeClass,includeFileModificationTimes)
    call self%cosmologicalMassVariance_%descriptor(parameters,includeClass,includeFileModificationTimes)
    call self%criticalOverdensity_     %descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine shethTormenPlusDescriptor
