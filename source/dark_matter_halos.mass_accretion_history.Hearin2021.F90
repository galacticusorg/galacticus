!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of dark matter halo mass accretion histories using the rolling power-law model of \cite{hearin_differentiable_2021}.
  !!}

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryHearin2021">
   <description>Dark matter halo mass accretion histories using the rolling power-law model of \cite{hearin_differentiable_2021}.</description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryHearin2021
     !!{
     A dark matter halo mass accretion history class using the rolling power-law model of \cite{hearin_differentiable_2021}.
     !!}
     private
     double precision :: powerLawIndexEarly, powerLawIndexLate, &
          &              rateRollOver      , timeMaximum
   contains
     !![
     <methods>
      <method description="Return the power law index at the given time."                      method="powerLawIndex"          />
      <method description="Return the derivative of the power law index with respect to time." method="powerLawIndexDerivative"/>
      <method description="Return the $\log_{10}(t_0)$ parameter."                             method="timeZeroLogarithmic"    />
      <method description="Return the maximum mass in the mass accretion history."             method="massMaximum"            />
      <method description="The sigmoid interpolation function."                                method="sigmoid"                />
      <method description="Return the early-time power law index."                             method="powerLawIndexEarly_"    />
      <method description="Return the late-time power law index."                              method="powerLawIndexLate_"     />
      <method description="Return the roll-over rate."                                         method="rateRollOver_"          />
      <method description="Return the time of maximum mass."                                   method="timeMaximum_"           />
     </methods>
     !!]
     procedure         :: mass                    => hearin2021Mass
     procedure         :: massAccretionRate       => hearin2021MassAccretionRate
     procedure         :: powerLawIndex           => hearin2021PowerLawIndex
     procedure         :: powerLawIndexDerivative => hearin2021PowerLawIndex
     procedure         :: massMaximum             => hearin2021MassMaximum
     procedure         :: timeZeroLogarithmic     => hearin2021TimeZeroLogarithmic
     procedure, nopass :: sigmoid                 => hearin2021Sigmoid
     procedure         :: powerLawIndexEarly_     => hearin2021PowerLawIndexEarly
     procedure         :: powerLawIndexLate_      => hearin2021PowerLawIndexLate
     procedure         :: rateRollOver_           => hearin2021RateRollOver
     procedure         :: timeMaximum_            => hearin2021TimeMaximum
  end type darkMatterHaloMassAccretionHistoryHearin2021

  interface darkMatterHaloMassAccretionHistoryHearin2021
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryHearin2021} dark matter halo mass accretion history class.
     !!}
     module procedure hearin2021ConstructorParameters
     module procedure hearin2021ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryHearin2021

contains

  function hearin2021ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassAccretionHistoryHearin2021} dark matter halo mass accretion history class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: powerLawIndexEarly, powerLawIndexLate, &
         &                                                                           rateRollOver      , timeMaximum

    !![
    <inputParameter>
      <name>powerLawIndexEarly</name>
      <description>The early time power law index.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>powerLawIndexLate</name>
      <description>The late time power law index.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rateRollOver</name>
      <description>The roll over rate parameter, $k$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeMaximum</name>
      <description>The time of the maximum mass.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=darkMatterHaloMassAccretionHistoryHearin2021(powerLawIndexEarly,powerLawIndexLate,rateRollOver,timeMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hearin2021ConstructorParameters

  function hearin2021ConstructorInternal(powerLawIndexEarly,powerLawIndexLate,rateRollOver,timeMaximum) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryHearin2021} dark matter halo mass accretion history class.
    !!}
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021)                :: self
    double precision                                              , intent(in   ) :: powerLawIndexEarly, powerLawIndexLate, &
         &                                                                           rateRollOver      , timeMaximum
    !![
    <constructorAssign variables="powerLawIndexEarly, powerLawIndexLate, rateRollOver, timeMaximum"/>
    !!]
    
    return
  end function hearin2021ConstructorInternal

  double precision function hearin2021PowerLawIndexEarly(self,node)
    !!{
    Return the early power law index for the given node.
    !!}
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    hearin2021PowerLawIndexEarly=self%powerLawIndexEarly
    return
  end function hearin2021PowerLawIndexEarly
  
  double precision function hearin2021PowerLawIndexLate(self,node)
    !!{
    Return the late power law index for the given node.
    !!}
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    hearin2021PowerLawIndexLate=self%powerLawIndexLate
    return
  end function hearin2021PowerLawIndexLate
  
  double precision function hearin2021RateRollOver(self,node)
    !!{
    Return the roll-over rate for the given node.
    !!}
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    hearin2021RateRollOver=self%rateRollOver
    return
  end function hearin2021RateRollOver
  
  double precision function hearin2021TimeMaximum(self,node)
    !!{
    Return the time of maximum mass for the given node.
    !!}
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    hearin2021TimeMaximum=self%timeMaximum
    return
  end function hearin2021TimeMaximum
  
  double precision function hearin2021Mass(self,node,time)
    !!{
    Compute the mass corresponding to {\normalfont \ttfamily time} in the mass accretion history of {\normalfont \ttfamily
    node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout), target :: node
    double precision                                              , intent(in   )         :: time
    class           (nodeComponentBasic                          ), pointer               :: basic
    
    basic                  =>  node%basic()
    hearin2021Mass         =  +10.0d0**(                                     &
         &                              +self%massMaximum        (node     ) &
         &                              +self%powerLawIndex      (node,time) &
         &                              *log10(                              &
         &                                     +     time                    &
         &                                     /self%timeMaximum_(node     ) &
         &                                    )                              &
         &                              )
    return
  end function hearin2021Mass

  double precision function hearin2021MassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of
    {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time
    class           (nodeComponentBasic                          ), pointer       :: basic

    basic                       =>  node%basic(         )
    hearin2021MassAccretionRate =  +self%mass (node,time)                     &
         &                         *(                                         &
         &                           +self%powerLawIndex          (node,time) &
         &                           /                                  time  &
         &                           +self%powerLawIndexDerivative(node,time) &
         &                           *log(                                    &
         &                                +     time                          &
         &                                /self%timeMaximum_      (node     ) &
         &                               )                                    &
         &                          )
    return
  end function hearin2021MassAccretionRate

  double precision function hearin2021MassMaximum(self,node)
    !!{
    Compute the maximum mass in the mass accretion history of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    class(nodeComponentBasic                          ), pointer       :: basic

    basic                 =>  node%basic()
    hearin2021MassMaximum =  +log10(                                            &
         &                          +                 basic%mass        (    )  &
         &                         )                                            &
         &                   -self%powerLawIndex(node,basic%time        (    )) &
         &                   *log10(                                            &
         &                          +                 basic%time        (    )  &
         &                          /                 self %timeMaximum_(node)  &
         &                         )
    return
  end function hearin2021MassMaximum

  double precision function hearin2021PowerLawIndex(self,node,time)
    !!{
    Compute the power-law index.
    !!}
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time

    hearin2021PowerLawIndex=self%sigmoid(                                &
         &                               log10(time)                   , &
         &                               self%timeZeroLogarithmic(node), &
         &                               self%rateRollOver_      (node), &
         &                               self%powerLawIndexEarly_(node), &
         &                               self%powerLawIndexLate_ (node)  &
         &                              )
    return
  end function hearin2021PowerLawIndex
  
  double precision function hearin2021PowerLawIndexDerivative(self,node,time)
    !!{
    Compute the derivative of the power-law index with respect to time.
    !!}
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time

    hearin2021PowerLawIndexDerivative=+(                                                                          &
         &                              +    self%powerLawIndexLate_ (node)                                       &
         &                              -    self%powerLawIndexEarly_(node)                                       &
         &                             )                                                                          &
         &                            *      self%rateRollOver_      (node)                                       &
         &                            *           time**(-1.0d0      +self%rateRollOver_      (node)/log(10.0d0)) &
         &                            *  exp(self%rateRollOver_(node)*self%timeZeroLogarithmic(node)            ) &
         &                            /(                                                                          &
         &                              +exp(self%rateRollOver_(node)*self%timeZeroLogarithmic(node)            ) &
         &                              +         time**(            +self%rateRollOver_      (node)/log(10.0d0)) &
         &                             )**2
    return
  end function hearin2021PowerLawIndexDerivative
  
  double precision function hearin2021TimeZeroLogarithmic(self,node)
    !!{
    Compute the $t_0$ parameter.
    !!}
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , parameter     :: timeZeroLogarithmicMinimum     =-0.175d0, timeZeroLogarithmicMaximum=+0.085d0, &
         &                                                                           timeZeroLogarithmicRolloverRate=+5.000d0, powerLawIndexEarlyZero    =+2.850d0

    hearin2021TimeZeroLogarithmic=+self%sigmoid(                                            &
         &                                      self%powerLawIndexEarly_            (node), &
         &                                           powerLawIndexEarlyZero               , &
         &                                           timeZeroLogarithmicRolloverRate      , &
         &                                           timeZeroLogarithmicMinimum           , &
         &                                           timeZeroLogarithmicMaximum             &
         &                                     )
    return
  end function hearin2021TimeZeroLogarithmic
  
  double precision function hearin2021Sigmoid(x,x0,k,yMinimum,yMaximum)
    !!{
    Sigmoid interpolation function.
    !!}
    implicit none
    double precision, intent(in   ) :: x       , x0      , &
         &                             yMinimum, yMaximum, &
         &                             k

    hearin2021Sigmoid=+  yMinimum  &
         &            +(           &
         &              +yMaximum  &
         &              -yMinimum  &
         &             )           &
         &            /(           &
         &              +1.0d0     &
         &              +exp(      &
         &                   -k    &
         &                   *(    &
         &                     +x  &
         &                     -x0 &
         &                    )    &
         &                  )      &
         &             )
    return
  end function hearin2021Sigmoid
