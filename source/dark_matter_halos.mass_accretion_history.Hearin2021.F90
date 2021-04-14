!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% An implementation of dark matter halo mass accretion histories using Andrew Hearin's rolling power-law model.

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryHearin2021">
  !#  <description>Dark matter halo mass accretion histories using Andrew Hearin's rolling power-law model.</description>
  !# </darkMatterHaloMassAccretionHistory>
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryHearin2021
     !% A dark matter halo mass accretion history class using Andrew Hearin's rolling power-law model.
     private
     double precision :: powerLawIndexEarly, powerLawIndexLate, &
          &              rateRollOver      , timeMaximum
   contains
     !# <methods>
     !#  <method description="Return the power law index at the given time."                      method="powerLawIndex"          />
     !#  <method description="Return the derivative of the power law index with respect to time." method="powerLawIndexDerivative"/>
     !#  <method description="Return the $\log_{10}(t_0)$ parameter."                             method="timeZeroLogarithmic"    />
     !#  <method description="Return the maximum mass in the mass accretion history."             method="massMaximum"            />
     !# </methods>
     procedure :: mass                    => hearin2021Mass
     procedure :: massAccretionRate       => hearin2021MassAccretionRate
     procedure :: powerLawIndex           => hearin2021PowerLawIndex
     procedure :: powerLawIndexDerivative => hearin2021PowerLawIndex
     procedure :: massMaximum             => hearin2021MassMaximum
     procedure :: timeZeroLogarithmic     => hearin2021TimeZeroLogarithmic
  end type darkMatterHaloMassAccretionHistoryHearin2021

  interface darkMatterHaloMassAccretionHistoryHearin2021
     !% Constructors for the {\normalfont \ttfamily hearin2021} dark matter halo mass accretion history class.
     module procedure hearin2021ConstructorParameters
     module procedure hearin2021ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryHearin2021

contains

  function hearin2021ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily hearin2021} dark matter halo mass accretion history class which takes a parameter
    !% set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: powerLawIndexEarly, powerLawIndexLate, &
         &                                                                           rateRollOver      , timeMaximum

    !# <inputParameter>
    !#   <name>powerLawIndexEarly</name>
    !#   <description>The early time power law index.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>powerLawIndexLate</name>
    !#   <description>The late time power law index.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rateRollOver</name>
    !#   <description>The roll over rate parameter, $k$.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeMaximum</name>
    !#   <description>The time of the maximum mass.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=darkMatterHaloMassAccretionHistoryHearin2021(powerLawIndexEarly, powerLawIndexLate, rateRollOver, timeMaximum)
    !# <inputParametersValidate source="parameters"/>
    return
  end function hearin2021ConstructorParameters

  function hearin2021ConstructorInternal(powerLawIndexEarly, powerLawIndexLate, rateRollOver, timeMaximum) result(self)
    !% Internal constructor for the {\normalfont \ttfamily hearin2021} dark matter halo mass accretion history class.
    implicit none
    type            (darkMatterHaloMassAccretionHistoryHearin2021)                :: self
    double precision                                              , intent(in   ) :: powerLawIndexEarly, powerLawIndexLate, &
         &                                                                           rateRollOver      , timeMaximum
    !# <constructorAssign variables="powerLawIndexEarly, powerLawIndexLate, rateRollOver, timeMaximum"/>
    
    return
  end function hearin2021ConstructorInternal

  double precision function hearin2021Mass(self,node,time)
    !% Compute the mass corresponding to {\normalfont \ttfamily time} in the mass accretion history of {\normalfont \ttfamily
    !% node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout), target :: node
    double precision                                              , intent(in   )         :: time
    class           (nodeComponentBasic                          ), pointer               :: basic
    
    basic                  =>  node%basic()
    hearin2021Mass         =  +10.0d0**(                          &
         &                              +self%massMaximum  (node) &
         &                              +self%powerLawIndex(time) &
         &                              *log10(                   &
         &                                     +     time         &
         &                                     /self%timeMaximum  &
         &                                    )                   &
         &                              )
    return
  end function hearin2021Mass

  double precision function hearin2021MassAccretionRate(self,node,time)
    !% Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of
    !% {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time
    class           (nodeComponentBasic                          ), pointer       :: basic

    basic                       =>  node%basic(         )
    hearin2021MassAccretionRate =  +self%mass (node,time)                &
         &                         *(                                    &
         &                           +self%powerLawIndex          (time) &
         &                           /                             time  &
         &                           +self%powerLawIndexDerivative(time) &
         &                           *log(                               &
         &                                +     time                     &
         &                                /self%timeMaximum              &
         &                               )                               &
         &                          )
    return
  end function hearin2021MassAccretionRate

  double precision function hearin2021MassMaximum(self,node)
    !% Compute the maximum mass in the mass accretion history of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    class(nodeComponentBasic                          ), pointer       :: basic

    basic                 =>  node%basic()
    hearin2021MassMaximum =  +log10(                                &
         &                          +            basic%mass()       &
         &                         )                                &
         &                   -self%powerLawIndex(basic%time()     ) &
         &                   *log10(                                &
         &                          +            basic%time()       &
         &                          /            self %timeMaximum  &
         &                         )
    return
  end function hearin2021MassMaximum

  double precision function hearin2021PowerLawIndex(self,time)
    !% Compute the power-law index.
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    hearin2021PowerLawIndex=sigmoid(                                      &
         &                          log10(time)                         , &
         &                          self%timeZeroLogarithmic()          , &
         &                          self%rateRollOver                   , &
         &                          self%powerLawIndexEarly             , &
         &                          self%powerLawIndexLate                &
         &                         )
    return
  end function hearin2021PowerLawIndex
  
  double precision function hearin2021PowerLawIndexDerivative(self,time)
    !% Compute the derivative of the power-law index with respect to time.
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    hearin2021PowerLawIndexDerivative=+(                                                       &
         &                              +    self%powerLawIndexLate                            &
         &                              -    self%powerLawIndexEarly                           &
         &                             )                                                       &
         &                            *      self%rateRollOver                                 &
         &                            *           time**(-1.0d0+self%rateRollOver/log(10.0d0)) &
         &                            *  exp(self%rateRollOver*self%timeZeroLogarithmic())     &
         &                            /(                                                       &
         &                              +exp(self%rateRollOver*self%timeZeroLogarithmic())     &
         &                              +         time**(      +self%rateRollOver/log(10.0d0)) &
         &                             )**2
    return
  end function hearin2021PowerLawIndexDerivative
  
  double precision function hearin2021TimeZeroLogarithmic(self)
    !% Compute the $t_0$ parameter.
    implicit none
    class           (darkMatterHaloMassAccretionHistoryHearin2021), intent(inout) :: self
    double precision                                              , parameter     :: timeZeroLogarithmicMinimum     =-0.175d0, timeZeroLogarithmicMaximum=+0.085d0, &
         &                                                                           timeZeroLogarithmicRolloverRate=+5.000d0, powerLawIndexEarlyZero    =+2.850d0

    hearin2021TimeZeroLogarithmic=+sigmoid(                                      &
         &                                 self%powerLawIndexEarly             , &
         &                                      powerLawIndexEarlyZero         , &
         &                                      timeZeroLogarithmicRolloverRate, &
         &                                      timeZeroLogarithmicMinimum     , &
         &                                      timeZeroLogarithmicMaximum       &
         &                                )
    return
  end function hearin2021TimeZeroLogarithmic
  
  double precision function sigmoid(x,x0,k,yMinimum,yMaximum)
    !% Sigmoid interpolation function.
    implicit none
    double precision, intent(in   ) :: x       , x0      , &
         &                             yMinimum, yMaximum, &
         &                             k

    sigmoid=+  yMinimum  &
         &  +(           &
         &    +yMaximum  &
         &    -yMinimum  &
         &   )           &
         &  /(           &
         &    +1.0d0     &
         &    +exp(      &
         &         -k    &
         &         *(    &
         &           +x  &
         &           -x0 &
         &          )    &
         &        )      &
         &   )
    return
  end function sigmoid
