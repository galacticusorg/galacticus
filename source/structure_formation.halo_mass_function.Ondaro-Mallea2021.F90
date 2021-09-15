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

!!{ Contains a module which implements a dark matter halo mass function class for non-universal primordial power spectra.
!!}
  use :: Cosmological_Density_Field    , only : cosmologicalMassVarianceClass
  use :: Excursion_Sets_First_Crossings, only : excursionSetFirstCrossingClass
  use :: Linear_Growth                 , only: linearGrowthClass

!![
 <haloMassFunction name="haloMassFunctionOndaroMallea2021">
    <description>
     A dark matter halo mass function class using the function given by \cite{press_formation_1974}. Specifically,
     \begin{equation}
     n(M,t) = 2 {\Omega_\mathrm{M} \rho_\mathrm{crit} \over M^2} \alpha \sigma^2(M) f[S(M,t)],
     \end{equation}
     where $\alpha = \mathrm{d}\ln\sigma/\mathrm{d}\ln M$ and $f[S]$ is the excursion set barrier first crossing distribution
     for variance $S(M)=\sigma^2(M)$, computed using the selected \refClass{excursionSetFirstCrossingClass}.
    </description>
   </haloMassFunction>
!!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionOndaroMallea2021
     !!{
      A halo mass function class modifying the primordial power spectrum.
     !!}
     private
     double precision                                    :: n_0, n_1, n_2, a_0, a_1
     class(cosmologicalMassVarianceClass ), pointer :: cosmologicalMassVariance_  => null()
     class(excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_ => null()
     class(linearGrowthClass             ), pointer :: linearGrowth_ => null()
     class           (haloMassFunctionClass   ), pointer :: haloMassFunctionShethTormen  => null()
    contains
     final     ::                 OndaroMallea2021Destructor
     procedure :: differential => OndaroMallea2021Differential
  end type haloMassFunctionOndaroMallea2021

  interface haloMassFunctionOndaroMallea2021
     !!{
 Constructors for the primordial power spectrum halo mass function class.
     !!}
     module procedure OndaroMallea2021ConstructorParameters
     module procedure OndaroMallea2021ConstructorInternal
  end interface haloMassFunctionOndaroMallea2021

contains

  function OndaroMallea2021ConstructorParameters(parameters) result(self)
    !!{
 Constructor for the primordial power halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionOndaroMallea2021)               :: self
    type (inputParameters               ), intent(inout) :: parameters
    class           (haloMassFunctionClass   ), pointer  :: haloMassFunctionShethTormen
    class(cosmologicalMassVarianceClass ), pointer       :: cosmologicalMassVariance_
    class(excursionSetFirstCrossingClass), pointer       :: excursionSetFirstCrossing_
    class(cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class(linearGrowthClass             ), pointer       :: linearGrowth_
    double precision                                     :: n_0, n_1, n_2, a_0, a_1

!![
 <objectBuilder class="cosmologyParameters"       name="cosmologyParameters_"       source="parameters"/>
     <objectBuilder class="cosmologicalMassVariance"  name="cosmologicalMassVariance_"  source="parameters"/>
     <objectBuilder class="excursionSetFirstCrossing" name="excursionSetFirstCrossing_" source="parameters"/>
     <objectBuilder class="linearGrowth"              name="linearGrowth_"              source="parameters"/>
     <objectBuilder class="haloMassFunction"    name="haloMassFunctionShethTormen"  source="parameters"/>
  <inputParameter>
       <name>n_0</name>
       <source>parameters</source>
       <defaultValue>0.707d0</defaultValue>
       <description>The parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     </inputParameter>
     <inputParameter>
       <name>n_1</name>
       <source>parameters</source>
       <defaultValue>0.3d0</defaultValue>
       <description>The parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     </inputParameter>
     <inputParameter>
       <name>n_2</name>
       <source>parameters</source>
       <defaultValue>0.3221836349d0</defaultValue>
       <description>The normalization parameter $A$ in the halo mass function fit.</description>
     </inputParameter>
     <inputParameter>
       <name>a_0</name>
       <source>parameters</source>
       <defaultValue>0.707d0</defaultValue>
       <description>The parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     </inputParameter>
     <inputParameter>
       <name>a_1</name>
       <source>parameters</source>
       <defaultValue>0.3d0</defaultValue>
       <description>The parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
     </inputParameter>
!!]
    self=haloMassFunctionOndaroMallea2021(haloMassFunctionShethTormen,n_0,n_1,n_2,a_0,a_1,cosmologyParameters_,cosmologicalMassVariance_,excursionSetFirstCrossing_,linearGrowth_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"      />
    !# <objectDestructor name="cosmologicalMassVariance_" />
    !# <objectDestructor name="excursionSetFirstCrossing_"/>
    !# <objectDestructor name="linearGrowth_"             />
    !# <objectDestructor name="haloMassFunctionShethTormen"/>
    return
  end function OndaroMallea2021ConstructorParameters

  function OndaroMallea2021ConstructorInternal(haloMassFunctionShethTormen,n_0,n_1,n_2,a_0,a_1,cosmologyParameters_,cosmologicalMassVariance_,excursionSetFirstCrossing_,linearGrowth_) result(self)
    !!{
 Internal constructor for the primordial power halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionOndaroMallea2021)                       :: self
    class(cosmologyParametersClass      ), target, intent(in   ) :: cosmologyParameters_
    class(cosmologicalMassVarianceClass ), target, intent(in   ) :: cosmologicalMassVariance_
    class(excursionSetFirstCrossingClass), target, intent(in   ) :: excursionSetFirstCrossing_
    class(linearGrowthClass             ), target, intent(in   ) :: linearGrowth_
    class(haloMassFunctionClass   ), target, intent(in   )       :: haloMassFunctionShethTormen
    double precision               , intent(in   )               :: n_0, n_1, n_2, a_0, a_1
!![
 <constructorAssign variables="*haloMassFunctionShethTormen,n_0,n_1,n_2,a_0,a_1,*cosmologyParameters_, *cosmologicalMassVariance_, *excursionSetFirstCrossing_,*linearGrowth_"/>
!!]

    return
  end function OndaroMallea2021ConstructorInternal

  subroutine OndaroMallea2021Destructor(self)
    !!{
 Destructor for the primordial power halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionOndaroMallea2021), intent(inout) :: self

    !![
     <objectDestructor name="self%cosmologyParameters_"       />
     <objectDestructor name="self%cosmologicalMassVariance_"  />
     <objectDestructor name="self%excursionSetFirstCrossing_" />
     <objectDestructor name="self%haloMassFunctionShethTormen" />
     <objectDestructor name="self%linearGrowth_" />
    !!]
    return
  end subroutine OndaroMallea2021Destructor

  double precision function OndaroMallea2021Differential(self,time,mass,node)
    !!{
 Return the differential halo mass function at the given time and mass.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (haloMassFunctionOndaroMallea2021), intent(inout)          :: self
    double precision                                , intent(in   )           :: time , mass
    type            (treeNode                      ), intent(inout), optional :: node
    double precision                                                          :: alpha, variance, neff, aeff, f_2, f_3

    if (.not.present(node)) call Galacticus_Error_Report('"node" must be present'//{introspection:location})
    alpha                     =abs(self%cosmologicalMassVariance_ %rootVarianceLogarithmicGradient(mass,time))
    variance                  =    self%cosmologicalMassVariance_ %rootVariance                   (mass,time) **2
    OndaroMallea2021Differential=+self%haloMassFunctionShethTormen%Differential(time,mass,node)
    aeff = self%linearGrowth_%logarithmicDerivativeExpansionFactor(time)
    neff = -3 - 6 * self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass, time)
    f_2 = self%n_0 * neff**2 + self%n_1 * neff + self%n_2
    f_3 = self%a_0 * aeff + self%a_1
    if (variance > 0.0d0) then
      OndaroMallea2021Differential = OndaroMallea2021Differential * f_2 * f_3     
    else
       OndaroMallea2021Differential=+0.0d0
    end if
    return
  end function OndaroMallea2021Differential
