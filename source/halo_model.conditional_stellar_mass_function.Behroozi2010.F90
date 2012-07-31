!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu

!% Contains a module which implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional stellar mass function.

module Conditional_Stellar_Mass_Functions_Behroozi2010
  !% Implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional stellar mass function.
  use FGSL
  private
  public :: Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize

  ! Parameters of the Behroozi et al. fitting function.
  double precision                            :: conditionalStellarMassFunctionBehrooziAlphaSatellite
  double precision                            :: conditionalStellarMassFunctionBehrooziLog10M1
  double precision                            :: conditionalStellarMassFunctionBehrooziLog10Mstar0
  double precision                            :: conditionalStellarMassFunctionBehrooziBeta
  double precision                            :: conditionalStellarMassFunctionBehrooziDelta
  double precision                            :: conditionalStellarMassFunctionBehrooziGamma
  double precision                            :: conditionalStellarMassFunctionBehrooziSigmaLogMstar
  double precision                            :: conditionalStellarMassFunctionBehrooziBCut
  double precision                            :: conditionalStellarMassFunctionBehrooziBSatellite
  double precision                            :: conditionalStellarMassFunctionBehrooziBetaCut
  double precision                            :: conditionalStellarMassFunctionBehrooziBetaSatellite
  double precision                            :: Mstar0

  ! Tablulation of stellar-halo mass relation.
  integer,          parameter                 :: fMassStellarTablePointsPerDecade=10
  integer                                     :: fMassStellarTableCount
  double precision                            :: fMassStellarTableMinimum=1.0d8,fMassStellarTableMaximum=1.0d13
  double precision                            :: fMassHaloTableMinimum,fMassHaloTableMaximum
  double precision, allocatable, dimension(:) :: fMassStellarTable,fMassHaloTable
  logical                                     :: interpolationReset=.false.
  type(fgsl_interp)                           :: interpolationObject
  type(fgsl_interp_accel)                     :: interpolationAccelerator

contains

  !# <conditionalStellarMassFunctionMethod>
  !#  <unitName>Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize</unitName>
  !# </conditionalStellarMassFunctionMethod>
  subroutine Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize(conditionalStellarMassFunctionMethod&
       &,Cumulative_Conditional_Stellar_Mass_Function_Get,Cumulative_Conditional_Stellar_Mass_Function_Var_Get)
    !% Initializes the ``Behroozi2010'' conditional stellar mass function method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: conditionalStellarMassFunctionMethod
    procedure(double precision), pointer, intent(inout) :: Cumulative_Conditional_Stellar_Mass_Function_Get,Cumulative_Conditional_Stellar_Mass_Function_Var_Get
    
    if (conditionalStellarMassFunctionMethod == 'Behroozi2010') then
       Cumulative_Conditional_Stellar_Mass_Function_Get     => Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010
       Cumulative_Conditional_Stellar_Mass_Function_Var_Get => Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010

       ! Get parameters of the fitting function.
       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziAlphaSatellite</name>
       !@   <defaultValue>1.0 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\alpha_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziAlphaSatellite',conditionalStellarMassFunctionBehrooziAlphaSatellite,defaultValue=1.0d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziLog10M1</name>
       !@   <defaultValue>12.520 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_1$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziLog10M1',conditionalStellarMassFunctionBehrooziLog10M1,defaultValue=12.520d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziLog10Mstar0</name>
       !@   <defaultValue>10.916 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_{\star,0}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziLog10Mstar0',conditionalStellarMassFunctionBehrooziLog10Mstar0,defaultValue=10.916d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBeta</name>
       !@   <defaultValue>0.457 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBeta',conditionalStellarMassFunctionBehrooziBeta,defaultValue=0.457d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziDelta</name>
       !@   <defaultValue>0.5666 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\delta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziDelta',conditionalStellarMassFunctionBehrooziDelta,defaultValue=0.5666d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziGamma</name>
       !@   <defaultValue>1.53 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\gamma$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziGamma',conditionalStellarMassFunctionBehrooziGamma,defaultValue=1.53d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziSigmaLogMstar</name>
       !@   <defaultValue>0.206 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\sigma_{\log M_\star}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziSigmaLogMstar',conditionalStellarMassFunctionBehrooziSigmaLogMstar,defaultValue=0.206d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBCut</name>
       !@   <defaultValue>1.47 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBCut',conditionalStellarMassFunctionBehrooziBCut,defaultValue=1.47d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBSatellite</name>
       !@   <defaultValue>10.62 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBSatellite',conditionalStellarMassFunctionBehrooziBSatellite,defaultValue=10.62d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBetaCut</name>
       !@   <defaultValue>$-0.13$ (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBetaCut',conditionalStellarMassFunctionBehrooziBetaCut,defaultValue=-0.13d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBetaSatellite</name>
       !@   <defaultValue>0.859 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBetaSatellite',conditionalStellarMassFunctionBehrooziBetaSatellite,defaultValue=0.859d0)

       ! Compute the actual value of the Mstar0 parameter.
       Mstar0=10.0d0**conditionalStellarMassFunctionBehrooziLog10Mstar0
    end if
    return
  end subroutine Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize

  double precision function Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010(massHalo,massStellar)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in) :: massHalo,massStellar
    double precision             :: numberCentrals,numberSatellites

    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellar,numberCentrals,numberSatellites)

    ! Sum central and satellite contributions.
    Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010=numberCentrals+numberSatellites
    return
  end function Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010

  double precision function Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010(massHalo,massStellarLow,massStellarHigh)
    !% Computes the variance in the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv
    !% \phi(M_\star|M_{\rm halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}. Assumes that the number of
    !% satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that
    !% the numbers of satellites and centrals are uncorrelated.
    implicit none
    double precision, intent(in) :: massHalo,massStellarLow,massStellarHigh
    double precision             :: numberCentrals,numberSatellites,numberCentralsLow,numberSatellitesLow,numberCentralsHigh&
         &,numberSatellitesHigh

    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellarLow ,numberCentralsLow ,numberSatellitesLow )
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellarHigh,numberCentralsHigh,numberSatellitesHigh)
    numberSatellites=max(numberSatellitesLow-numberSatellitesHigh,0.0d0)
    numberCentrals  =max(numberCentralsLow  -numberCentralsHigh  ,0.0d0)

    ! Compute the variance (using the Bienaymé formula for uncorrelated variables).
    Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010= &
         &  numberSatellites                                       & ! Satellites are Poisson distributed, so the variance is just their number.
         & +numberCentrals*(1.0d0-numberCentrals)                    ! Centrals are Bernoulli distributed.
    return
  end function Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010

  subroutine Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellar,numberCentrals,numberSatellites)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    implicit none
    double precision, intent(in ) :: massHalo,massStellar
    double precision, intent(out) :: numberCentrals,numberSatellites
    double precision, parameter   :: massNormalization=1.0d12
    double precision              :: fMassHalo,fMassStellar,massCut,massSatellite

    ! Ensure that the stellar-halo mass relation is tabulated over a sufficient range.
    do while (massHalo < fMassHaloTableMinimum .or. massHalo > fMassHaloTableMaximum)
       if (massHalo < fMassHaloTableMinimum) fMassStellarTableMinimum=0.5d0*fMassStellarTableMinimum
       if (massHalo > fMassHaloTableMaximum) fMassStellarTableMaximum=2.0d0*fMassStellarTableMaximum
       fMassStellarTableCount=int(log10(fMassStellarTableMaximum/fMassStellarTableMinimum)*fMassStellarTablePointsPerDecade)+1
       if (allocated(fMassStellarTable)) call Dealloc_Array(fMassStellarTable)
       if (allocated(fMassHaloTable   )) call Dealloc_Array(fMassHaloTable   )
       call Alloc_Array(fMassStellarTable,[fMassStellarTableCount])
       call Alloc_Array(fMassHaloTable   ,[fMassStellarTableCount])
       fMassStellarTable    =Make_Range(fMassStellarTableMinimum,fMassStellarTableMaximum,fMassStellarTableCount,rangeType=rangeTypeLogarithmic)
       fMassHaloTable       =fSHMRInverse(fMassStellarTable)
       fMassHaloTableMinimum=fMassHaloTable(                     1)
       fMassHaloTableMaximum=fMassHaloTable(fMassStellarTableCount)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
       interpolationReset=.true.
    end do

    ! Compute the forward and inverse stellar-halo mass relation.
    fMassHalo   =fSHMRInverse(massStellar)
    fMassStellar=Interpolate(fMassStellarTableCount,fMassHaloTable,fMassStellarTable,interpolationObject,interpolationAccelerator,massHalo&
       &,extrapolationType=extrapolationTypeNone,reset=interpolationReset)

    ! Compute mass scales.
    massSatellite=massNormalization*conditionalStellarMassFunctionBehrooziBSatellite*(fMassHalo/massNormalization)**conditionalStellarMassFunctionBehrooziBetaSatellite
    massCut      =massNormalization*conditionalStellarMassFunctionBehrooziBCut      *(fMassHalo/massNormalization)**conditionalStellarMassFunctionBehrooziBetaCut
    ! Compute the number of central galaxies.
    numberCentrals=                                                            &
         &          0.5d0                                                      &
         &         *(                                                          &
         &           1.0d0                                                     &
         &           -erf(                                                     &
         &                 log10(massStellar/fMassStellar)                     &
         &                /dsqrt(2.0d0)                                        &
         &                /conditionalStellarMassFunctionBehrooziSigmaLogMstar &
         &               )                                                     &
         &          )

    ! Compute the number of satellites.
    numberSatellites=                                                                                &
         &            numberCentrals                                                                 &
         &           *(massHalo/massSatellite)**conditionalStellarMassFunctionBehrooziAlphaSatellite &
         &           *exp(-massCut/massHalo)
    return
  end subroutine Cumulative_Conditional_Stellar_Mass_Function_Compute

  elemental double precision function fSHMRInverse(massStellar)
    !% The median stellar mass vs. halo mass relation functional form from \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in) :: massStellar
    double precision, parameter  :: logHaloMassTransition=20.0d0
    double precision             :: argument

    ! Compute the logarithmic halo mass for the given stellar mass.
    argument=                                                                                                                   &
         &    conditionalStellarMassFunctionBehrooziLog10M1                                                                     &
         &   +conditionalStellarMassFunctionBehrooziBeta*(log10(massStellar)-conditionalStellarMassFunctionBehrooziLog10Mstar0) &
         &   +(massStellar/Mstar0)**conditionalStellarMassFunctionBehrooziDelta                                                 &
         &   /(1.0d0+1.0d0/(massStellar/Mstar0)**conditionalStellarMassFunctionBehrooziGamma)                                   &
         &   -0.5d0
    ! For some parameter choices, teh halo mass can grow unreasonably large. Therefore, above a transition value, allow the
    ! logarithmic halo mass to grow only logarithmically. Halo masses this high are irrelevant anyway (since the halo mass
    ! function will be so suppressed, while allowing the mass to continue to slowly grow allows for any tabulation to remain
    ! monotonically growing.
    if (argument > logHaloMassTransition) argument=logHaloMassTransition+log(argument/logHaloMassTransition)
    ! Compute the halo mass.
    fSHMRInverse=10.0d0**min(argument,100.0d0)
    return
  end function fSHMRInverse

end module Conditional_Stellar_Mass_Functions_Behroozi2010
