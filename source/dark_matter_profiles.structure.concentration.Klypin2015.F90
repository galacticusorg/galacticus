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
  An implementation of dark matter halo profile concentrations using the \cite{klypin_multidark_2014} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Tables                    , only : table1DGeneric

  ! Labels for virial density contrast definition.
  !![
  <enumeration>
   <name>klypin2015DensityContrast</name>
   <description>Enumeration of density contrasts in the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.</description>
   <visibility>private</visibility>
   <entry label="fixed" />
   <entry label="virial"/>
  </enumeration>
  !!]

  ! Labels for fitting function type.
  !![
  <enumeration>
   <name>klypin2015FittingFunction</name>
   <description>Enumeration of fitting functions in the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.</description>
   <visibility>private</visibility>
   <entry label="eqn24"/>
   <entry label="eqn25"/>
  </enumeration>
  !!]

  ! Labels for sample selection.
  !![
  <enumeration>
   <name>klypin2015Sample</name>
   <description>Enumeration of sample choices available in the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.</description>
   <visibility>private</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="planck200CritRelaxedMass"   />
   <entry label="planck200CritAllMass"       />
   <entry label="planck200CritRelaxedVmax"   />
   <entry label="planck200CritAllVmax"       />
   <entry label="planckVirialRelaxedMass"    />
   <entry label="planckVirialAllMass"        />
   <entry label="planckVirialRelaxedVmax"    />
   <entry label="planckVirialAllVmax"        />
   <entry label="wmap7200CritRelaxedMass"    />
   <entry label="wmap7200CritAllMass"        />
   <entry label="wmap7200CritRelaxedVmax"    />
   <entry label="wmap7VirialRelaxedMass"     />
   <entry label="wmap7VirialAllMass"         />
   <entry label="wmap7VirialRelaxedVmax"     />
   <entry label="planck200CritAllMassUni"    />
   <entry label="planck200CritRelaxedMassUni"/>
   <entry label="planckVirialAllMassUni"     />
   <entry label="planckVirialRelaxedMassUni" />
  </enumeration>
  !!]

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationKlypin2015">
   <description>Dark matter halo concentrations are computed using the algorithm of \cite{klypin_multidark_2014}.</description>
   <deepCopy>
    <functionClass variables="darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationKlypin2015
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{klypin_multidark_2014}.
     !!}
     private
     class(cosmologyFunctionsClass                 ), pointer :: cosmologyFunctions_              => null()
     class(cosmologyParametersClass                ), pointer :: cosmologyParameters_             => null()
     class(cosmologicalMassVarianceClass           ), pointer :: cosmologicalMassVariance_        => null()
     class(virialDensityContrastClass              ), pointer :: virialDensityContrastDefinition_ => null()
     type (darkMatterProfileDMONFW                 ), pointer :: darkMatterProfileDMODefinition_  => null()
     type (enumerationKlypin2015DensityContrastType)          :: virialDensityContrast
     type (enumerationKlypin2015FittingFunctionType)          :: fittingFunction
     type (enumerationKlypin2015SampleType         )          :: sample
     type (table1DGeneric                          )          :: fitParameters
   contains
     final     ::                                   klypin2015Destructor
     procedure :: concentration                  => klypin2015Concentration
     procedure :: densityContrastDefinition      => klypin2015DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => klypin2015DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationKlypin2015

  interface darkMatterProfileConcentrationKlypin2015
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationKlypin2015} dark matter halo profile concentration class.
     !!}
     module procedure klypin2015ConstructorParameters
     module procedure klypin2015ConstructorInternal
  end interface darkMatterProfileConcentrationKlypin2015

contains

  function klypin2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.
    !!}
    implicit none
    type (darkMatterProfileConcentrationKlypin2015)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass           ), pointer       :: cosmologicalMassVariance_
    type (varying_string                          )                :: sample

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>sample</name>
      <source>parameters</source>
      <defaultValue>var_str('planck200CritRelaxedMass')</defaultValue>
      <description>The sample to use for the halo concentration algorithm of \cite{klypin_multidark_2014}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    ! Construct the object.
    self=darkMatterProfileConcentrationKlypin2015(enumerationKlypin2015SampleEncode(char(sample),includesPrefix=.false.),cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function klypin2015ConstructorParameters

  function klypin2015ConstructorInternal(sample,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationKlypin2015} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Table_Labels           , only : extrapolationTypeFix
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type (darkMatterProfileConcentrationKlypin2015          )                        :: self
    type (enumerationKlypin2015SampleType                   ), intent(in   )         :: sample
    class(cosmologyParametersClass                          ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass                           ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass                     ), intent(in   ), target :: cosmologicalMassVariance_
    type (darkMatterHaloScaleVirialDensityContrastDefinition), pointer               :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="sample, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_"/>
    !!]

    select case (sample%ID)
    case (klypin2015SamplePlanck200CritRelaxedMass%ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.750d0,6.700d0,6.250d0,5.020d0,4.190d0,3.300d0,3.000d0,2.720d0,2.400d0,2.100d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d-1,0.950d-1,0.920d-1,0.880d-1,0.850d-1,0.830d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[4.500d5,2.000d4,8.000d3,7.800d2,1.600d2,2.700d1,1.400d1,6.800d0,1.600d0,0.300d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllMass    %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.400d0,6.250d0,5.650d0,4.300d0,3.530d0,2.700d0,2.420d0,2.200d0,1.920d0,1.650d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.120d0,0.117d0,0.115d0,0.110d0,0.095d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,1.000d5,2.000d4,9.000d2,3.000d2,4.200d1,1.700d1,8.500d0,2.000d0,0.300d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritRelaxedVmax%ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[8.000d0,6.820d0,6.400d0,5.200d0,4.350d0,3.500d0,3.120d0,2.850d0,2.550d0,2.160d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d-1,0.950d-1,0.920d-1,0.880d-1,0.850d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d5,9.000d3,4.500d3,6.000d2,1.500d2,2.700d1,1.100d1,5.500d0,1.500d0,0.220d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllVmax    %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.750d0,6.500d0,5.950d0,4.550d0,3.680d0,2.750d0,2.500d0,2.250d0,2.050d0,1.760d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.115d0,0.115d0,0.115d0,0.110d0,0.105d0,1.000d-1,0.950d-1,0.900d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,1.800d4,6.000d3,6.000d2,1.500d2,2.000d1,1.000d1,5.000d0,1.500d0,0.250d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedMass %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.020d1,7.850d0,7.160d0,5.450d0,4.550d0,3.550d0,3.240d0,2.920d0,2.600d0,2.300d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.100d0,0.095d0,0.092d0,0.088d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d5,1.200d4,5.500d3,7.000d2,1.800d2,3.000d1,1.500d1,7.000d0,1.900d0,0.360d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllMass     %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.750d0,7.250d0,6.500d0,4.750d0,3.800d0,3.000d0,2.650d0,2.420d0,2.100d0,1.860d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.107d0,0.105d0,0.100d0,0.095d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.000d5,2.200d4,1.000d4,1.000d3,2.100d2,4.300d1,1.800d1,9.000d0,1.900d0,0.420d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedVmax %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.070d1,8.100d0,7.330d0,5.650d0,4.650d0,3.700d0,2.350d0,2.980d0,2.700d0,2.350d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.100d0,0.100d0,0.088d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.400d4,5.000d3,2.200d3,5.200d2,1.200d2,2.500d1,1.200d1,5.000d0,1.400d0,0.260d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllVmax     %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.030d1,7.600d0,6.830d0,4.960d0,3.960d0,3.000d0,2.730d0,2.450d0,2.240d0,2.030d0], &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.115d0,0.115d0,0.115d0,0.110d0,0.105d0,0.100d0,0.095d0,0.090d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[4.800d4,5.000d3,2.700d3,3.900d2,1.100d2,1.800d1,1.000d1,5.000d0,1.400d0,0.360d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritRelaxedMass  %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[6.900d0,5.700d0,4.550d0,3.750d0,2.900d0,2.600d0,2.400d0,2.200d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,6.000d3,5.000d2,1.000d2,2.000d1,1.000d1,6.000d0,3.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritAllMass      %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[6.600d0,5.250d0,3.850d0,3.000d0,2.100d0,1.800d0,1.600d0,1.400d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.105d0,0.103d0,0.097d0,0.095d0,0.095d0,0.095d0,0.095d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d6,6.000d4,8.000d2,1.100d2,1.300d1,6.000d0,3.000d0,1.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritRelaxedVmax  %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.200d0,5.900d0,4.700d0,3.850d0,3.000d0,2.700d0,2.500d0,2.300d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d5,4.000d3,4.000d2,8.000d1,1.300d1,7.000d0,3.500d0,2.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialRelaxedMass   %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.500d0,6.750d0,5.000d0,4.050d0,3.100d0,2.800d0,2.450d0,2.200d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[3.000d5,5.000d3,4.500d2,9.000d1,1.500d1,8.000d0,3.500d0,1.500d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialAllMass       %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.000d0,6.000d0,4.300d0,3.300d0,2.300d0,2.100d0,1.850d0,1.700d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.100d0,0.100d0,0.100d0,0.100d0,0.095d0,0.095d0,0.095d0,0.095d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d6,7.000d3,5.500d2,9.000d1,1.100d1,6.000d0,2.500d0,1.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialRelaxedVmax    %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn24
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.750d0,7.020d0,5.230d0,4.250d0,3.200d0,2.900d0,2.500d0,2.350d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.300d5,4.000d3,4.000d2,8.000d1,1.100d1,6.000d0,2.500d0,1.200d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllMassUni    %ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn25
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,2.890d0,5.410d0]                , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.400d0,0.650d0,0.820d0,1.080d0,1.230d0,1.600d0,1.680d0,1.700d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.278d0,0.375d0,0.411d0,0.436d0,0.426d0,0.375d0,0.360d0,0.351d0]                , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritRelaxedMassUni%ID)
       self%virialDensityContrast=klypin2015DensityContrastFixed
       self%fittingFunction      =klypin2015FittingFunctionEqn25
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,2.890d0,5.410d0]                , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.950d0,1.060d0,1.150d0,1.280d0,1.390d0,1.660d0,1.700d0,1.720d0]                , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.522d0,0.550d0,0.562d0,0.562d0,0.541d0,0.480d0,0.464d0,0.450d0]                , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllMassUni     %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn25
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,5.500d0]                        , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.750d0,0.900d0,0.970d0,1.120d0,1.280d0,1.520d0,1.620d0]                        , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.567d0,0.541d0,0.529d0,0.496d0,0.474d0,0.421d0,0.393d0]                        , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedMassUni %ID)
       self%virialDensityContrast=klypin2015DensityContrastVirial
       self%fittingFunction      =klypin2015FittingFunctionEqn25
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,5.500d0]                        , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=[extrapolationTypeFix,extrapolationTypeFix]                                        &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.990d0,1.100d0,1.160d0,1.290d0,1.410d0,1.650d0,1.720d0]                        , &
            &            table            =1                                                                                  &
            &           )
       call self                                                                                                              &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.716d0,0.673d0,0.660d0,0.615d0,0.585d0,0.518d0,0.476d0]                        , &
            &            table            =2                                                                                  &
            &           )
    end select
    ! Build definitions.
    select case (self%virialDensityContrast%ID)
    case (klypin2015DensityContrastFixed%ID)
       allocate(virialDensityContrastFixed                                      :: self%virialDensityContrastDefinition_)
       select type (virialDensityContrastDefinition_ => self%virialDensityContrastDefinition_)
       type is (virialDensityContrastFixed)
          !![
          <referenceConstruct object="virialDensityContrastDefinition_">
           <constructor>
            virialDensityContrastFixed                                    (                                                                &amp;
             &amp;                                                         densityContrastValue                =200.0d0                  , &amp;
             &amp;                                                         densityType                         =fixedDensityTypeCritical , &amp;
             &amp;                                                         turnAroundOverVirialRadius          =2.0d0                    , &amp;
             &amp;                                                         cosmologyParameters_                =self%cosmologyParameters_, &amp;
             &amp;                                                         cosmologyFunctions_                 =self%cosmologyFunctions_   &amp;
             &amp;                                                        )
           </constructor>
          </referenceConstruct>
          !!]
       end select
    case (klypin2015DensityContrastVirial%ID)
       allocate(virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt :: self%virialDensityContrastDefinition_)
       select type (virialDensityContrastDefinition_ => self%virialDensityContrastDefinition_)
       type is (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)
          !![
          <referenceConstruct object="virialDensityContrastDefinition_">
           <constructor>
            virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                &amp;
             &amp;                                                         tableStore                          =.true.                   , &amp;
             &amp;                                                         cosmologyFunctions_                 =self%cosmologyFunctions_   &amp;
             &amp;                                                        )
           </constructor>
          </referenceConstruct>
          !!]
       end select
    end select
    allocate(self%darkMatterProfileDMODefinition_)
    allocate(     darkMatterHaloScaleDefinition_ )
    !![
    <referenceConstruct              object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition(                                                                            &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                             virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                           (                                                                            &amp;
       &amp;                                             velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function klypin2015ConstructorInternal

  subroutine klypin2015Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationKlypin2015} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationKlypin2015), intent(inout) :: self

    call self%fitParameters%destroy()
    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine klypin2015Destructor

  double precision function klypin2015Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{klypin_multidark_2014} algorithm.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , treeNode
    implicit none
    class           (darkMatterProfileConcentrationKlypin2015), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBasic                      )               , pointer :: basic
    double precision                                          , parameter              :: massReference            =1.0d12
    double precision                                                                   :: massLittleH                     , concentration0, &
         &                                                                                mass0                           , gamma         , &
         &                                                                                redshift                        , a0            , &
         &                                                                                b0                              , sigma

    ! Get the basic component, and find the halo mass in "little-h" units.
    basic       =>  node                     %basic         (                        )
    massLittleH =  +basic                    %mass          (                        ) &
         &         *self%cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
    ! Find redshift.
    redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(              &
         &   self%cosmologyFunctions_%expansionFactor             (             &
         &                                                         basic%time() &
         &                                                        )             &
         &                                                       )
    ! Determine which fitting function to use.
    select case (self%fittingFunction%ID)
    case (klypin2015FittingFunctionEqn24%ID)
       ! Evaluate fitting function parameters.
       concentration0=self%fitParameters%interpolate(redshift,table=1)
       gamma         =self%fitParameters%interpolate(redshift,table=2)
       mass0         =self%fitParameters%interpolate(redshift,table=3)*massReference
       ! Evaluate the concentration.
       klypin2015Concentration=+concentration0                       &
            &                  *(                                    &
            &                    +1.0d0                              &
            &                    +(massLittleH/mass0        )**0.4d0 &
            &                  )                                     &
            &                  /  (massLittleH/massReference)**gamma
    case (klypin2015FittingFunctionEqn25%ID)
       ! Find σ(M).
       sigma         =self%cosmologicalMassVariance_%rootVariance(basic%mass(),self%cosmologyFunctions_%cosmicTime(1.0d0))
       ! Evaluate fitting function parameters.
       a0            =self%fitParameters%interpolate(redshift,table=1)
       b0            =self%fitParameters%interpolate(redshift,table=2)
       ! Evaluate the concentration.
       klypin2015Concentration=+b0          &
            &                  *(           &
            &                    +1.00d0    &
            &                    +7.37d0    &
            &                    *(         &
            &                      +sigma   &
            &                      /a0      &
            &                     )**0.75d0 &
            &                   )           &
            &                  *(           &
            &                    +1.00d0    &
            &                    +0.14d0    &
            &                    /(         &
            &                      +sigma   &
            &                      /a0      &
            &                     )**2      &
            &                   )
    case default
       klypin2015Concentration=0.0d0
       call Error_Report('unknown fit type'//{introspection:location})
    end select
    return
  end function klypin2015Concentration

  function klypin2015DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{klypin_multidark_2014} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass              ), pointer       :: klypin2015DensityContrastDefinition
    class(darkMatterProfileConcentrationKlypin2015), intent(inout) :: self

    klypin2015DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function klypin2015DensityContrastDefinition

  function klypin2015DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{klypin_multidark_2014} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass               ), pointer       :: klypin2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationKlypin2015), intent(inout) :: self

    klypin2015DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function klypin2015DarkMatterProfileDefinition
