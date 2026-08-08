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

!+ This file was generated with Claude and verified by Andrew Benson

!!{RST
Implements a weight operator which applies the observational selection function for Milky Way satellite galaxies.
!!}

  use :: Galactic_Filters         , only : enumerationPositionTypeType               , positionTypeOrbital               , positionTypePosition
  use :: ISO_Varying_String       , only : varying_string
  use :: Node_Property_Extractors , only : nodePropertyExtractorRadiusHalfMassStellar, nodePropertyExtractorRadiusOrbital
  use :: Numerical_Interpolation  , only : interpolator

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorLocalGroupDetection" docformat="rst">
   <description>
   Multiplies the weight of each node by the probability that it would be detected as a Milky Way satellite galaxy in a
   census of the Local Group.

   The probability is read from a tabulation of the observational selection function as a function of galactocentric radius,
   absolute magnitude, and projected half-light radius, and is interpolated trilinearly (linearly in absolute magnitude, and
   in the logarithms of the two radii). The default tabulation is that of the DELVE Milky Way Census I
   :cite:p:`tan_delve_2025`, which combines searches of the Dark Energy Survey Year 6, DECam Local Volume Exploration Survey
   Data Release 3, and Pan-STARRS1 Data Release 1 datasets covering :math:`\sim 91\%` of the sky at high Galactic latitude.

   Since the satellite positions predicted by Galacticus carry no meaningful orientation with respect to the Galactic plane
   or to the survey footprints, the tabulated selection function is averaged over all directions of the satellite from the
   Galactic center---which also accounts for the offset of the Sun from the Galactic center, and for the Zone of Avoidance.
   The tabulation therefore assumes an isotropic satellite distribution, and does not capture any anisotropy such as that
   associated with satellites clustered around the Large Magellanic Cloud. It is built by
   ``scripts/aux/localGroupSelectionFunction.py``.

   The absolute magnitude is computed from the total stellar luminosity of the node in the filter ``[filterName]``, which
   must therefore be included in the set of luminosities computed by the model. The projected half-light radius is taken to
   be ``[radiusHalfLightRatio]`` times the three-dimensional stellar half-mass radius of the node.
   </description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorLocalGroupDetection
     !!{RST
     A weight operator class which applies the observational selection function for Milky Way satellite galaxies.
     !!}
     private
     type            (varying_string                            )                              :: fileName                   , filterName            , &
          &                                                                                       filterType
     double precision                                                                          :: redshift                   , radiusHalfLightRatio
     double precision                                                                          :: radiusMaximum
     double precision                                                                          :: radiusTableMinimum         , radiusTableMaximum    , &
          &                                                                                       magnitudeMinimum           , magnitudeMaximum      , &
          &                                                                                       radiusHalfLightMinimum     , radiusHalfLightMaximum
     integer                                                                                   :: indexLuminosity
     type            (enumerationPositionTypeType               )                              :: positionType
     double precision                                            , allocatable, dimension(:,:,:) :: detectionProbability
     type            (interpolator                              )                              :: interpolatorRadius         , interpolatorMagnitude , &
          &                                                                                       interpolatorRadiusHalfLight
     type            (nodePropertyExtractorRadiusOrbital        )                              :: radiusOrbital_
     type            (nodePropertyExtractorRadiusHalfMassStellar)                              :: radiusHalfMassStellar_
   contains
     procedure :: operate => localGroupDetectionOperate
  end type outputAnalysisWeightOperatorLocalGroupDetection

  interface outputAnalysisWeightOperatorLocalGroupDetection
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisWeightOperatorLocalGroupDetection` output analysis weight operator class.
     !!}
     module procedure localGroupDetectionConstructorParameters
     module procedure localGroupDetectionConstructorInternal
  end interface outputAnalysisWeightOperatorLocalGroupDetection

contains

  function localGroupDetectionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisWeightOperatorLocalGroupDetection` output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Galactic_Filters  , only : enumerationPositionTypeEncode
    use :: Input_Parameters  , only : inputParameter               , inputParameters
    use :: Input_Paths       , only : inputPath                    , pathTypeDataStatic
    use :: ISO_Varying_String, only : var_str                      , operator(//)
    implicit none
    type            (outputAnalysisWeightOperatorLocalGroupDetection)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    type            (varying_string                                 )                :: fileName     , filterName          , &
         &                                                                              filterType   , positionType
    double precision                                                                 :: redshift     , radiusHalfLightRatio, &
         &                                                                              radiusMaximum

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <source>parameters</source>
      <defaultValue>inputPath(pathTypeDataStatic)//'observations/localGroup/localGroupSelectionFunctionDELVE.hdf5'</defaultValue>
      <description>
      The name of the file containing the tabulated detection probability for Milky Way satellite galaxies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>filterName</name>
      <source>parameters</source>
      <defaultValue>var_str('Buser_V')</defaultValue>
      <description>
      The name of the filter in which absolute magnitudes are computed. The detection probability is tabulated as a function of
      absolute :math:`V`-band magnitude, so this should be a :math:`V`-band filter.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>filterType</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>
      The type (``rest`` or ``observed``) of the luminosity used to compute absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The redshift of the luminosity used to compute absolute magnitudes.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusHalfLightRatio</name>
      <source>parameters</source>
      <defaultValue>0.75d0</defaultValue>
      <defaultSource>
      Appropriate for a Plummer profile, for which the projected half-light radius is :math:`0.766` times the three-dimensional
      half-mass radius.
      </defaultSource>
      <description>
      The ratio of the projected, azimuthally-averaged half-light radius to the three-dimensional stellar half-mass radius.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusMaximum</name>
      <source>parameters</source>
      <defaultValue>0.3d0</defaultValue>
      <description>
      The maximum galactocentric radius at which a satellite is considered to belong to the sample. This should match the radius
      used to select the observed sample against which the model is being compared---it defines the sample, and so is imposed in
      addition to the tabulated selection function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('orbital')</defaultValue>
      <description>
      The type of position used to determine the galactocentric radius of the satellite.
      </description>
    </inputParameter>
    !!]
    self=outputAnalysisWeightOperatorLocalGroupDetection(fileName,filterName,filterType,redshift,radiusHalfLightRatio,radiusMaximum,enumerationPositionTypeEncode(positionType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function localGroupDetectionConstructorParameters

  function localGroupDetectionConstructorInternal(fileName,filterName,filterType,redshift,radiusHalfLightRatio,radiusMaximum,positionType) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisWeightOperatorLocalGroupDetection` output analysis weight operator class.
    !!}
    use :: Error                         , only : Error_Report
    use :: HDF5_Access                   , only : hdf5Access
    use :: IO_HDF5                       , only : hdf5File
    use :: ISO_Varying_String            , only : char
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: Table_Labels                  , only : extrapolationTypeFix   , extrapolationTypeZero
    implicit none
    type            (outputAnalysisWeightOperatorLocalGroupDetection)                              :: self
    type            (varying_string                                 ), intent(in   )               :: fileName            , filterName          , &
         &                                                                                            filterType
    double precision                                                 , intent(in   )               :: redshift            , radiusHalfLightRatio, &
         &                                                                                            radiusMaximum
    type            (enumerationPositionTypeType                    ), intent(in   )               :: positionType
    double precision                                                 , allocatable  , dimension(:) :: radius              , magnitude           , &
         &                                                                                            radiusHalfLight
    type            (hdf5File                                       )                              :: selectionFunctionFile
    !![
    <constructorAssign variables="fileName, filterName, filterType, redshift, radiusHalfLightRatio, radiusMaximum, positionType"/>
    !!]

    ! Read the tabulated detection probability.
    !$ call hdf5Access%set()
    selectionFunctionFile=hdf5File(char(fileName),readOnly=.true.)
    call selectionFunctionFile%readDataset('radiusGalactocentric',radius                   )
    call selectionFunctionFile%readDataset('magnitudeAbsoluteV'  ,magnitude                )
    call selectionFunctionFile%readDataset('radiusHalfLight'     ,radiusHalfLight          )
    call selectionFunctionFile%readDataset('detectionProbability',self%detectionProbability)
    !$ call hdf5Access%unset()
    ! Validate the shape of the table. The three axes are of different lengths, so this will detect a table which has been
    ! written with its dimensions transposed.
    if     (                                                                &
         &   size(self%detectionProbability,dim=1) /= size(radius         ) &
         &  .or.                                                            &
         &   size(self%detectionProbability,dim=2) /= size(magnitude      ) &
         &  .or.                                                            &
         &   size(self%detectionProbability,dim=3) /= size(radiusHalfLight) &
         & ) call Error_Report('shape of tabulated detection probability in file "'//fileName//'" does not match its axes'//{introspection:location})
    ! Build interpolators in the tabulated coordinates. Both radii are interpolated logarithmically. Outside of the tabulated
    ! range of galactocentric radius the satellite is undetectable---it is either inside the range excluded by the census, or
    ! beyond the range over which the selection function was determined---so zero extrapolation is used. For the absolute
    ! magnitude and half-light radius the tabulated range spans all plausible galaxies, so the boundary values are used.
    self%interpolatorRadius          =interpolator(log10(radius         ),extrapolationType=extrapolationTypeZero)
    self%interpolatorMagnitude       =interpolator(      magnitude       ,extrapolationType=extrapolationTypeFix )
    self%interpolatorRadiusHalfLight =interpolator(log10(radiusHalfLight),extrapolationType=extrapolationTypeFix )
    self%radiusTableMinimum          =      radius         (                   1 )
    self%radiusTableMaximum          =      radius         (size(radius         ))
    self%magnitudeMinimum            =      magnitude      (                   1 )
    self%magnitudeMaximum            =      magnitude      (size(magnitude      ))
    self%radiusHalfLightMinimum      =      radiusHalfLight(                   1 )
    self%radiusHalfLightMaximum      =      radiusHalfLight(size(radiusHalfLight))
    ! Find the luminosity corresponding to the requested filter.
    self%indexLuminosity=unitStellarLuminosities%index(char(filterName),char(filterType),redshift)
    ! Construct property extractors.
    self%radiusOrbital_        =nodePropertyExtractorRadiusOrbital        ()
    self%radiusHalfMassStellar_=nodePropertyExtractorRadiusHalfMassStellar()
    return
  end function localGroupDetectionConstructorInternal

  double precision function localGroupDetectionOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{RST
    Implement a ``localGroupDetection`` output analysis weight operator.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Error                         , only : Error_Report
    use            :: Galactic_Structure_Options    , only : massTypeStellar       , weightByLuminosity
    use            :: Galacticus_Nodes              , only : nodeComponentBasic    , nodeComponentPosition
    use            :: Mass_Distributions            , only : massDistributionClass
    use            :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use            :: Vectors                       , only : Vector_Magnitude
    implicit none
    class           (outputAnalysisWeightOperatorLocalGroupDetection), intent(inout)  :: self
    type            (treeNode                                       ), intent(inout)  :: node
    double precision                                                 , intent(in   )  :: propertyValue         , propertyValueIntrinsic, &
         &                                                                               weightValue
    type            (enumerationOutputAnalysisPropertyTypeType      ), intent(in   )  :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType  ), intent(in   )  :: propertyQuantity
    integer         (c_size_t                                       ), intent(in   )  :: outputIndex
    class           (nodeComponentBasic                             ), pointer        :: basic
    class           (nodeComponentPosition                          ), pointer        :: position
    class           (massDistributionClass                          ), pointer        :: massDistribution_
    double precision                                                 , dimension(0:1) :: factorsRadius         , factorsMagnitude      , &
         &                                                                               factorsRadiusHalfLight
    integer         (c_size_t                                       )                 :: indexRadius           , indexMagnitude        , &
         &                                                                               indexRadiusHalfLight
    integer                                                                           :: i                     , j                     , &
         &                                                                               k
    double precision                                                                  :: radius                , radiusHalfLight       , &
         &                                                                               magnitude             , luminosity            , &
         &                                                                               probability
    !$GLC attributes unused :: propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity, outputIndex

    ! Galaxies which are not being analyzed at an output time for which this luminosity is available can not be assigned a
    ! detection probability.
    basic => node%basic()
    if (.not.unitStellarLuminosities%isOutput(self%indexLuminosity,basic%time())) then
       localGroupDetectionOperate=0.0d0
       return
    end if
    ! Find the galactocentric radius of the satellite.
    select case (self%positionType%ID)
    case (positionTypePosition%ID)
       position => node    %position()
       radius   =  Vector_Magnitude(position%position())
    case (positionTypeOrbital %ID)
       radius   =  self%radiusOrbital_%extract(node)
    case default
       radius   =  0.0d0
       call Error_Report('unknown position type'//{introspection:location})
    end select
    ! Satellites beyond the radius which defines the sample, or outside of the tabulated range of galactocentric radius (where
    ! they are undetectable), receive zero weight.
    if (radius > self%radiusMaximum .or. radius < self%radiusTableMinimum .or. radius > self%radiusTableMaximum) then
       localGroupDetectionOperate=0.0d0
       return
    end if
    ! Find the absolute magnitude of the satellite. Galaxies with no stellar luminosity are undetectable.
    massDistribution_ => node             %massDistribution(massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=self%indexLuminosity)
    luminosity        =  massDistribution_%massTotal       (                                                                                     )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (luminosity <= 0.0d0) then
       localGroupDetectionOperate=0.0d0
       return
    end if
    magnitude=-2.5d0*log10(luminosity)
    ! Find the projected half-light radius of the satellite.
    radiusHalfLight=+self%radiusHalfLightRatio                 &
         &          *self%radiusHalfMassStellar_%extract(node)
    ! Clamp the magnitude and half-light radius to the tabulated ranges. The tabulation spans all plausible galaxies, so this
    ! affects only pathological cases (for example, a galaxy with an unresolved size).
    magnitude      =min(max(magnitude      ,self%magnitudeMinimum      ),self%magnitudeMaximum      )
    radiusHalfLight=min(max(radiusHalfLight,self%radiusHalfLightMinimum),self%radiusHalfLightMaximum)
    ! Interpolate in the tabulated detection probability.
    call self%interpolatorRadius         %linearFactors(log10(radius         ),indexRadius         ,factorsRadius         )
    call self%interpolatorMagnitude      %linearFactors(      magnitude       ,indexMagnitude      ,factorsMagnitude      )
    call self%interpolatorRadiusHalfLight%linearFactors(log10(radiusHalfLight),indexRadiusHalfLight,factorsRadiusHalfLight)
    probability=0.0d0
    do i=0,1
       do j=0,1
          do k=0,1
             probability=+probability                                                                      &
                  &      +factorsRadius            (            i                                        ) &
                  &      *factorsMagnitude         (                             j                       ) &
                  &      *factorsRadiusHalfLight   (                                                    k) &
                  &      *self%detectionProbability(indexRadius+i,indexMagnitude+j,indexRadiusHalfLight+k)
          end do
       end do
    end do
    localGroupDetectionOperate=+weightValue                       &
         &                     *min(max(probability,0.0d0),1.0d0)
    return
  end function localGroupDetectionOperate
