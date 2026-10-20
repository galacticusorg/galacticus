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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which reads the SKIRTOR library of models of the dusty torus around an AGN.
!!}

module Dust_AGN_Torus_SKIRTOR
  !!{RST
  Reads the SKIRTOR library of radiative transfer models of the dusty torus around an AGN
  :cite:p:`stalevski_3d_2012,stalevski_dust_2016`, as converted to HDF5 by ``dataPreparation/convertSKIRTOR.py`` of the
  ``galacticusDevTools`` repository (https://github.com/galacticusorg/galacticusDevTools), and provides the weights for
  averaging its tabulations over orientation.

  The library tabulates models on a grid of edge-on optical depth at :math:`9.7\,\mu\hbox{m}` (``opticalDepth9.7``),
  exponents of the radial (``exponentRadial``, :math:`p`) and polar (``exponentPolar``, :math:`q`) gradients of dust
  density, ratio of outer to inner radius (``ratioRadii``), and angle between the equatorial plane and the edge of the
  torus (``angleTorus``), each seen at ten inclinations from face-on to edge-on. The first four select a set of models,
  which differ only in ``angleTorus``; :galacticus-module:`Dust_AGN_Torus_SKIRTOR` returns that set, ordered by
  ``angleTorus``.
  !!}
  private
  public :: skirtorRead, skirtorWeightsInclination

contains

  subroutine skirtorRead(fileName,url,depthOptical,exponentRadial,exponentPolar,ratioRadii,angleTorus,wavelength,inclination,transmissionDirect,transmissionScattered,emissionDustMean)
    !!{RST
    Read the set of SKIRTOR models with the given ``depthOptical``, ``exponentRadial``, ``exponentPolar`` and
    ``ratioRadii``, ordered by the angle between the equatorial plane and the edge of the torus, returned in
    ``angleTorus`` (in radians). The tabulation is found in the ``dust/agnTorus`` directory of the dynamic datasets path,
    and downloaded from ``url`` if it is not there. Wavelengths are returned in Å, and inclinations in radians. The
    transmissions are indexed (wavelength, inclination, model), and the orientation-averaged dust emission (wavelength,
    model).
    !!}
    use            :: Error                           , only : Error_Report
    use            :: File_Utilities                  , only : Directory_Make    , File_Exists        , File_Lock     , File_Remove , &
         &                                                     File_Unlock       , lockDescriptor
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: Input_Paths                     , only : inputPath         , pathTypeDataDynamic
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: IO_HDF5                         , only : hdf5File
    use            :: ISO_Varying_String              , only : assignment(=)     , char               , operator(//)  , operator(/=), &
         &                                                     operator(==)      , var_str            , varying_string
    use            :: Numerical_Constants_Astronomical, only : degreesToRadians
    use            :: Numerical_Constants_Units       , only : micronsToAngstroms
    use            :: System_Download                 , only : download
    implicit none
    type            (varying_string), intent(in   )                                          :: fileName                    , url
    double precision                , intent(in   )                                          :: depthOptical                , exponentRadial          , &
         &                                                                                      exponentPolar               , ratioRadii
    double precision                , intent(  out), allocatable, dimension(:    )           :: angleTorus                  , wavelength              , &
         &                                                                                      inclination
    double precision                , intent(  out), allocatable, dimension(:,:,:), optional :: transmissionDirect          , transmissionScattered
    double precision                , intent(  out), allocatable, dimension(:,:  ), optional :: emissionDustMean
    ! Tolerance on the match between a requested parameter value and a tabulated one.
    double precision                , parameter                                              :: toleranceParameter   =1.0d-6
    double precision                               , allocatable, dimension(:    )           :: depthOpticalModel           , exponentRadialModel     , &
         &                                                                                      exponentPolarModel          , ratioRadiiModel         , &
         &                                                                                      angleTorusModel
    double precision                               , allocatable, dimension(:,:,:)           :: transmissionDirectAll       , transmissionScatteredAll
    double precision                               , allocatable, dimension(:,:  )           :: emissionDustMeanAll
    integer                                        , allocatable, dimension(:    )           :: indexModel
    logical                                        , allocatable, dimension(:    )           :: selected
    type            (lockDescriptor)                                                         :: lock
    type            (varying_string)                                                         :: pathFile                    , pathDirectory
    integer         (c_size_t      )                                                         :: sizeFile
    integer                                                                                  :: status                      , i                       , &
         &                                                                                      j                           , indexSwap

    ! Locate the tabulation, downloading it if we do not already have it. The lock makes concurrent processes wait for
    ! the first of them to finish the download rather than each starting one of their own.
    pathDirectory=inputPath(pathTypeDataDynamic)//'dust/agnTorus'
    pathFile     =pathDirectory//'/'//fileName
    if (.not.File_Exists(pathFile)) then
       if (url == 'none') call Error_Report(                                                                                                                      &
            &                               var_str('the SKIRTOR tabulation `')//pathFile//'` is not present, and no `url` was given from which to download it'// &
            &                               ' - it can be made from the SKIRTOR library with `dataPreparation/convertSKIRTOR.py` of galacticusDevTools'        // &
            &                               {introspection:location}                                                                                              &
            &                              )
       call Directory_Make(pathDirectory                          )
       call File_Lock     (pathFile     ,lock,lockIsShared=.false.)
       if (.not.File_Exists(pathFile)) then
          call download(char(url),char(pathFile),retries=3,retryWait=10,status=status)
          ! Judge the download by what actually landed on disk, not by the reported status alone.
          sizeFile=0_c_size_t
          if (File_Exists(pathFile)) inquire(file=char(pathFile),size=sizeFile)
          if (status /= 0 .or. sizeFile <= 0_c_size_t) then
             if (File_Exists(pathFile)) call File_Remove(pathFile)
             call File_Unlock(lock)
             call Error_Report('unable to download the SKIRTOR tabulation from `'//url//'`'//{introspection:location})
          end if
       end if
       call File_Unlock(lock)
    end if
    !$ call hdf5Access%set()
    hdf5ReadScope : block
      type(hdf5File) :: file
      file=hdf5File(char(pathFile),readOnly=.true.)
      call file%readDataset('wavelength'     ,wavelength         )
      call file%readDataset('inclination'    ,inclination        )
      call file%readDataset('opticalDepth9.7',depthOpticalModel  )
      call file%readDataset('exponentRadial' ,exponentRadialModel)
      call file%readDataset('exponentPolar'  ,exponentPolarModel )
      call file%readDataset('ratioRadii'     ,ratioRadiiModel    )
      call file%readDataset('angleTorus'     ,angleTorusModel    )
      if (present(transmissionDirect   )) call file%readDataset('transmissionDirect'   ,transmissionDirectAll   )
      if (present(transmissionScattered)) call file%readDataset('transmissionScattered',transmissionScatteredAll)
      if (present(emissionDustMean     )) call file%readDataset('emissionDustMean'     ,emissionDustMeanAll     )
    end block hdf5ReadScope
    !$ call hdf5Access%unset()
    wavelength =wavelength *micronsToAngstroms
    inclination=inclination*degreesToRadians
    ! Select the models with the requested parameters.
    selected= abs(depthOpticalModel  -depthOptical  ) <= toleranceParameter*max(1.0d0,abs(depthOptical  )) &
         &   .and.                                                                                         &
         &    abs(exponentRadialModel-exponentRadial) <= toleranceParameter*max(1.0d0,abs(exponentRadial)) &
         &   .and.                                                                                         &
         &    abs(exponentPolarModel -exponentPolar ) <= toleranceParameter*max(1.0d0,abs(exponentPolar )) &
         &   .and.                                                                                         &
         &    abs(ratioRadiiModel    -ratioRadii    ) <= toleranceParameter*max(1.0d0,abs(ratioRadii    ))
    if (count(selected) < 1)                                                                                                       &
         & call Error_Report(                                                                                                      &
         &                   var_str('no SKIRTOR model has the requested parameters - tabulated values are: `opticalDepth9.7` ')// &
         &                   valuesList(depthOpticalModel  )//'; `exponentRadial` '                                             // &
         &                   valuesList(exponentRadialModel)//'; `exponentPolar` '                                              // &
         &                   valuesList(exponentPolarModel )//'; `ratioRadii` '                                                 // &
         &                   valuesList(ratioRadiiModel    )                                                                    // &
         &                   {introspection:location}                                                                              &
         &                  )
    indexModel=pack([(i,i=1,size(selected))],selected)
    ! Order the models by the angle of the torus.
    do i=2,size(indexModel)
       j=i
       do while (j > 1)
          if (angleTorusModel(indexModel(j-1)) <= angleTorusModel(indexModel(j))) exit
          indexSwap      =indexModel(j-1)
          indexModel(j-1)=indexModel(j  )
          indexModel(j  )=indexSwap
          j              =           j-1
       end do
    end do
    angleTorus=angleTorusModel(indexModel)*degreesToRadians
    if (size(angleTorus) > 1) then
       if (any(angleTorus(2:size(angleTorus)) <= angleTorus(1:size(angleTorus)-1))) &
            & call Error_Report('the SKIRTOR tabulation has more than one model with the same parameters'//{introspection:location})
    end if
    if (present(transmissionDirect   )) transmissionDirect   =transmissionDirectAll   (:,:,indexModel)
    if (present(transmissionScattered)) transmissionScattered=transmissionScatteredAll(:,:,indexModel)
    if (present(emissionDustMean     )) emissionDustMean     =emissionDustMeanAll     (:  ,indexModel)
    return

  contains

    function valuesList(values) result(list)
      !!{RST
      Return a list of the distinct ``values``, in increasing order.
      !!}
      implicit none
      type            (varying_string)                              :: list
      double precision                , intent(in   ), dimension(:) :: values
      double precision                                              :: valueLast, valueNext
      character       (len=16        )                              :: label

      list     =''
      valueLast=-huge(1.0d0)
      do while (any(values > valueLast))
         valueNext=minval(values,mask=values > valueLast)
         write (label,'(g0.4)') valueNext
         if (list /= '') list=list//', '
         list     =list//trim(adjustl(label))
         valueLast=valueNext
      end do
      return
    end function valuesList

  end subroutine skirtorRead

  function skirtorWeightsInclination(inclination) result(weights)
    !!{RST
    Return weights, :math:`w_k`, such that :math:`\sum_k w_k T(i_k)` is the average over orientation,
    :math:`\int_0^{\pi/2} T(i) \sin i\,\mathrm{d}i`, of a quantity :math:`T` tabulated at ``inclinations`` :math:`i_k`
    (in radians, spanning :math:`0` to :math:`\pi/2`) and interpolated linearly in inclination between them. Each weight
    is the integral of the interpolating basis function of its inclination, which is found analytically.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:                ) :: inclination
    double precision               , dimension(size(inclination)) :: weights
    double precision                                              :: inclinationLower, inclinationUpper, &
         &                                                           width
    integer                                                       :: k

    weights=0.0d0
    do k=1,size(inclination)-1
       inclinationLower=inclination(k  )
       inclinationUpper=inclination(k+1)
       width           =inclinationUpper-inclinationLower
       weights(k  )    =+weights(k  )                  &
            &           +(                             &
            &             +width*cos(inclinationLower) &
            &             +      sin(inclinationLower) &
            &             -      sin(inclinationUpper) &
            &            )                             &
            &           /width
       weights(k+1)    =+weights(k+1)                  &
            &           +(                             &
            &             -width*cos(inclinationUpper) &
            &             +      sin(inclinationUpper) &
            &             -      sin(inclinationLower) &
            &            )                             &
            &           /width
    end do
    return
  end function skirtorWeightsInclination

end module Dust_AGN_Torus_SKIRTOR
