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
Contains a module which provides various enumerations for the galactic structure functions.
!!}

module Galactic_Structure_Options
  !!{
  Provides various internal option codes for the galactic structure functions.
  !!}
  implicit none
  public

  !![
  <enumeration>
   <name>massType</name>
   <description>Used to specify the mass type(s) to be queried in galactic structure functions.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="all"       description="All mass types"                       />
   <entry label="dark"      description="Dark matter mass"                     />
   <entry label="baryonic"  description="Baryonic mass"                        />
   <entry label="galactic"  description="Galactic mass (i.e. mass in a galaxy)"/>
   <entry label="gaseous"   description="Gaseous mass"                         />
   <entry label="stellar"   description="Stellar mass"                         />
   <entry label="blackHole" description="Mass in black holes"                  />
   <entry label="unknown"   description="Unknown mass type"                    />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>componentType</name>
   <description>Used to specify the component(s) to be queried in galactic structure functions.</description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <entry label="all"                description="All components"                />
   <entry label="disk"               description="The disk component"            />
   <entry label="spheroid"           description="The spheroid component"        />
   <entry label="hotHalo"            description="The hot halo (CGM) component"  />
   <entry label="nuclearStarCluster" description="The \gls{nsc} component"       />
   <entry label="coldHalo"           description="The cold halo (CGM) component" />
   <entry label="darkHalo"           description="The dark matter halo component"/>
   <entry label="blackHole"          description="The black hole component"      />
   <entry label="darkMatterOnly"     description="The dark matter only component"/>
   <entry label="none"               description="No component"                  />
   <entry label="unknown"            description="Unknown components"            />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>coordinateSystem</name>
   <description>Used to specify the coordinate system of the input coordinates in galactic structure functions.</description>
   <entry label="spherical"   />
   <entry label="cylindrical" />
   <entry label="cartesian"   />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>weightBy</name>
   <description>Used to specify by which quantity to weight the results in galactic structure functions.</description>
   <entry label="mass"       />
   <entry label="luminosity" />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>structureErrorCode</name>
   <description>Error codes for galactic structure functions.</description>
   <entry label="success"     description="Successful completion"/>
   <entry label="infinite"    description="Result is ±∞"         />
   <entry label="integration" description="Integration failed"   />
  </enumeration>
  !!]

  ! Null value to use when no weighting is to be applied.
  integer         , parameter, public :: weightIndexNull=0

  ! Suitably large value to represent infinite radius.
  double precision, parameter, public :: radiusLarge    =1.0d10

end module Galactic_Structure_Options
