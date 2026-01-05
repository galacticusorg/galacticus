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
Contains a module which defines a labels for the {\normalfont \ttfamily table} class.
!!}

module Table_Labels
  !!{
  Defines labels for the {\normalfont \ttfamily table} class.
  !!}
  private

  ! Enumeration for extrapolation options in tables.
  !![
  <enumeration>
   <name>extrapolationType</name>
   <description>Used to specify the type of extrapolation to use when interpolating in tables.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="extrapolate" description="Extrapolate beyond the range of the tabulated data"/>
   <entry label="fix"         description="Fix the value to that at the last tabulated point" />
   <entry label="abort"       description="Abort if extrapolation would be required"          />
   <entry label="zero"        description="Return zero beyond the range of the tabulated data"/>
  </enumeration>
  !!]

end module Table_Labels
