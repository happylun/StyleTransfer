//============================================================================
//
// This file is part of the Style Transfer project.
//
// Copyright (c) 2016
// -Zhaoliang Lun, Evangelos Kalogerakis (authors of the code) / UMass-Amherst
//
// This is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software.  If not, see <http://www.gnu.org/licenses/>.
//
//============================================================================

#pragma once

/*
	Modified (to make it compatible with VS2012):
		libigl/include/igl/sort.cpp#L177&L178
			swap(a,b);
			swap(ai, bi);
					=>
			std::swap(a,b);
			std::swap(ai, bi);
		libigl/include/igl/unique_simplices.cpp#L40
			for(size_t i = 0;i<mff;i++)
					=>
			for(int i = 0;i<mff;i++)
		libigl/include/igl/cgal/SelfIntersectMesh.h#L493&L700
			swap(i,j);
					=>
			std::swap(i,j);
	Modified (to prevent numerical overflow)
		libigl/include/igl/min_quad_with_fixed.cpp#L99
			assert(Auu.size() > 0 && "There should be at least one unknown.");
					=>
			assert(Auu.rows() > 0 && "There should be at least one unknown.");
	Modified (seems buggy)
		libigl/include/igl/cotmatrix_entries.cpp#L69
			volume(l,vol);
					=>
			//volume(l,vol);
			volume(V,F,vol);
		
*/

#pragma warning(push)
#pragma warning(disable:4018)
#pragma warning(disable:4101)
#pragma warning(disable:4129)
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4305)
#pragma warning(disable:4503)
#pragma warning(disable:4667)
#pragma warning(disable:4800)
#pragma warning(disable:4819)
#pragma warning(disable:4996)