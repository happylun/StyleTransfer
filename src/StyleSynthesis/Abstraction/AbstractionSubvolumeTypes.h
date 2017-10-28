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

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	namespace AbstractionSubvolumeNamespace {

		struct TPrimitive {

			static const int TYPE_REC_PRISM  = 0;
			static const int TYPE_TRI_PRISM  = 1;
			static const int TYPE_CYLINDER   = 2;
			static const int TYPE_TRUNC_CONE = 3;
			static const int NUM_TYPES       = 4;

			TPrimitive(int type) {
				mArea = 0;
				mType = type;
			}

			virtual bool computeArea() = 0;
			virtual bool fit(TPointSet *points, vec4i samples) = 0;
			virtual bool refit(TPointSet &points, int index) = 0;
			virtual float distance(vec3 point) = 0;
			virtual bool inside(vec3 point) = 0;
			virtual bool tesselate(TTriangleMesh &mesh) = 0;

			// common properties
			float mArea;
			int mType;
		};

		struct TRecPrism : TPrimitive {

			TRecPrism() : TPrimitive(TYPE_REC_PRISM) {}

			bool computeArea();
			bool fit(TPointSet *points, vec4i samples);
			bool refit(TPointSet &points, int index);
			float distance(vec3 point);
			bool inside(vec3 point);
			bool tesselate(TTriangleMesh &mesh);

			vec3 mCorner[4];
		};

		struct TTriPrism : TPrimitive {

			TTriPrism() : TPrimitive(TYPE_TRI_PRISM) {}

			bool computeArea();
			bool fit(TPointSet *points, vec4i samples);
			bool refit(TPointSet &points, int index);
			float distance(vec3 point);
			bool inside(vec3 point);
			bool tesselate(TTriangleMesh &mesh);

			vec3 mCorner[4];
		};

		struct TCylinder : TPrimitive {

			TCylinder() : TPrimitive(TYPE_CYLINDER) {}

			bool computeArea();
			bool fit(TPointSet *points, vec4i samples);
			bool refit(TPointSet &points, int index);
			float distance(vec3 point);
			bool inside(vec3 point);
			bool tesselate(TTriangleMesh &mesh);

			vec3 mCenter[2];
			float mRadius;
		};

		struct TTruncCone : TPrimitive {

			TTruncCone() : TPrimitive(TYPE_TRUNC_CONE) {}

			bool computeArea();
			bool fit(TPointSet *points, vec4i samples);
			bool refit(TPointSet &points, int index);
			float distance(vec3 point);
			bool inside(vec3 point);
			bool tesselate(TTriangleMesh &mesh);

			vec3 mCenter[3];
			float mAngle;
		};
	}
}
