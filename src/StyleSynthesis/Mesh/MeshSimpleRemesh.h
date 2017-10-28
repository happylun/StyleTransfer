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

#include <queue>

#include "Data/StyleSynthesisTypes.h"

namespace StyleSynthesis {

	class MeshSimpleRemesh {

	private:

		// make it non-instantiable
		MeshSimpleRemesh() {}
		~MeshSimpleRemesh() {}

	public:

		static bool remesh(TTriangleMesh &mesh);

	public:

		// type defs

		struct TMeshEdge {
			vec2i edgeVertices;
			vec2i neighborVertices;
			float edgeLength;
		};

		typedef map<vec2i, TMeshEdge> TMeshEdgeSet;

	public:

		static bool buildEdgeStructure(
			TTriangleMesh &inMesh,
			TMeshEdgeSet &outMeshEdges);

	private:

		static bool edgeSplit(
			TTriangleMesh &mesh,
			TMeshEdgeSet &edges);

		static bool edgeCollapse(
			TTriangleMesh &mesh,
			TMeshEdgeSet &edges);

		static bool edgeFlip(
			TTriangleMesh &mesh,
			TMeshEdgeSet &edges);

		static bool edgeSmooth(TTriangleMesh &mesh);

		static bool splitNonManifoldEdges(TTriangleMesh &mesh);

		inline static vec2i edgeKey(int a, int b) { return a < b ? vec2i(a, b) : vec2i(b, a); }
		inline static vec2i edgeKey(vec2i e) { return edgeKey(e[0], e[1]); }
		static float edgeFlipBenefit(TTriangleMesh &mesh, TMeshEdge &edge);
		static float faceCosMaxAngle(vec3 p1, vec3 p2, vec3 p3);
	};
}