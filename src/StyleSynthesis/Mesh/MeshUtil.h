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

#include <Eigen/Dense>

#include "Data/StyleSynthesisTypes.h"

#include "Library/TheaKDTreeHelper.h"

namespace StyleSynthesis {

	class MeshUtil {

	private:

		// make it non-instantiable
		MeshUtil() {}
		~MeshUtil() {}

	public:

		static bool saveMesh(string fileName, TTriangleMesh &mesh, bool ascii = false);
		static bool loadMesh(string fileName, TTriangleMesh &mesh);

		static bool mesh2mat(TTriangleMesh &mesh, Eigen::MatrixXd &matV, Eigen::MatrixXi &matF);
		static bool mat2mesh(TTriangleMesh &mesh, Eigen::MatrixXd &matV, Eigen::MatrixXi &matF);
		static bool tet2mat(TTetrahedralMesh &tet, Eigen::MatrixXd &matV, Eigen::MatrixXi &matT);
		static bool mat2tet(TTetrahedralMesh &tet, Eigen::MatrixXd &matV, Eigen::MatrixXi &matT);
		static bool tet2mesh(TTetrahedralMesh &tet, TTriangleMesh &mesh);
		static bool tet2meshSimple(TTetrahedralMesh &tet, TTriangleMesh &mesh);

		static bool removeDuplicateVertices(
			vector<vec3> &inVertices,
			vector<vec3> &outVertices,
			vector<int> &outVertexIndices, // vertex ID in outVertices : # of vertices in inVertices
			double eps = -1);

		static bool removeDuplicateVertices(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outVertexIndices, // vertex ID in outMesh : # of vertices in inMesh
			double eps = -1);

		static bool removeDuplicateFaces(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices, // face ID in outMesh : # of faces in inMesh
			double eps = -1);

		static bool removeDegeneratedFaces(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices, // face ID in outMesh : # of faces in inMesh
			double eps = -1);

		static bool cleanUp(TTriangleMesh &mesh);
		static bool cleanUp(TTetrahedralMesh &mesh);

		static bool buildVertexGraph(TTriangleMesh &mesh, vector<vector<int>> &graph);
		static bool buildFaceGraph(TTriangleMesh &mesh, vector<vector<int>> &graph);

		static bool subdivideMesh(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices,
			double radius = 0);

		static bool subdivideMeshKeepTopology(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices,
			double radius = 0);

		static bool subdivideMeshMidPoint(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices,
			double radius = 0);

		static bool splitMeshIntoConnectedComponents(
			TTriangleMesh &inMesh,
			vector<TTriangleMesh> &outMesh,
			vector<vector<int>> &outFaceIndices);

		static bool splitMeshIntoCleanedConnectedComponents(
			TTriangleMesh &inMesh,
			vector<TTriangleMesh> &outMesh);

		static bool splitMeshAlongCrease(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> *outLabels = 0);

		static bool extractSubMesh(
			TTriangleMesh &inMesh,
			vector<int> &subMeshFaces,
			TTriangleMesh &outMesh);

		static bool extractMeshFromSamplePoints(
			TTriangleMesh &inMesh,
			TSampleSet &inSamples,
			vector<bool> &inFlags,
			TTriangleMesh &outMesh);

		static bool weldMeshFaces(
			TTriangleMesh &inMesh,
			TTriangleMesh &outMesh,
			vector<int> &outFaceIndices, // face ID in outMesh : # of faces in inMesh
			double eps = -1);

		static bool reorientComponentFaces(TTriangleMesh &inoutMesh);
		static bool reorientMeshFaces(TTriangleMesh &inoutMesh);
		static bool reorientAnyMesh(TTriangleMesh &inoutMesh);
		static bool reorientPatch(
			TKDTree &inMeshTree,
			TTriangleMesh &inMesh,
			vector<int> &inPatch,
			bool &outFlipFlag,
			double eps = -1);

		static bool recomputeNormals(TTriangleMesh &mesh);
		static bool computeDihedralAngle(vec3 center1, vec3 normal1, vec3 center2, vec3 normal2, double &angle);
		static bool computeAABB(TTriangleMesh &mesh, vec3 &bbMin, vec3 &bbMax);
		static bool computeFaceArea(TTriangleMesh &mesh, double &area);
		static bool computeVolume(TTriangleMesh &mesh, double &volume);
		
		static bool buildKdTree(TTriangleMesh &mesh, TKDTree &tree, TKDTreeData &treeData); // kd tree of face triangles
		static bool buildKdTree(TTriangleMesh &mesh, SKDTree &tree, SKDTreeData &treeData); // kd tree of face center points

		static bool distanceToMesh(TTriangleMesh &mesh, SKDTree &faceTree, vec3 &point, double &distance);
		static bool checkInsideMesh(TKDTree &tree, vec3 &point, bool &insideFlag, double eps = 1e-5);
		
	};
}