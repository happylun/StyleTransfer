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

#include "MeshUtil.h"

#include <fstream>
#include <sstream>
#include <set>
#include <unordered_set>

#include "Data/StyleSynthesisConfig.h"
#include "Utility/PlyExporter.h"
#include "Utility/PlyLoader.h"

#include "Sample/SampleSimplePoissonDisk.h"
#include "Sample/SampleUtil.h"

using namespace StyleSynthesis;

bool MeshUtil::saveMesh(string fileName, TTriangleMesh &mesh, bool ascii) {

	PlyExporter pe;
	if (!pe.addMesh(&mesh.indices, &mesh.positions, mesh.normals.empty() ? 0 : &mesh.normals)) return false;
	if (!pe.output(fileName, ascii)) return false;

	return true;
}

bool MeshUtil::loadMesh(string fileName, TTriangleMesh &mesh) {

	string ext = fileName.substr(fileName.find_last_of(".") + 1);

	if (ext == "ply") {
		if (!PlyLoader::loadMesh(fileName, &mesh.indices, &mesh.positions, &mesh.normals)) return false;
		mesh.amount = (int)mesh.positions.size();
	} else {
		cout << "Unable to load mesh with extension \"" << ext << "\"" << endl;
		return false;
	}

	return true;
}

bool MeshUtil::mesh2mat(TTriangleMesh &mesh, Eigen::MatrixXd &matV, Eigen::MatrixXi &matF) {

	matV.resize(mesh.amount, 3);
	for (int vertID = 0; vertID < mesh.amount; vertID++) {
		vec3 &v = mesh.positions[vertID];
		matV.row(vertID) = Eigen::RowVector3d((double)v[0], (double)v[1], (double)v[2]);
	}

	matF.resize(mesh.indices.size(), 3);
	for (int faceID = 0; faceID < (int)mesh.indices.size(); faceID++) {
		vec3i &f = mesh.indices[faceID];
		matF.row(faceID) = Eigen::RowVector3i(f[0], f[1], f[2]);
	}

	return true;
}

bool MeshUtil::mat2mesh(TTriangleMesh &mesh, Eigen::MatrixXd &matV, Eigen::MatrixXi &matF) {

	mesh.positions.resize((int)matV.rows());
	mesh.amount = (int)mesh.positions.size();

	for (int vertID = 0; vertID < mesh.amount; vertID++) {
		mesh.positions[vertID] = vec3d(matV(vertID, 0), matV(vertID, 1), matV(vertID, 2));
	}

	mesh.indices.resize((int)matF.rows());
	for (int faceID = 0; faceID < (int)mesh.indices.size(); faceID++) {
		mesh.indices[faceID] = vec3i(matF(faceID, 0), matF(faceID, 1), matF(faceID, 2));
	}

	if (!recomputeNormals(mesh)) return false;

	return true;
}

bool MeshUtil::tet2mat(TTetrahedralMesh &tet, Eigen::MatrixXd &matV, Eigen::MatrixXi &matT) {

	matV.resize(tet.amount, 3);
	for (int vertID = 0; vertID < tet.amount; vertID++) {
		vec3 &v = tet.positions[vertID];
		matV.row(vertID) = Eigen::RowVector3d((double)v[0], (double)v[1], (double)v[2]);
	}

	matT.resize(tet.indices.size(), 4);
	for (int tetID = 0; tetID < (int)tet.indices.size(); tetID++) {
		vec4i &t = tet.indices[tetID];
		matT.row(tetID) = Eigen::RowVector4i(t[0], t[1], t[2], t[3]);
	}

	return true;
}

bool MeshUtil::mat2tet(TTetrahedralMesh &tet, Eigen::MatrixXd &matV, Eigen::MatrixXi &matT) {

	tet.positions.resize((int)matV.rows());
	tet.normals.clear();
	tet.amount = (int)tet.positions.size();

	for (int vertID = 0; vertID < tet.amount; vertID++) {
		tet.positions[vertID] = vec3d(matV(vertID, 0), matV(vertID, 1), matV(vertID, 2));
	}

	tet.indices.resize((int)matT.rows());
	for (int faceID = 0; faceID < (int)tet.indices.size(); faceID++) {
		tet.indices[faceID] = vec4i(matT(faceID, 0), matT(faceID, 1), matT(faceID, 2), matT(faceID, 3));
	}

	return true;
}

bool MeshUtil::tet2mesh(TTetrahedralMesh &tet, TTriangleMesh &mesh) {

	int numVerts = tet.amount;
	int numTets = (int)tet.indices.size();

	mesh.positions = tet.positions;
	mesh.normals = tet.normals;
	mesh.amount = numVerts;

	set<vec3i> triangleSet;
	for (int tetID = 0; tetID < numTets; tetID++) {
		auto t = tet.indices[tetID];
		for (int k = 0; k < 4; k++) {
			vec3i key(t[k], t[(k + 1) % 4], t[(k + 2) % 4]);
			if (key[0] > key[1]) swap(key[0], key[1]);
			if (key[0] > key[2]) swap(key[0], key[2]);
			if (key[1] > key[2]) swap(key[1], key[2]);
			triangleSet.insert(key);
		}
	}
	mesh.indices.clear();
	for (vec3i key : triangleSet) {
		mesh.indices.push_back(key);
	}

	return true;
}

bool MeshUtil::tet2meshSimple(TTetrahedralMesh &tet, TTriangleMesh &mesh) {

	int numVerts = tet.amount;
	int numTets = (int)tet.indices.size();

	mesh.positions = tet.positions;
	mesh.normals = tet.normals;
	mesh.amount = numVerts;

	mesh.indices.resize(numTets * 4);
	for (int tetID = 0; tetID < numTets; tetID++) {
		auto t = tet.indices[tetID];
		mesh.indices[tetID * 4 + 0] = vec3i(t[0], t[1], t[3]);
		mesh.indices[tetID * 4 + 1] = vec3i(t[0], t[2], t[1]);
		mesh.indices[tetID * 4 + 2] = vec3i(t[3], t[2], t[0]);
		mesh.indices[tetID * 4 + 3] = vec3i(t[1], t[2], t[3]);
	}

	return true;
}

bool MeshUtil::removeDuplicateVertices(
	vector<vec3> &inVertices,
	vector<vec3> &outVertices,
	vector<int> &outVertexIndices,
	double eps)
{
	// compute eps

	if (eps < 0) {
		vec3 bbMin(FLT_MAX, FLT_MAX, FLT_MAX);
		vec3 bbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		for (vec3 &v : inVertices) {
			bbMin.minimize(v);
			bbMax.maximize(v);
		}
		float bbLen = (bbMax - bbMin).length();
		eps = bbLen * 1e-5; // UNDONE: param eps
	}

	// merge duplicated vertices

	SKDTree vertexTree;
	SKDTreeData vertexTreeData;
	if (!SampleUtil::buildKdTree(inVertices, vertexTree, vertexTreeData)) return false;

	outVertexIndices.resize(inVertices.size());
#pragma omp parallel for
	for (int vertID = 0; vertID < (int)inVertices.size(); vertID++) {
		vec3 p = inVertices[vertID];
		SKDT::NamedPoint queryPoint(p[0], p[1], p[2]);
		Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(100);
		vertexTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult, eps);
		int mapID = vertID; // by default map to itself
		for (int qID = 0; qID < queryResult.size(); qID++) {
			int nbID = (int)vertexTree.getElements()[queryResult[qID].getIndex()].id;
			if (nbID < mapID) mapID = nbID;
		}
		outVertexIndices[vertID] = mapID;
	}

	vector<vec3> newVertices;
	for (int vertID = 0; vertID < (int)inVertices.size(); vertID++) {
		if (outVertexIndices[vertID] == vertID) {
			outVertexIndices[vertID] = (int)newVertices.size();
			newVertices.push_back(inVertices[vertID]);
		} else {
			outVertexIndices[vertID] = outVertexIndices[outVertexIndices[vertID]];
		}
	}

	outVertices.swap(newVertices);

	return true;
}

bool MeshUtil::removeDuplicateVertices(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outVertexIndices,
	double eps)
{
	// compute eps

	if (eps < 0) {
		vec3 bbMin, bbMax;
		if (!computeAABB(inMesh, bbMin, bbMax)) return false;
		float bbLen = (bbMax - bbMin).length();
		eps = bbLen * 1e-5; // UNDONE: param eps
	}

	// merge duplicated vertices

	vector<vec3> newPositions;
	vector<int> newVertexIndices;
	if (!removeDuplicateVertices(inMesh.positions, newPositions, newVertexIndices, eps)) return false;

	// update vertex normals

	vector<vec3> newNormals(newPositions.size());
	for (int vertID = 0; vertID < (int)inMesh.amount; vertID++) {
		newNormals[newVertexIndices[vertID]] = inMesh.normals[vertID];
	}

	// update face indices

	vector<vec3i> newIndices = inMesh.indices;
	for (int faceID = 0; faceID < (int)newIndices.size(); faceID++) {
		vec3i &faceIdx = newIndices[faceID];
		for (int k = 0; k < 3; k++) faceIdx[k] = newVertexIndices[faceIdx[k]];
	}

	outMesh.positions.swap(newPositions);
	outMesh.normals.swap(newNormals);
	outMesh.indices.swap(newIndices);
	outMesh.amount = (int)outMesh.positions.size();

	return true;
}

bool MeshUtil::removeDuplicateFaces(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double eps)
{
	// compute eps

	if (eps < 0) {
		vec3 bbMin, bbMax;
		if (!computeAABB(inMesh, bbMin, bbMax)) return false;
		float bbLen = (bbMax - bbMin).length();
		eps = bbLen * 1e-5; // UNDONE: param eps
	}

	// find unique faces

	vector<vec3i> newFaces(0);
	map<vec3i, int> faceMap; // face point indices => new face ID
	outFaceIndices.resize(inMesh.indices.size());
	for (int faceID = 0; faceID < (int)inMesh.indices.size(); faceID++) {
		vec3i key = inMesh.indices[faceID];
		if (key[0] > key[1]) swap(key[0], key[1]);
		if (key[0] > key[2]) swap(key[0], key[2]);
		if (key[1] > key[2]) swap(key[1], key[2]);
		auto it = faceMap.find(key);
		if (it == faceMap.end()) {
			int newFaceID = (int)newFaces.size();
			outFaceIndices[faceID] = newFaceID;
			faceMap[key] = newFaceID;
			newFaces.push_back(inMesh.indices[faceID]);
		} else {
			outFaceIndices[faceID] = it->second;
		}
	}

	outMesh = inMesh;
	outMesh.indices.swap(newFaces);

	return true;
}

bool MeshUtil::removeDegeneratedFaces(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double eps)
{
	// assumes no duplicated vertices

	// compute eps

	if (eps < 0) {
		vec3 bbMin, bbMax;
		if (!computeAABB(inMesh, bbMin, bbMax)) return false;
		float bbLen = (bbMax - bbMin).length();
		eps = cml::sqr(bbLen * 1e-5); // UNDONE: param eps
	}

	// find degenerated faces & edges to flip

	int numFaces = (int)inMesh.indices.size();
	vector<bool> faceFlags(numFaces, true);
	map<vec2i, int> flipEdgeSet; // flip edge key => inner point ID
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = inMesh.indices[faceID];
		vec3d faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = vec3d(inMesh.positions[faceIdx[k]]);
		double dblArea = cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]).length();
		if (dblArea < eps) {
			faceFlags[faceID] = false;
			double flipEdgeLenSq = 0;
			int flipEdgeID = -1;
			for (int k = 0; k < 3; k++) {
				double edgeLenSq = (faceP[(k + 1) % 3] - faceP[k]).length_squared();
				if (edgeLenSq > flipEdgeLenSq) {
					flipEdgeLenSq = edgeLenSq;
					flipEdgeID = k;
				}
			}
			vec2i flipEdgeKey(faceIdx[flipEdgeID], faceIdx[(flipEdgeID + 1) % 3]);
			if (flipEdgeKey[0] > flipEdgeKey[1]) swap(flipEdgeKey[0], flipEdgeKey[1]);
			int flipPointID = faceIdx[(flipEdgeID + 2) % 3];
			flipEdgeSet[flipEdgeKey] = flipPointID;
		}
	}

	// flip edges to prevent T-vertices

	vector<vec3i> extraFaces;
	map<int, int> extraFaceIndices;
	for (int faceID = 0; faceID < numFaces; faceID++) {
		if (!faceFlags[faceID]) continue; // skip degenerated faces
		vec3i faceIdx = inMesh.indices[faceID];
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
			auto it = flipEdgeSet.find(edgeKey);
			if (it != flipEdgeSet.end()) {
				faceFlags[faceID] = false;
				extraFaceIndices[faceID] = (int)extraFaces.size();
				extraFaces.push_back(vec3i(faceIdx[k], it->second, faceIdx[(k + 2) % 3]));
				extraFaces.push_back(vec3i(it->second, faceIdx[(k + 1) % 3], faceIdx[(k + 2) % 3]));
				break;
			}
		}
	}

	// update mesh faces

	vector<vec3i> newFaces(0);
	outFaceIndices.resize(numFaces);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		if (faceFlags[faceID]) {
			outFaceIndices[faceID] = (int)newFaces.size();
			newFaces.push_back(inMesh.indices[faceID]);
		} else {
			outFaceIndices[faceID] = -1;
		}
	}
	for (auto &it : extraFaceIndices) {
		outFaceIndices[it.first] = it.second + (int)newFaces.size();
	}
	newFaces.insert(newFaces.end(), extraFaces.begin(), extraFaces.end());

	outMesh = inMesh;
	outMesh.indices.swap(newFaces);

	return true;
}

bool MeshUtil::cleanUp(TTriangleMesh &mesh) {

	// remove zero area faces
	set<vec3i> newIndicesSet;
	vector<bool> verticesFlag(mesh.amount, false);
	for (int faceID = 0; faceID < (int)mesh.indices.size(); faceID++) {
		vec3i idx = mesh.indices[faceID];
		if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2]) continue;
		vec3 pos[3];
		for (int k = 0; k < 3; k++) pos[k] = mesh.positions[idx[k]];
		vec3d v1 = pos[1] - pos[0];
		vec3d v2 = pos[2] - pos[0];
		if (cml::dot(v1, v2) < 0) v2 = -v2;
		if (cml::unsigned_angle(v1, v2) < 1e-5) continue;
		newIndicesSet.insert(idx);
		for(int k=0; k<3; k++) verticesFlag[idx[k]] = true;
	}
	vector<vec3i> newIndices(newIndicesSet.begin(), newIndicesSet.end());

	// remove unreferenced vertices
	vector<vec3> newPositions;
	vector<vec3> newNormals;
	vector<int> verticesMap(mesh.amount);
	for (int vertID = 0; vertID < mesh.amount; vertID++) {
		if (verticesFlag[vertID]) {
			verticesMap[vertID] = (int)newPositions.size();
			newPositions.push_back(mesh.positions[vertID]);
			newNormals.push_back(mesh.normals[vertID]);
		}
	}

	// update face indices
	for (vec3i &idx : newIndices) {
		for (int k = 0; k < 3; k++) idx[k] = verticesMap[idx[k]];
	}


	mesh.positions.swap(newPositions);
	mesh.normals.swap(newNormals);
	mesh.indices.swap(newIndices);
	mesh.amount = (int)mesh.positions.size();

	return true;
}

bool MeshUtil::cleanUp(TTetrahedralMesh &mesh) {

	// compute average edge length
	double totalLength = 0;
	for (vec4i idx : mesh.indices) {
		for (int i = 0; i < 3; i++) {
			vec3d p = mesh.positions[idx[i]];
			for (int j = i + 1; j < 4; j++) {
				vec3d q = mesh.positions[idx[j]];
				totalLength += (p - q).length();
			}
		}
	}
	double averageLength = totalLength / ((int)mesh.indices.size() * 6);
	double volumeEps = pow(averageLength*1e-6, 3.0);

	// remove zero volume tetrahedra
	set<vec4i> newIndicesSet;
	vector<bool> verticesFlag(mesh.amount, false);
	for (vec4i idx : mesh.indices) {
		
		vec3d side[3];
		for (int k = 0; k < 3; k++) side[k] = mesh.positions[idx[k + 1]] - mesh.positions[idx[0]];
		double volume = fabs(cml::dot(side[0], cml::cross(side[1], side[2]))) / 6;
		if (volume < volumeEps) continue;
		newIndicesSet.insert(idx);
		for (int k = 0; k<4; k++) verticesFlag[idx[k]] = true;
	}
	vector<vec4i> newIndices(newIndicesSet.begin(), newIndicesSet.end());

	// remove unreferenced vertices
	vector<vec3> newPositions;
	vector<int> verticesMap(mesh.amount);
	int firstOne = -1;
	for (int vertID = 0; vertID < mesh.amount; vertID++) {
		if (verticesFlag[vertID]) {
			verticesMap[vertID] = (int)newPositions.size();
			newPositions.push_back(mesh.positions[vertID]);
		}
		else {
			if (firstOne < 0) firstOne = vertID;
		}
	}

	// update face indices
	for (vec4i &idx : newIndices) {
		for (int k = 0; k < 4; k++) idx[k] = verticesMap[idx[k]];
	}

	// make consistant index ordering
	for (vec4i &idx : newIndices) {
		vec3d v[3];
		for (int k = 0; k < 3; k++) v[k] = newPositions[idx[k + 1]] - newPositions[idx[0]];
		if (cml::dot(cml::cross(v[0], v[1]), v[2]) < 0) swap(idx[0], idx[1]);
	}

	cout << "Tetrahedral mesh clean up:" << endl;
	cout << "\tremoved " << (mesh.positions.size() - newPositions.size()) << " vertices ";
	cout << " and " << (mesh.indices.size() - newIndices.size()) << " tetrahedra" << endl;

	mesh.normals.clear();
	mesh.positions.swap(newPositions);
	mesh.indices.swap(newIndices);
	mesh.amount = (int)mesh.positions.size();

	return true;
}

bool MeshUtil::buildVertexGraph(TTriangleMesh &mesh, vector<vector<int>> &graph) {

	int numVertices = mesh.amount;

	vector<set<int>> vertNeighbors(numVertices);

	for (vec3i &faceIdx : mesh.indices) {
		for (int k = 0; k < 3; k++) {
			vertNeighbors[faceIdx[k]].insert(faceIdx[(k + 1) % 3]);
			vertNeighbors[faceIdx[k]].insert(faceIdx[(k + 2) % 3]);
		}
	}

	graph.resize(numVertices);
	for (int vertID = 0; vertID < numVertices; vertID++) {
		auto &nb = vertNeighbors[vertID];
		graph[vertID].assign(nb.begin(), nb.end());
	}

	return true;
}

bool MeshUtil::buildFaceGraph(TTriangleMesh &mesh, vector<vector<int>> &graph) {

	int numFaces = (int)mesh.indices.size();

	// extract all edges and faces touching them

	map<vec2i, set<int>> edgeMap; // (vertex ID, vertex ID) => face ID set
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i faceIdx = mesh.indices[faceID];
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
			edgeMap[edgeKey].insert(faceID);
		}
	}

	// extract neighbor faces

	vector<set<int>> faceNeighbors(numFaces);
	for (auto &it : edgeMap) {
		vector<int> faceList(it.second.begin(), it.second.end());
		int n = (int)faceList.size();
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				faceNeighbors[faceList[i]].insert(faceList[j]);
				faceNeighbors[faceList[j]].insert(faceList[i]);
			}
		}
	}

	// organize graph

	graph.resize(numFaces);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		graph[faceID].assign(faceNeighbors[faceID].begin(), faceNeighbors[faceID].end());
	}

	return true;
}

bool MeshUtil::subdivideMesh(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double radius)
{

	// clone existing vertices
	outMesh.positions.assign(inMesh.positions.begin(), inMesh.positions.end());
	outMesh.normals.assign(inMesh.normals.begin(), inMesh.normals.end());

	// compute face area
	int numFaces = (int)inMesh.indices.size();

	float maxLength;
	if (radius) {
		maxLength = (float)(radius * 2);
	} else {
		// compute radius from default configuration
		float totalArea = 0;
		for (int faceID = 0; faceID < numFaces; faceID++) {
			vec3i faceIdx = inMesh.indices[faceID];
			vec3 facePos[3];
			for (int j = 0; j < 3; j++) facePos[j] = inMesh.positions[faceIdx[j]];
			float area = cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]).length();
			totalArea += area;
		}
		int numSamples = StyleSynthesisConfig::mSample_WholeMeshSampleNumber;
		maxLength = sqrt(totalArea / (sqrt(3.0f)*numSamples / 2)) * 0.76f * 2;
	}

	// process each triangle
	outMesh.indices.clear();
	outFaceIndices.clear();
	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = inMesh.indices[faceID];
		vec3 facePos[3];
		for (int j = 0; j < 3; j++) facePos[j] = inMesh.positions[faceIdx[j]];
		float length = 0;
		for (int j = 0; j < 3; j++) length = max(length, (facePos[(j + 1) % 3] - facePos[j]).length());

		if(length > maxLength) {
			// subdivide long triangle
			int div = (int)ceil(length / maxLength);

			// add new vertices
			vector<vector<int>> divPointIdx(div+1, vector<int>(div+1));
			for (int j = 0; j <= div; j++) {
				for (int v = 0; v <= j; v++) {
					int u = j - v;
					if (u == 0 && v == 0) divPointIdx[u][v] = faceIdx[0];
					else if (u == div && v == 0) divPointIdx[u][v] = faceIdx[1];
					else if (u == 0 && v == div) divPointIdx[u][v] = faceIdx[2];
					else {
						divPointIdx[u][v] = (int)outMesh.positions.size();						
						outMesh.positions.push_back(facePos[0] + ((facePos[1] - facePos[0]) * u + (facePos[2] - facePos[0]) * v) / div);
						outMesh.normals.push_back(vec3()); // will calculate later
					}

				}
			}

			// add new faces
			for (int j = 0; j <= div; j++) {
				for (int v = 0; v <= j; v++) {
					int u = j - v;
					if (j < div) {
						outMesh.indices.push_back(vec3i(
							divPointIdx[u][v],
							divPointIdx[u+1][v],
							divPointIdx[u][v+1]));
						outFaceIndices.push_back(faceID);
						if (u == 0) continue;
						outMesh.indices.push_back(vec3i(
							divPointIdx[u][v],
							divPointIdx[u][v+1],
							divPointIdx[u-1][v+1]));
						outFaceIndices.push_back(faceID);
					}
				}
			}
		} else {
			if (length == 0) continue; // skip zero area faces
			// retain original triangle
			outMesh.indices.push_back(faceIdx);
			outFaceIndices.push_back(faceID);
		}
	}
	outMesh.amount = (int)outMesh.positions.size();

	// compute normal
	if (!recomputeNormals(outMesh)) return false;

	return true;
}

bool MeshUtil::subdivideMeshKeepTopology(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double radius)
{

	// clone existing vertices
	outMesh.positions.assign(inMesh.positions.begin(), inMesh.positions.end());
	outMesh.normals.assign(inMesh.normals.begin(), inMesh.normals.end());

	int numFaces = (int)inMesh.indices.size();

	// compute max length of edge

	float maxLength;
	if (radius) {
		maxLength = (float)(radius * 2);
	} else {
		// compute radius from default configuration
		float totalArea = 0;
		for (int faceID = 0; faceID < numFaces; faceID++) {
			vec3i faceIdx = inMesh.indices[faceID];
			vec3 facePos[3];
			for (int j = 0; j < 3; j++) facePos[j] = inMesh.positions[faceIdx[j]];
			float area = cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]).length();
			totalArea += area;
		}
		int numSamples = StyleSynthesisConfig::mSample_WholeMeshSampleNumber;
		maxLength = sqrt(totalArea / (sqrt(3.0f)*numSamples / 2)) * 0.76f;
	}

	// divide each edge (add new vertices)

	vector<vector<int>> edgePointList; // point ID : # points on edge : # of edges
	map<vec2i, int> edgeMap; // edge ID : point ID pair (small first)

	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = inMesh.indices[faceID];
		for (int j = 0; j < 3; j++) {
			int pID1 = faceIdx[j];
			int pID2 = faceIdx[(j + 1) % 3];
			if (pID1 > pID2) swap(pID1, pID2);
			vec2i edgeKey(pID1, pID2);
			auto it = edgeMap.find(edgeKey);
			if (it == edgeMap.end()) {
				vec3 p1 = inMesh.positions[pID1];
				vec3 n1 = inMesh.normals[pID1];
				vec3 p2 = inMesh.positions[pID2];				
				vec3 n2 = inMesh.normals[pID2];
				float length = (p2 - p1).length();
				vector<int> pointList;
				pointList.push_back(pID1);
				if (length > maxLength) {
					// divide edge (add divide points)
					int div = (int)(ceil(length / maxLength));					
					for (int k = 1; k < div; k++) {
						float r = k / (float)div;
						vec3 pNew = p2*r + p1*(1 - r);
						vec3 nNew = cml::normalize(n2*r + n1*(1 - r));
						int iNew = (int)outMesh.positions.size();
						outMesh.positions.push_back(pNew);
						outMesh.normals.push_back(nNew);
						pointList.push_back(iNew);
					}					
				}
				pointList.push_back(pID2);
				edgeMap[edgeKey] = (int)edgePointList.size();
				edgePointList.push_back(pointList);
			}
		}
	}

	// first pass: divide triangle along one edge (the one to the vertex with smallest ID)

	vector<vec3i> tmpFaces(0);
	vector<int> tmpFaceIndices(0);

	for (int faceID = 0; faceID < numFaces; faceID++) {

		vec3i faceIdx = inMesh.indices[faceID];

		// rotate indices to make first index smallest
		if (faceIdx[1] <= faceIdx[0] && faceIdx[1] <= faceIdx[2]) {
			faceIdx = vec3i(faceIdx[1], faceIdx[2], faceIdx[0]);
		} else if (faceIdx[2] <= faceIdx[0] && faceIdx[2] <= faceIdx[1]) {
			faceIdx = vec3i(faceIdx[2], faceIdx[0], faceIdx[1]);
		}

		int topPointID = faceIdx[0];
		vec2i bottomEdgeKey = vec2i(faceIdx[1], faceIdx[2]);
		bool isSwapped = false;
		if (bottomEdgeKey[0] > bottomEdgeKey[1]) {
			swap(bottomEdgeKey[0], bottomEdgeKey[1]);
			isSwapped = true;
		}
		int bottomEdgeID = edgeMap[bottomEdgeKey];
		auto bottomEdgePointList = edgePointList[bottomEdgeID];
		for (int k = 0; k < (int)bottomEdgePointList.size() - 1; k++) {
			if (k > 0) {
				// add internal edge
				int pID1 = topPointID;
				int pID2 = bottomEdgePointList[k];
				vec2i edgeKey(pID1, pID2);
				vec3 p1 = outMesh.positions[pID1];
				vec3 n1 = outMesh.normals[pID1];
				vec3 p2 = outMesh.positions[pID2];
				vec3 n2 = outMesh.normals[pID2];
				float length = (p2 - p1).length();
				vector<int> pointList;
				pointList.push_back(pID1);
				if (length > maxLength) {
					// divide edge (add divide points)
					int div = (int)(ceil(length / maxLength));
					for (int k = 1; k < div; k++) {
						float r = k / (float)div;
						vec3 pNew = p2*r + p1*(1 - r);
						vec3 nNew = cml::normalize(n2*r + n1*(1 - r));
						int iNew = (int)outMesh.positions.size();
						outMesh.positions.push_back(pNew);
						outMesh.normals.push_back(nNew);
						pointList.push_back(iNew);
					}
				}
				pointList.push_back(pID2);
				edgeMap[edgeKey] = (int)edgePointList.size();
				edgePointList.push_back(pointList);
			}
			vec3i face = vec3i(topPointID, bottomEdgePointList[k], bottomEdgePointList[k + 1]);
			if (isSwapped) swap(face[1], face[2]);
			tmpFaces.push_back(face);
			tmpFaceIndices.push_back(faceID);
		}
	}

	// second pass: divide triangle along two divided edge

	outMesh.indices.clear();
	outFaceIndices.clear();

	for (int faceID = 0; faceID < (int)tmpFaces.size(); faceID++) {

		vec3i faceIdx = tmpFaces[faceID];
		int edgeID1 = edgeMap[vec2i(faceIdx[0], faceIdx[1])];
		int edgeID2 = edgeMap[vec2i(faceIdx[0], faceIdx[2])];
		auto &pointList1 = edgePointList[edgeID1];
		auto &pointList2 = edgePointList[edgeID2];
		if (pointList1[0] != faceIdx[0] || pointList2[0] != faceIdx[0]) {
			cout << "Error: incorrect point list" << endl;
			return false;
		}
		int n1 = (int)pointList1.size();
		int n2 = (int)pointList2.size();
		int p1 = 1;
		int p2 = 1;
		float l1 = (outMesh.positions[pointList1[p1]] - outMesh.positions[pointList1[0]]).length_squared();
		float l2 = (outMesh.positions[pointList2[p2]] - outMesh.positions[pointList2[0]]).length_squared();
		outMesh.indices.push_back(vec3i(faceIdx[0], pointList1[p1], pointList2[p2]));
		outFaceIndices.push_back(tmpFaceIndices[faceID]);
		while (true) {
			int oldP1 = p1;
			int oldP2 = p2;
			if (p1 == n1 - 1) p2++;
			else if (p2 == n2 - 1) p1++;
			else if (l1 < l2) p1++;
			else p2++;
			if (p1 >= n1 || p2 >= n2) {
				break;
			}
			if (oldP1 != p1) {
				l1 = (outMesh.positions[pointList1[p1]] - outMesh.positions[pointList1[0]]).length_squared();
				outMesh.indices.push_back(vec3i(pointList1[oldP1], pointList1[p1], pointList2[oldP2]));
				outFaceIndices.push_back(tmpFaceIndices[faceID]);
			} else {
				l2 = (outMesh.positions[pointList2[p2]] - outMesh.positions[pointList2[0]]).length_squared();
				outMesh.indices.push_back(vec3i(pointList1[oldP1], pointList2[p2], pointList2[oldP2]));
				outFaceIndices.push_back(tmpFaceIndices[faceID]);
			}
		}
	}

	outMesh.amount = (int)outMesh.positions.size();

	// compute normal
	if (!recomputeNormals(outMesh)) return false;

	return true;
}

bool MeshUtil::subdivideMeshMidPoint(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double radius)
{
	TTriangleMesh tmpMesh;
	tmpMesh.positions = inMesh.positions;
	tmpMesh.normals = inMesh.normals;

	outFaceIndices.clear();
	map<vec2i, int> edgeMap;
	vec4i caseCount(0, 0, 0, 0);
	for (int faceID = 0; faceID < (int)inMesh.indices.size(); faceID++) {
		vec3i facePoints = inMesh.indices[faceID];

		/*
		// skip dividing faces with short edges
		bool canSubdivide = true;
		for (int k = 0; k < 3; k++) {
			float edgeLen = (inMesh.positions[facePoints[k]] - inMesh.positions[facePoints[(k + 1) % 3]]).length();
			if (edgeLen < threshold) {
				canSubdivide = false;
				break;
			}
		}
		if (!canSubdivide) {
			tmpMesh.indices.push_back(facePoints);
			outFaceIndices.push_back(faceID);
			continue;
		}
		*/

		// compute mid point on edges
		vec3i midPoints(-1,-1,-1);
		for (int k = 0; k < 3; k++) {
			vec2i edgeID(facePoints[k], facePoints[(k + 1) % 3]);
			if (edgeID[0] > edgeID[1]) swap(edgeID[0], edgeID[1]);
			auto it = edgeMap.find(edgeID);
			if (it == edgeMap.end()) {
				vec3 p0 = inMesh.positions[edgeID[0]];
				vec3 p1 = inMesh.positions[edgeID[1]];
				float edgeLen = (p0 - p1).length();
				if (edgeLen < radius) {
					edgeMap[edgeID] = -1;
					midPoints[k] = -1;
				} else {
					vec3 midPosition = (p0 + p1) / 2;
					vec3 midNormal = inMesh.normals[edgeID[0]] + inMesh.normals[edgeID[1]];
					if (midNormal.length_squared() > 0) midNormal.normalize();
					int midIndex = (int)tmpMesh.positions.size();
					tmpMesh.positions.push_back(midPosition);
					tmpMesh.normals.push_back(midNormal);
					edgeMap[edgeID] = midIndex;
					midPoints[k] = midIndex;
				}
			} else {
				midPoints[k] = it->second;
			}
		}
		// handle different cases
		int collapseCount = 0;
		for (int k = 0; k < 3; k++) if (midPoints[k] < 0) collapseCount++;
		if (collapseCount == 0) {
			tmpMesh.indices.push_back(vec3i(facePoints[0], midPoints[0], midPoints[2]));
			tmpMesh.indices.push_back(vec3i(midPoints[0], facePoints[1], midPoints[1]));
			tmpMesh.indices.push_back(vec3i(midPoints[2], midPoints[1], facePoints[2]));
			tmpMesh.indices.push_back(vec3i(midPoints[0], midPoints[1], midPoints[2]));
			for (int k = 0; k < 4; k++) outFaceIndices.push_back(faceID);
			caseCount[0]++;
		} else if (collapseCount == 1) {
			int collapseID = -1;
			for (int k = 0; k < 3; k++) if (midPoints[k] < 0) collapseID = k;
			vec3i i(collapseID, (collapseID + 1) % 3, (collapseID + 2) % 3);
			tmpMesh.indices.push_back(vec3i(facePoints[i[2]], midPoints[i[2]], midPoints[i[1]]));
			float l1 = (tmpMesh.positions[facePoints[i[0]]] - tmpMesh.positions[midPoints[i[1]]]).length();
			float l2 = (tmpMesh.positions[facePoints[i[1]]] - tmpMesh.positions[midPoints[i[2]]]).length();
			if (l1 < l2) {
				tmpMesh.indices.push_back(vec3i(facePoints[i[0]], midPoints[i[1]], midPoints[i[2]]));
				tmpMesh.indices.push_back(vec3i(facePoints[i[0]], facePoints[i[1]], midPoints[i[1]]));
			} else {
				tmpMesh.indices.push_back(vec3i(facePoints[i[0]], facePoints[i[1]], midPoints[i[2]]));
				tmpMesh.indices.push_back(vec3i(facePoints[i[1]], midPoints[i[1]], midPoints[i[2]]));
			}
			for (int k = 0; k < 3; k++) outFaceIndices.push_back(faceID);
			caseCount[1]++;
		} else if (collapseCount == 2) {
			int subID = -1;
			for (int k = 0; k < 3; k++) if (midPoints[k] >= 0) subID = k;
			vec3i i(subID, (subID + 1) % 3, (subID + 2) % 3);
			tmpMesh.indices.push_back(vec3i(facePoints[i[0]], midPoints[i[0]], facePoints[i[2]]));
			tmpMesh.indices.push_back(vec3i(midPoints[i[0]], facePoints[i[1]], facePoints[i[2]]));
			for (int k = 0; k < 2; k++) outFaceIndices.push_back(faceID);
			caseCount[2]++;
		} else if (collapseCount == 3) {
			tmpMesh.indices.push_back(facePoints);
			outFaceIndices.push_back(faceID);
			caseCount[3]++;
		}
	}

	outMesh.positions.swap(tmpMesh.positions);
	outMesh.normals.swap(tmpMesh.normals);
	outMesh.indices.swap(tmpMesh.indices);
	outMesh.amount = (int)outMesh.positions.size();

	return true;
}

bool MeshUtil::splitMeshIntoConnectedComponents(
	TTriangleMesh &inMesh,
	vector<TTriangleMesh> &outMeshes,
	vector<vector<int>> &outFaceIndices)
{
	int numVertices = inMesh.amount;
	int numFaces = (int)inMesh.indices.size();

	// build graph

	vector<unordered_set<int>> vertexGraph(numVertices);
	for (vec3i idx : inMesh.indices) {
		for (int k = 0; k < 3; k++) {
			vertexGraph[idx[k]].insert(idx[(k + 1) % 3]);
			vertexGraph[idx[k]].insert(idx[(k + 2) % 3]);
		}
	}

	// label connected vertices

	int numComponents = 0;
	vector<int> vertexLabels(numVertices, -1);
	vector<int> vertexMap(numVertices, -1);
	outMeshes.clear();

	for (int vertID = 0; vertID < numVertices; vertID++) {

		if (vertexLabels[vertID] >= 0) continue;
		vertexLabels[vertID] = numComponents;
		vertexMap[vertID] = 0;

		// BFS

		vector<int> queue(1, vertID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			for (int neighborID : vertexGraph[currentID]) {
				if (vertexLabels[neighborID] >= 0) continue;
				vertexLabels[neighborID] = numComponents;
				vertexMap[neighborID] = (int)queue.size();
				queue.push_back(neighborID);
			}
			head++;
		}

		// add component

		outMeshes.push_back(TTriangleMesh());
		TTriangleMesh &mesh = outMeshes.back();
		mesh.positions.clear();
		mesh.normals.clear();
		mesh.indices.clear();
		for (int currentID : queue) {
			mesh.positions.push_back(inMesh.positions[currentID]);
			mesh.normals.push_back(inMesh.normals[currentID]);
		}
		mesh.amount = (int)mesh.positions.size();

		numComponents++;
	}

	// label faces

	outFaceIndices.assign(numComponents, vector<int>(0));
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i oldFaceIdx = inMesh.indices[faceID];
		int compID = vertexLabels[oldFaceIdx[0]];
		vec3i newFaceIdx = vec3i(vertexMap[oldFaceIdx[0]], vertexMap[oldFaceIdx[1]], vertexMap[oldFaceIdx[2]]);
		outMeshes[compID].indices.push_back(newFaceIdx);
		outFaceIndices[compID].push_back(faceID);
	}

	return true;
}

bool MeshUtil::splitMeshIntoCleanedConnectedComponents(
	TTriangleMesh &inMesh,
	vector<TTriangleMesh> &outMeshes)
{

	// get raw connected components

	vector<TTriangleMesh> rawComponents;
	vector<vector<int>> rawComponentFaceIndices;
	if (!splitMeshIntoConnectedComponents(inMesh, rawComponents, rawComponentFaceIndices)) return false;
	int numRawComponents = (int)rawComponents.size();

	// sort connected components by face number

	vector<int> sortedComponentOrder(numRawComponents);
	for (int k = 0; k < numRawComponents; k++) sortedComponentOrder[k] = k;
	sort(sortedComponentOrder.begin(), sortedComponentOrder.end(),
		[&rawComponents](int lhs, int rhs) {
			return rawComponents[lhs].indices.size() > rawComponents[rhs].indices.size();
		});

	// weld mesh

	TTriangleMesh weldedMesh;
	vector<int> weldedFaceIndices;
	if (true) {
		// NOTE: don't really weld it -- hard to handle degenerated faces
		//if (!weldMeshFaces(inMesh, weldedMesh, weldedFaceIndices)) return false;
		TTriangleMesh tmpMesh;
		vector<int> vertexIndices;
		if (!removeDuplicateVertices(inMesh, tmpMesh, vertexIndices)) return false;
		if (!removeDuplicateFaces(tmpMesh, weldedMesh, weldedFaceIndices)) return false;
	}
	int numWeldedFaces = (int)weldedMesh.indices.size();

	// label welded mesh faces by raw connected components
	
	vector<int> weldedFaceLabel(numWeldedFaces, -1);

	int numWeldedComponents = 0;
	for (int orderID = 0; orderID < numRawComponents; orderID++) { // check large components first
		int rawCompID = sortedComponentOrder[orderID];
		bool hasNewComponent = false;
		for (int inFaceID : rawComponentFaceIndices[rawCompID]) {
			int weldedFaceID = weldedFaceIndices[inFaceID];
			if (weldedFaceID < 0) {
				cout << "Error: incorrect welded face index" << endl;
				return false;
			}
			int &label = weldedFaceLabel[weldedFaceID];
			if (label < 0) {
				label = numWeldedComponents;
				hasNewComponent = true;
			}
		}
		if (hasNewComponent) numWeldedComponents++;
	}

	// output components

	vector<vector<int>> weldedComponentFaces(numWeldedComponents, vector<int>(0));
	for (int faceID = 0; faceID < numWeldedFaces; faceID++) {
		int label = weldedFaceLabel[faceID];
		if (label < 0) continue;
		weldedComponentFaces[label].push_back(faceID);
	}

	outMeshes.resize(numWeldedComponents);
	int numOutMeshes = 0;
	for (int compID = 0; compID < numWeldedComponents; compID++) {
		if (!extractSubMesh(weldedMesh, weldedComponentFaces[compID], outMeshes[numOutMeshes])) return false;
		double meshArea;
		if (!computeFaceArea(outMeshes[numOutMeshes], meshArea)) return false;
		if (meshArea) numOutMeshes++; // skip empty component
	}
	outMeshes.resize(numOutMeshes);
	
	return true;
}

bool MeshUtil::splitMeshAlongCrease(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> *outLabels)
{
	// build edge set

	map<vec2i, vec2i> edgeMap; // edge vertex ID pair => face ID pair
	for (int faceID = 0; faceID < (int)inMesh.indices.size(); faceID++) {
		vec3i faceIdx = inMesh.indices[faceID];
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
			auto it = edgeMap.find(edgeKey);
			if (it == edgeMap.end()) {
				edgeMap[edgeKey] = vec2i(faceID, -1);
			} else {
				if (it->second[0] < 0) {
					// already marked as invalid - skip
				} else if (it->second[1] < 0) {
					it->second[1] = faceID;
				} else {
					// non-manifold edge - mark it as invalid
					it->second = vec2i(-1, -1);
				}
			}
		}
	}

	// check crease edge

	for (auto &it : edgeMap) {
		if (it.second[0] < 0 || it.second[1] < 0) continue;
		vec3 faceCenters[2];
		vec3 faceNormals[2];
		for (int alterID = 0; alterID < 2; alterID++) {
			int faceID = it.second[alterID];
			vec3i faceIdx = inMesh.indices[faceID];
			vec3 facePos[3];
			for (int k = 0; k < 3; k++) facePos[k] = inMesh.positions[faceIdx[k]];
			faceNormals[alterID] = cml::normalize(cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]));
			faceCenters[alterID] = (facePos[0] + facePos[1] + facePos[2]) / 3;
		}
		double dihAngle;
		if (!computeDihedralAngle(faceCenters[0], faceNormals[0], faceCenters[1], faceNormals[1], dihAngle)) return false;

		if (dihAngle >= cml::rad(269.0)) { // UNDONE: param crease angle threshold
			it.second = vec2i(-1, -1);
		}
	}

	// build graph

	vector<set<int>> faceNeighbors(inMesh.indices.size());
	for (auto &it : edgeMap) {
		if (it.second[0] >= 0 && it.second[1] >= 0) {
			faceNeighbors[it.second[0]].insert(it.second[1]);
			faceNeighbors[it.second[1]].insert(it.second[0]);
		}
	}

	// extract connected components

	vector<int> newVertexIndices(0);
	outMesh.indices.resize(inMesh.indices.size());

	int numComponents = 0;
	vector<int> visitedFlag(inMesh.indices.size(), false);
	if (outLabels) outLabels->resize(inMesh.indices.size());
	for (int faceID = 0; faceID < (int)inMesh.indices.size(); faceID++) {
		if (visitedFlag[faceID]) continue;
		visitedFlag[faceID] = true;

		// BFS
		vector<int> queue(1, faceID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentFace = queue[head];
			for (int neighborFace : faceNeighbors[currentFace]) {
				if (!visitedFlag[neighborFace]) {
					visitedFlag[neighborFace] = true;
					queue.push_back(neighborFace);
				}
			}
			head++;
		}

		// organize face data
		map<int, int> vertexMap;
		for (int faceID : queue) {
			vec3i oldIdx = inMesh.indices[faceID];
			vec3i newIdx;
			for (int k = 0; k < 3; k++) {
				auto it = vertexMap.find(oldIdx[k]);
				if (it == vertexMap.end()) {
					newIdx[k] = (int)newVertexIndices.size();
					vertexMap[oldIdx[k]] = newIdx[k];
					newVertexIndices.push_back(oldIdx[k]);
				} else {
					newIdx[k] = it->second;
				}
			}
			outMesh.indices[faceID] = newIdx;
			if (outLabels) (*outLabels)[faceID] = numComponents;
		}

		numComponents++;
	}

	// organize vertex data
	outMesh.positions.resize(newVertexIndices.size());
	outMesh.normals.resize(newVertexIndices.size());
	for (int vertID = 0; vertID < (int)newVertexIndices.size(); vertID++) {
		int newID = newVertexIndices[vertID];
		outMesh.positions[vertID] = inMesh.positions[newID];
		outMesh.normals[vertID] = inMesh.normals[newID];
	}
	outMesh.amount = (int)outMesh.positions.size();

	return true;
}

bool MeshUtil::extractSubMesh(
	TTriangleMesh &inMesh,
	vector<int> &subMeshFaces,
	TTriangleMesh &outMesh)
{

	// get all vertices

	set<int> newVertexSet;
	for (int faceID : subMeshFaces) {
		vec3i faceIdx = inMesh.indices[faceID];
		for (int k = 0; k < 3; k++) newVertexSet.insert(faceIdx[k]);
	}

	vector<int> newVertexList(newVertexSet.begin(), newVertexSet.end());
	int numNewVertices = (int)newVertexList.size();

	// output vertices

	vector<int> vertexMap(inMesh.amount, -1);
	outMesh.positions.resize(numNewVertices);
	outMesh.normals.resize(numNewVertices);
	for (int newID = 0; newID < numNewVertices; newID++) {
		int oldID = newVertexList[newID];
		vertexMap[oldID] = newID;
		outMesh.positions[newID] = inMesh.positions[oldID];
		outMesh.normals[newID] = inMesh.normals[oldID];
	}

	// output faces

	outMesh.indices.resize(subMeshFaces.size());
	for (int faceID = 0; faceID < (int)subMeshFaces.size(); faceID++) {
		vec3i &oldIdx = inMesh.indices[subMeshFaces[faceID]];
		vec3i &newIdx = outMesh.indices[faceID];
		for (int k = 0; k < 3; k++) newIdx[k] = vertexMap[oldIdx[k]];
	}

	outMesh.amount = numNewVertices;

	return true;
}

bool MeshUtil::extractMeshFromSamplePoints(
	TTriangleMesh &inMesh, // subdivided mesh
	TSampleSet &inSamples, // all sample points
	vector<bool> &inFlags, // valid sample flag : # sample points
	TTriangleMesh &outMesh) // subset of inMesh
{
	int numFaces = (int)inMesh.indices.size();

	// build tree
	SKDTree tree;
	SKDTreeData treeData;
	if (!SampleUtil::buildKdTree(inSamples.positions, tree, treeData)) return false;
	
	// build query point sets (center of faces)
	Eigen::Matrix3Xd faceCenters(3, numFaces);
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i idx = inMesh.indices[faceID];
		vec3 center = (inMesh.positions[idx[0]] + inMesh.positions[idx[1]] + inMesh.positions[idx[2]]) / 3;
		faceCenters.col(faceID) = Eigen::Vector3d(center[0], center[1], center[2]);
	}

	// create sub-mesh
	Eigen::VectorXi faceSamples;
	if (!SampleUtil::findNearestNeighbors(tree, faceCenters, faceSamples)) return false;
	outMesh = inMesh;
	outMesh.indices.clear();
	for (int faceID = 0; faceID < numFaces; faceID++) {
		if (inFlags[faceSamples[faceID]]) outMesh.indices.push_back(inMesh.indices[faceID]);
	}

	return true;
}

bool MeshUtil::weldMeshFaces(
	TTriangleMesh &inMesh,
	TTriangleMesh &outMesh,
	vector<int> &outFaceIndices,
	double eps)
{
	// compute eps

	if (eps < 0) {
		vec3 bbMin, bbMax;
		if (!computeAABB(inMesh, bbMin, bbMax)) return false;
		float bbLen = (bbMax - bbMin).length();
		eps = bbLen * 1e-5; // UNDONE: param eps
	}

	// remove duplicate vertices/faces

	TTriangleMesh tmpMesh;
	vector<int> vertexIndices;
	vector<int> dupliFaceIndices;
	vector<int> degenFaceIndices;
	if (!removeDuplicateVertices(inMesh, tmpMesh, vertexIndices, eps)) return false;
	if (!removeDuplicateFaces(tmpMesh, outMesh, dupliFaceIndices, eps)) return false;
	if (!removeDegeneratedFaces(outMesh, outMesh, degenFaceIndices)) return false;

	outFaceIndices.resize(dupliFaceIndices.size());
	for (int k = 0; k < (int)outFaceIndices.size(); k++) {
		outFaceIndices[k] = degenFaceIndices[dupliFaceIndices[k]];
	}

	return true;
}

bool MeshUtil::reorientComponentFaces(TTriangleMesh &inoutMesh) {

	// assumes welded and connected mesh component

	vec3 bbMin, bbMax;
	if (!MeshUtil::computeAABB(inoutMesh, bbMin, bbMax)) return false;
	float eps = (bbMax - bbMin).length() * 1e-4f;

	int numVertices = inoutMesh.amount;
	int numFaces = (int)inoutMesh.indices.size();

	TKDTree tree;
	TKDTreeData treeData;
	if (!MeshUtil::buildKdTree(inoutMesh, tree, treeData)) return false;

	// build face neighboring graph

	map<vec2i, vector<int>> edgeFaces; // edge key => face ID list
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i &faceIdx = inoutMesh.indices[faceID];
		for (int k = 0; k < 3; k++) {
			vec2i edgeKey(faceIdx[k], faceIdx[(k + 1) % 3]);
			if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
			auto it = edgeFaces.find(edgeKey);
			if (it == edgeFaces.end()) {
				edgeFaces[edgeKey].assign(1, faceID);
			} else {
				it->second.push_back(faceID);
			}
		}
	}

	// compute face visibility

	int numRays = 16;

	vector<float> faceVisibility(numFaces);
	vector<vec3> faceNormal(numFaces);

#pragma omp parallel for
	for (int faceID = 0; faceID < numFaces; faceID++) {
		vec3i &faceIdx = inoutMesh.indices[faceID];
		vec3 faceP[4];
		for (int k = 0; k < 3; k++) {
			faceP[k] = inoutMesh.positions[faceIdx[k]];
		}
		faceP[3] = (faceP[0] + faceP[1] + faceP[2]) / 3;
		for (int k = 0; k < 3; k++) faceP[k] = (faceP[k] + faceP[3]) / 2;
		vec3 faceN = cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]);
		if (faceN.length_squared()) faceN.normalize();
		faceNormal[faceID] = faceN; // direction does not matter

		float maxVisibility = 0;
		bool flipFace = false;
		for (int k = 0; k < 4; k++) {
			int visibleCountP = 0;
			int visibleCountN = 0;
			for (int rayID = 0; rayID < numRays; rayID++) {
				double r1 = cml::random_unit();
				double r2 = cml::random_unit();
				vec3 rayDirP = vec3d(
					2.0 * cos(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
					2.0 * sin(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
					1.0 - 2.0*r2); // random direction on unit sphere
				if (cml::dot(rayDirP, faceN) < 0) rayDirP = -rayDirP;
				vec3 rayDirN = -rayDirP;

				vec3 rayOriginP = faceP[k] + rayDirP * eps; // add eps for offset
				vec3 rayOriginN = faceP[k] + rayDirN * eps;
				Thea::Ray3 rayP(G3D::Vector3(rayOriginP.data()), G3D::Vector3(rayDirP.data()));
				Thea::Ray3 rayN(G3D::Vector3(rayOriginN.data()), G3D::Vector3(rayDirN.data()));
				if (tree.rayIntersectionTime(rayP) < -0.5) visibleCountP++;
				if (tree.rayIntersectionTime(rayN) < -0.5) visibleCountN++;
			}
			float visibility = max(visibleCountP, visibleCountN) / (float)numRays;
			if (visibility > maxVisibility) {
				maxVisibility = visibility;
				flipFace = visibleCountN > visibleCountP;
			}
		}

		faceVisibility[faceID] = maxVisibility;
		if (flipFace) swap(faceIdx[1], faceIdx[2]);
	}

	// iteratively flip faces based on neighboring faces

	bool finish = false;
	int iteration = 0;
	while (!finish) {
		finish = true;
		vector<vec3i> newIndices = inoutMesh.indices;
		vector<bool> flipFlags(numFaces, false);
		int flipCount = 0;
#pragma omp parallel for shared(flipCount, flipFlags)
		for (int faceID = 0; faceID < numFaces; faceID++) {
			vec3i faceIdx = inoutMesh.indices[faceID];
			float flipVote = 0;
			float unflipVote = 0;
			bool canFlip = true;
			for (int k = 0; k < 3; k++) {
				vec2i edge(faceIdx[k], faceIdx[(k + 1) % 3]);
				vec2i edgeKey = edge;
				if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
				auto &nbList = edgeFaces[edgeKey];
				for (int nbID : nbList) {
					if (nbID == faceID) continue;
					if (flipFlags[nbID]) canFlip = false;
					float normalCosine = fabs(cml::dot(faceNormal[faceID], faceNormal[nbID]));
					vec3i &nbIdx = inoutMesh.indices[nbID];
					int ptr = 0;
					for (ptr = 0; ptr < 3; ptr++) if (nbIdx[ptr] == edge[0]) break;
					if (nbIdx[(ptr + 1) % 3] == edge[1]) {
						flipVote += (faceVisibility[nbID] + 0.1f);// *normalCosine;
					}
					else {
						unflipVote += (faceVisibility[nbID] + 0.1f);// *normalCosine;
					}
				}
			}
			if (canFlip && flipVote > unflipVote) {
				swap(newIndices[faceID][1], newIndices[faceID][2]);
				flipFlags[faceID] = true;
				finish = false;
#pragma omp atomic
				flipCount++;
			}
		}

		if (!finish) inoutMesh.indices.swap(newIndices);

		//cout << "Iteration " << iteration << " : flipped " << flipCount << " faces" << endl;
		iteration++;
		if (iteration > 100) break; // UNDONE: param maximum iteration

		//if (!MeshUtil::saveMesh(outMeshName, mesh)) return false;
		//system("pause");
	}

	return true;
	
	// sort face visibility

	vector<int> faceOrder(numFaces);
	for (int k = 0; k < numFaces; k++) faceOrder[k] = k;
	sort(faceOrder.begin(), faceOrder.end(),
		[&faceVisibility](int id1, int id2){return faceVisibility[id1] > faceVisibility[id2]; });

	// determine face orientations by chaining

	vector<bool> faceProcessedFlag(numFaces, false);
	for (int orderID = 0; orderID < numFaces; orderID++) {
		int faceID = faceOrder[orderID];
		if (faceProcessedFlag[faceID]) continue;
		faceProcessedFlag[faceID] = true;

		// process all adjacent faces (BFS)

		vector<int> queue(1, faceID);
		int head = 0;
		while (head < (int)queue.size()) {
			int currentID = queue[head];
			vec3i faceIdx = inoutMesh.indices[currentID];
			for (int k = 0; k < 3; k++) {
				vec2i edge(faceIdx[k], faceIdx[(k + 1) % 3]);
				vec2i edgeKey = edge;
				if (edgeKey[0] > edgeKey[1]) swap(edgeKey[0], edgeKey[1]);
				auto &nbList = edgeFaces[edgeKey];
				for (int nbID : nbList) {
					if (faceProcessedFlag[nbID]) continue;
					faceProcessedFlag[nbID] = true;
					queue.push_back(nbID);
					vec3i &nbIdx = inoutMesh.indices[nbID];
					int ptr = 0;
					for (ptr = 0; ptr < 3; ptr++) if (nbIdx[ptr] == edge[0]) break;
					if (nbIdx[(ptr + 1) % 3] == edge[1]) {
						swap(nbIdx[1], nbIdx[2]); // flip face
					}
				}
			}
			head++;
		}
	}	
	
	return true;
}

bool MeshUtil::reorientMeshFaces(TTriangleMesh &inoutMesh) {

	// assumes welded mesh

	// compute face normals

	vector<vec3> faceNormals(inoutMesh.indices.size());
	for (int faceID = 0; faceID < (int)inoutMesh.indices.size(); faceID++) {
		vec3i faceIdx = inoutMesh.indices[faceID];
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = inoutMesh.positions[faceIdx[k]];
		faceNormals[faceID] = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
	}

	// split in connected components

	vector<TTriangleMesh> components;
	vector<vector<int>> componentFaceIndices;
	if (!splitMeshIntoConnectedComponents(inoutMesh, components, componentFaceIndices)) return false;

	// re-orient each component

	int numComponent = (int)components.size();
	for (int compID = 0; compID < numComponent; compID++) {
		auto &component = components[compID];
		auto &faceIndices = componentFaceIndices[compID];
		if (!reorientComponentFaces(component)) return false;

		// determine flipping on original face
		for (int compFaceID = 0; compFaceID < (int)faceIndices.size(); compFaceID++) {
			int origFaceID = faceIndices[compFaceID];
			vec3i compFaceIdx = component.indices[compFaceID];
			vec3 compFaceP[3];
			for (int k = 0; k < 3; k++) compFaceP[k] = component.positions[compFaceIdx[k]];
			vec3 compFaceN = cml::normalize(cml::cross(compFaceP[1] - compFaceP[0], compFaceP[2] - compFaceP[0]));
			if (cml::dot(compFaceN, faceNormals[origFaceID]) < 0) {
				vec3i &origFaceIdx = inoutMesh.indices[origFaceID];
				swap(origFaceIdx[1], origFaceIdx[2]);
			}
		}
	}

	return true;
}

bool MeshUtil::reorientAnyMesh(TTriangleMesh &inoutMesh) {

	// no assumption on mesh

	// compute face normals

	vector<vec3> faceNormals(inoutMesh.indices.size());
	for (int faceID = 0; faceID < (int)inoutMesh.indices.size(); faceID++) {
		vec3i faceIdx = inoutMesh.indices[faceID];
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = inoutMesh.positions[faceIdx[k]];
		faceNormals[faceID] = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
	}

	// weld & re-orient mesh

	TTriangleMesh weldedMesh;
	vector<int> weldedFaceIndices;
	if (!weldMeshFaces(inoutMesh, weldedMesh, weldedFaceIndices)) return false;
	if (!reorientMeshFaces(weldedMesh)) return false;

	// compute welded mesh face normals

	vector<vec3> weldedFaceNormals(weldedMesh.indices.size());
	for (int faceID = 0; faceID < (int)weldedMesh.indices.size(); faceID++) {
		vec3i faceIdx = weldedMesh.indices[faceID];
		vec3 faceP[3];
		for (int k = 0; k < 3; k++) faceP[k] = weldedMesh.positions[faceIdx[k]];
		weldedFaceNormals[faceID] = cml::normalize(cml::cross(faceP[1] - faceP[0], faceP[2] - faceP[0]));
	}

	// re-orient original mesh

	for (int faceID = 0; faceID < (int)inoutMesh.indices.size(); faceID++) {
		int weldedFaceID = weldedFaceIndices[faceID];
		if (cml::dot(faceNormals[faceID], weldedFaceNormals[weldedFaceID]) < 0) {
			vec3i &faceIdx = inoutMesh.indices[faceID];
			swap(faceIdx[1], faceIdx[2]);
		}
	}

	return true;
}

bool MeshUtil::reorientPatch(
	TKDTree &inMeshTree,
	TTriangleMesh &inMesh,
	vector<int> &inPatch,
	bool &outFlipFlag,
	double eps)
{

	int numRays = 16;

	if (eps < 0) {
		vec3 bbMin, bbMax;
		if (!MeshUtil::computeAABB(inMesh, bbMin, bbMax)) return false;
		float eps = (bbMax - bbMin).length() * 1e-4f;
	}

	int numPatchFaces = (int)inPatch.size();

	vector<float> faceAreas(numPatchFaces, 0.0);
	vector<bool> faceFlipFlags(numPatchFaces, false);
#pragma omp parallel for
	for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
		int faceID = inPatch[patchFaceID];
		vec3i &faceIdx = inMesh.indices[faceID];
		vec3 facePos[3];
		for (int k = 0; k < 3; k++) facePos[k] = inMesh.positions[faceIdx[k]];
		vec3 faceCenter = (facePos[0] + facePos[1] + facePos[2]) / 3;
		vec3 faceNormal = cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]);
		float area = faceNormal.length();
		if (area) faceNormal /= area;
		faceAreas[patchFaceID] = area;

		int visibleCountP = 0;
		int visibleCountN = 0;
		for (int rayID = 0; rayID < numRays; rayID++) {
			double r1 = cml::random_unit();
			double r2 = cml::random_unit();
			vec3 rayDirP = vec3d(
				2.0 * cos(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
				2.0 * sin(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
				1.0 - 2.0*r2); // random direction on unit sphere
			if (cml::dot(rayDirP, faceNormal) < 0) rayDirP = -rayDirP;
			vec3 rayDirN = -rayDirP;

			vec3 rayOriginP = faceCenter + rayDirP * eps; // add eps for offset
			vec3 rayOriginN = faceCenter + rayDirN * eps;
			Thea::Ray3 rayP(G3D::Vector3(rayOriginP.data()), G3D::Vector3(rayDirP.data()));
			Thea::Ray3 rayN(G3D::Vector3(rayOriginN.data()), G3D::Vector3(rayDirN.data()));
			if (inMeshTree.rayIntersectionTime(rayP) < -0.5) visibleCountP++;
			if (inMeshTree.rayIntersectionTime(rayN) < -0.5) visibleCountN++;
		}
		if (visibleCountN > visibleCountP + 1) { // HACK: if both are similar, prefer not to flip
			faceFlipFlags[patchFaceID] = true;
		}
	}

	double flipWeight = 0;
	double unflipWeight = 0;
	for (int patchFaceID = 0; patchFaceID < numPatchFaces; patchFaceID++) {
		if (faceFlipFlags[patchFaceID]) flipWeight += faceAreas[patchFaceID];
		else unflipWeight += faceAreas[patchFaceID];
	}
	outFlipFlag = flipWeight > unflipWeight;

	return true;
}

bool MeshUtil::recomputeNormals(TTriangleMesh &mesh) {

	int numVertices = mesh.amount;

	vector<vec3d> normals(numVertices, vec3d(0.0,0.0,0.0));
	vector<vec3d> binormals(numVertices, vec3d(0.0,0.0,0.0));
	vector<vec3d> rawnormals(numVertices, vec3d(0.0,0.0,0.0));
	vector<double> weights(numVertices, 0.0);

	for(vec3i faceIdx : mesh.indices) {
		vec3d v[3];
		for(int j=0; j<3; j++) v[j] = mesh.positions[faceIdx[j]];
		vec3d bn[3];
		for (int j = 0; j<3; j++) bn[j] = cml::normalize(v[j] - (v[(j + 1) % 3] + v[(j + 2) % 3]) / 2);

		vec3d vv1 = v[1] - v[0];
		vec3d vv2 = v[2] - v[0];
		double l1 = vv1.length();
		double l2 = vv2.length();
		vec3d n = cml::cross(vv1, vv2);
		double area = n.length() / 2;
		if (n.length_squared()) n.normalize();

		//double w = cml::unsigned_angle(vv1, vv2); // angle weight
		//double w = area; // area weight
		double w = area / cml::sqr(l1*l2); // spherical mesh weight
		
		for(int j=0; j<3; j++) {
			rawnormals[faceIdx[j]] += n;
			normals[faceIdx[j]] += n * w;
			binormals[faceIdx[j]] += bn[j] * w;
			weights[faceIdx[j]] += w;
		}
	}

	bool warnFlag = false;
	mesh.normals.resize(numVertices);
	for(int j=0; j<numVertices; j++) {		
		if (normals[j].length_squared() > 0) {
			mesh.normals[j] = (vec3)cml::normalize(normals[j]);
		} else if (binormals[j].length_squared() > 0) {
			mesh.normals[j] = (vec3)cml::normalize(binormals[j]);
		} else if (rawnormals[j].length_squared() > 0) {
			mesh.normals[j] = (vec3)cml::normalize(rawnormals[j]);
		} else {
			if (!warnFlag) {
				//cout << "Warning: detected zero area faces or unreferenced vertices" << endl;
				warnFlag = true;
			}
			mesh.normals[j] = vec3(0.0f, 1.0f, 0.0f);
		}
	}

	return true;
}

bool MeshUtil::computeDihedralAngle(vec3 center1, vec3 normal1, vec3 center2, vec3 normal2, double &angle) {

	double dihedralAngle = cml::constantsd::pi() - cml::acos_safe((double)cml::dot(normal1, normal2));

	vec3 centerDir = center1 - center2;
	float cosAngle1 = cml::dot(centerDir, normal2); // unnormalized
	float cosAngle2 = cml::dot(-centerDir, normal1); // unnormalized
	bool flag1 = (cosAngle1 >= 0);
	bool flag2 = (cosAngle2 >= 0);

	if (flag1 && flag2) {
		// concave
	} else if (!flag1 && !flag2) {
		// convex
		dihedralAngle = cml::constantsd::two_pi() - dihedralAngle;
	} else {
		if (cml::dot(normal1, normal2) > 0.9f) {
			// prevent numerical error on flat plane
		} else {
			// face orientation incompatible
			dihedralAngle = -cml::constantsd::pi();
		}
	}

	angle = dihedralAngle;

	return true;
}

bool MeshUtil::computeAABB(TTriangleMesh &mesh, vec3 &bbMin, vec3 &bbMax) {

	if (mesh.indices.empty()) {
		// empty mesh
		bbMin = vec3(0.0f, 0.0f, 0.0f);
		bbMax = vec3(0.0f, 0.0f, 0.0f);
		return true;
	}

	bbMin = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	bbMax = vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (vec3i idx : mesh.indices) {
		for (int j = 0; j < 3; j++) {
			bbMin.minimize(mesh.positions[idx[j]]);
			bbMax.maximize(mesh.positions[idx[j]]);
		}
	}

	return true;
}

bool MeshUtil::computeFaceArea(TTriangleMesh &mesh, double &area) {

	area = 0;
	for (vec3i &faceIdx : mesh.indices) {
		vec3 facePos[3];
		for (int k = 0; k < 3; k++) facePos[k] = mesh.positions[faceIdx[k]];
		area += cml::cross(facePos[1] - facePos[0], facePos[2] - facePos[0]).length() * 0.5;
	}

	return true;
}

bool MeshUtil::computeVolume(TTriangleMesh &mesh, double &volume) {

	TTriangleMesh thisMesh = mesh;
	if (!cleanUp(thisMesh)) return false;

	double volumeAABB;
	if (true) {

		// compute AABB volume

		vec3 bbMin, bbMax;
		if (!computeAABB(thisMesh, bbMin, bbMax)) return false;
		vec3 bbSize = bbMax - bbMin;
		volumeAABB = (double)(bbSize[0] * bbSize[1] * bbSize[2]);

		if (volumeAABB == 0) { // empty mesh
			volume = 0;
			return true;
		}
	}

	double volumeOBB;
	if (true) {

		// compute OBB volume

		TSampleSet samples;
		SampleSimplePoissonDisk sspd(&thisMesh);
		if (!sspd.runSampling(100)) return false;
		if (!sspd.exportSample(samples)) return false;

		Eigen::Matrix3Xd sampleMat;
		if (!SampleUtil::buildMatrices(samples, sampleMat)) return false;
		sampleMat = sampleMat.colwise() - sampleMat.rowwise().mean();
		Eigen::JacobiSVD< Eigen::Matrix3Xd > svd(sampleMat, Eigen::ComputeThinU);
		Eigen::Matrix3d transform = svd.matrixU().transpose(); // transformation from global CS to local CS

		Eigen::Matrix3Xd meshMat;
		if (!SampleUtil::buildMatrices(thisMesh, meshMat)) return false;
		meshMat = transform * meshMat;

		Eigen::Vector3d bbMin = meshMat.rowwise().minCoeff();
		Eigen::Vector3d bbMax = meshMat.rowwise().maxCoeff();
		Eigen::Vector3d bbSize = bbMax - bbMin;
		volumeOBB = bbSize[0] * bbSize[1] * bbSize[2];
	}

	volume = min(volumeAABB, volumeOBB); // use smaller one

	return true;
}

bool MeshUtil::buildKdTree(TTriangleMesh &mesh, TKDTree &tree, TKDTreeData &treeData) {

	treeData.resize(mesh.indices.size());

#pragma omp parallel for
	for (int faceID = 0; faceID<(int)mesh.indices.size(); faceID++) {
		vec3i idx = mesh.indices[faceID];
		G3D::Vector3 v0(mesh.positions[idx[0]].data());
		G3D::Vector3 v1(mesh.positions[idx[1]].data());
		G3D::Vector3 v2(mesh.positions[idx[2]].data());
		treeData[faceID].set(TKDT::NamedTriangle(v0, v1, v2, faceID));
	}

	tree.init(treeData.begin(), treeData.end());

	return true;
}

bool MeshUtil::buildKdTree(TTriangleMesh &mesh, SKDTree &tree, SKDTreeData &treeData) {

	vector<vec3> faceCenters(mesh.indices.size());
	for (int faceID = 0; faceID < (int)mesh.indices.size(); faceID++) {
		vec3i idx = mesh.indices[faceID];
		faceCenters[faceID] = (mesh.positions[idx[0]] + mesh.positions[idx[1]] + mesh.positions[idx[2]]) / 3;
	}

	return SampleUtil::buildKdTree(faceCenters, tree, treeData);
}

bool MeshUtil::distanceToMesh(TTriangleMesh &mesh, SKDTree &faceTree, vec3 &point, double &distance) {

	// better to parallelize outside this function

	int numQueryFace = 10;

	SKDT::NamedPoint queryPoint(point[0], point[1], point[2]);
	Thea::BoundedSortedArray<SKDTree::Neighbor> queryResult(numQueryFace);
	faceTree.kClosestElements<Thea::MetricL2>(queryPoint, queryResult);

	int minFace;
	int minType;
	vec3 minPoint;

	double minDistance = DBL_MAX;
	vec3d p = point;
	for (int qID = 0; qID < queryResult.size(); qID++) {
		int faceID = (int)faceTree.getElements()[queryResult[qID].getIndex()].id;
		vec3i idx = mesh.indices[faceID];
		vec3d v[3];
		for (int k = 0; k < 3; k++) v[k] = mesh.positions[idx[k]];
		vec3d n = cml::normalize(cml::cross(v[1] - v[0], v[2] - v[0]));
		double distFace = cml::dot(p - v[0], n);
		vec3d projP = p - n*distFace;
		bool insideFace = true;
		for (int k = 0; k < 3; k++) {
			if (cml::dot(n, cml::cross(projP - v[k], projP - v[(k + 1) % 3])) < 0) {
				insideFace = false;
				break;
			}
		}
		if (insideFace) {
			distFace = fabs(distFace);
			if (distFace < minDistance) {
				minDistance = distFace;
				minFace = faceID;
				minType = 0;
				minPoint = projP;
			}
			
		} else {
			bool insideEdge = false;
			for (int k = 0; k < 3; k++) {
				vec3d &v1 = v[k];
				vec3d &v2 = v[(k + 1) % 3];
				vec3d edge = v2 - v1;
				double sign1 = cml::dot(p - v1, edge);
				double sign2 = cml::dot(p - v2, edge);
				if (sign1 > 0 && sign2 < 0) {
					insideEdge = true;
					vec3 edgeP = v1 + edge * sign1 / edge.length_squared();
					double distEdge = (p - edgeP).length();
					if (distEdge < minDistance) {
						minDistance = distEdge;
						minFace = faceID;
						minType = 1;
						minPoint = edgeP;
					}					
				}
			}
			if (!insideEdge) {
				for (int k = 0; k < 3; k++) {
					double distVertex = (p - v[k]).length();
					if (distVertex < minDistance) {
						minDistance = distVertex;
						minFace = faceID;
						minType = 2;
						minPoint = v[k];
					}
				}
			}
		}
	}

	/*
	if (minDistance > 1e-3f) {
		vector<vec3i> colors(mesh.indices.size(), vec3i(127, 127, 127));
		colors[minFace] = vec3i(255, 0, 0);
		PlyExporter pe;
		if (!pe.addMesh(&mesh.positions, &mesh.normals, &mesh.indices, &colors)) return false;
		if (!pe.output("face.ply")) return false;

		vector<vec3> points;
		points.push_back(point);
		points.push_back(minPoint);
		pe.clearUp();
		if (!pe.addLine(&points)) return false;
		if (!pe.output("line.ply")) return false;

		cout << "Point to mesh distance type: " << minType << endl;
		system("pause");
	}
	*/
	distance = minDistance;

	return true;
}

bool MeshUtil::checkInsideMesh(TKDTree &tree, vec3 &point, bool &insideFlag, double eps) {

	int numQueryRay = 5;

	int insideCount = 0;
	for (int qID = 0; qID < numQueryRay; qID++) {

		vec3d rayOrigin = point;

		double r1 = cml::random_unit();
		double r2 = cml::random_unit();
		vec3d rayDir(
			2.0 * cos(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
			2.0 * sin(cml::constantsd::two_pi()*r1) * cml::sqrt_safe(r2*(1 - r2)),
			1.0 - 2.0*r2); // random direction on unit sphere

		int hitCount = 0;
		while (true) {
			Thea::Ray3 ray(G3D::Vector3(rayOrigin.data()), G3D::Vector3(rayDir.data()));
			double dist = tree.rayIntersectionTime(ray);
			if (dist < 0) break;
			hitCount++;
			rayOrigin += rayDir * (dist + eps);
		}

		insideCount += hitCount % 2;
	}

	insideFlag = (insideCount * 2 > numQueryRay);

	return true;
}
