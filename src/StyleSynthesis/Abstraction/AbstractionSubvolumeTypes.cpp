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

#include "AbstractionSubvolumeTypes.h"

#include <Eigen/Eigen>

#include "Mesh/MeshUtil.h"

#include "Segment/SegmentUtil.h"
#include "Utility/PlyExporter.h"

using namespace StyleSynthesis;
using namespace AbstractionSubvolumeNamespace;

#define DEBUG_VISUALIZATION

bool TRecPrism::computeArea() {
	
	float e1 = (mCorner[1] - mCorner[0]).length();
	float e2 = (mCorner[2] - mCorner[0]).length();
	float e3 = (mCorner[3] - mCorner[0]).length();
	mArea = (e1*e2 + e1*e3 + e2*e3) * 2;

	return true;
}

bool TTriPrism::computeArea() {

	float h = (mCorner[1] - mCorner[0]).length();
	float e1 = (mCorner[2] - mCorner[0]).length();
	float e2 = (mCorner[3] - mCorner[0]).length();
	float e3 = (mCorner[2] - mCorner[3]).length();
	mArea = cml::cross(mCorner[2] - mCorner[0], mCorner[3] - mCorner[0]).length() + (e1 + e2 + e3)*h;

	return true;
}

bool TCylinder::computeArea() {

	float h = (mCenter[1] - mCenter[0]).length();
	float r = mRadius;
	mArea = cml::constantsf::two_pi()*r*(h + r);

	return true;
}

bool TTruncCone::computeArea() {

	float d1 = (mCenter[1] - mCenter[0]).length();
	float d2 = (mCenter[2] - mCenter[0]).length();
	float tangent = tan(mAngle);
	float r1 = d1*tangent;
	float r2 = d2*tangent;
	mArea = cml::constantsf::pi()*((r1 + r2)*(d2 - d1) + r1*r1 + r2*r2);

	return true;
}

bool TRecPrism::fit(TPointSet *points, vec4i samples) {

	auto &p = mCorner; // shorter name

	TTriPrism triP; // triangular prism
	if (!triP.fit(points, samples)) return false;
	for (int k = 0; k < 4; k++) {
		p[k] = triP.mCorner[k];
	}

	float eps = 0;
	for (int k = 1; k < 4; k++) {
		eps = max(eps, (p[k] - p[0]).length());
	}
	eps *= 1e-3f;

	// orthogonalize base
	vec3 newAxis = cml::cross(p[1] - p[0], p[2] - p[0]);
	float l1 = newAxis.length();
	if (l1 < eps) return false;
	newAxis /= l1;
	vec3 oldAxis = p[3] - p[0];
	//float l2 = cml::dot(oldAxis, newAxis); // projected length
	float l2 = oldAxis.length(); // keep original length
	if (l2 < eps) return false;
	p[3] = p[0] + newAxis * l2;

	if (!computeArea()) return false;

	return true;
}

bool TTriPrism::fit(TPointSet *points, vec4i samples) {

	auto &p = mCorner; // shorter name
	for (int k = 0; k < 4; k++) {
		p[k] = points->positions[samples[k]];
	}

	float eps = 0;
	for (int k = 1; k < 4; k++) {
		eps = max(eps, (p[k] - p[0]).length());
	}
	eps *= 1e-3f;

	if (!(fabs(cml::dot(cml::normalize(p[1] - p[0]), cml::normalize(cml::cross(p[2] - p[0], p[3] - p[0])))) > eps)) return false; // co-planar

	// project to base plane
	vec3 axis = cml::normalize(p[1] - p[0]);
	for (int k = 2; k < 4; k++) {
		p[k] = p[k] - axis*cml::dot(p[k] - p[0], axis);
	}
	if (cml::dot(cml::cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0]) < 0) swap(p[2], p[3]); // rhs cs

	if (!computeArea()) return false;

	return true;
}

bool TCylinder::fit(TPointSet *points, vec4i samples) {

	if (samples[0] == samples[1]) return false;

	vec3 v1 = points->positions[samples[0]];
	vec3 v2 = points->positions[samples[1]];
	vec3 n1 = points->normals[samples[0]];
	vec3 n2 = points->normals[samples[1]];

	float eps = (v1 - v2).length() * 1e-3f;

	vec3 axis = cml::cross(n1, n2);
	float al = axis.length();
	if (al < eps) return false;                 //    ^
	axis /= al;                                 //    | n1
	vec3 vp1 = v1 - axis * cml::dot(v1, axis);  //    |
	vec3 vp2 = v2 - axis * cml::dot(v2, axis);  //    |
	vec3 cnv = cml::cross(n2, vp1 - vp2);       //    o vp1	                                                
	float t = cnv.length() / al;                //        vp2  n2
	if (cml::dot(cnv, axis) < 0) t = -t;        //    o   o-------->
	vec3 vo = vp1 + n1 * t;                     //  vo

	if (fabs(cml::dot(cml::normalize(v1 - v2), axis)) < eps) return false;

	mCenter[0] = vo + axis*cml::dot(v1 - vo, axis);
	mCenter[1] = vo + axis*cml::dot(v2 - vo, axis);
	mRadius = ((vp1 - vo).length() + (vp2 - vo).length()) / 2;

	if (!computeArea()) return false;

	return true;
}

bool TTruncCone::fit(TPointSet *points, vec4i samples) {

	if (samples[0] == samples[1] || samples[0] == samples[2] || samples[1] == samples[2]) return false;

	vec3 v[3], n[3];
	for (int k = 0; k < 3; k++) {
		v[k] = points->positions[samples[k]];
		n[k] = points->normals[samples[k]];
	}

	// tip point position (intersection point of 3 planes)
	// ref: http://mathworld.wolfram.com/Plane-PlaneIntersection.html

	float eps = 0;
	for (int k = 0; k < 3; k++) {
		eps = max(eps, (v[k] - v[(k + 1) % 3]).length());
	}
	eps *= 1e-3f;
	float angleEps = 1e-3f;

	vec3 tip = cml::dot(v[0], n[0])*cml::cross(n[1], n[2]) +
		cml::dot(v[1], n[1])*cml::cross(n[2], n[0]) +
		cml::dot(v[2], n[2])*cml::cross(n[0], n[1]);
	Eigen::Matrix3f mat;
	mat << n[0][0], n[1][0], n[2][0],
		n[0][1], n[1][1], n[2][1],
		n[0][2], n[1][2], n[2][2];
	float denom = mat.determinant();
	if (denom < eps) return false;
	tip /= denom;

	// cone axis (normal of the plane formed by 3 normalized point vectors)
	vec3 vt[3]; // normalized point vector
	for (int k = 0; k < 3; k++) {
		vt[k] = cml::normalize(v[k] - tip) + tip;
	}
	vec3 axis = cml::cross(vt[1] - vt[0], vt[2] - vt[0]);
	float axisLen = axis.length();
	if (axisLen < eps) return false;
	axis /= axisLen;
	if (cml::dot(n[0], axis) + cml::dot(n[1], axis) + cml::dot(n[2], axis) > 0) axis = -axis;
	float angle = cml::acos_safe(cml::dot(axis, vt[0] - tip));
	if (angle < angleEps) return false;
	if (fabs(cml::constantsf::pi_over_2() - angle) < angleEps) return false;

	vec3 apMin, apMax;
	float alMin = FLT_MAX;
	float alMax = -FLT_MAX;
	for (int k = 0; k < 3; k++) {
		float al = cml::dot(v[k] - tip, axis); // distance to tip point along axis direction
		vec3 ap = tip + axis*al; // projected point on axis
		if (al < alMin) {
			alMin = al;
			apMin = ap;
		}
		if (al > alMax) {
			alMax = al;
			apMax = ap;
		}
	}
	if (alMin < 0) return false;

	mCenter[0] = tip;
	mCenter[1] = apMin;
	mCenter[2] = apMax;
	mAngle = angle;

	if (!computeArea()) return false;

	return true;
}

bool TRecPrism::refit(TPointSet &points, int index) {

	if (index >= 6) return false;

	// expand the prism along axis

	if(true) {
		int axisID = (index % 3) + 1;
		vec3 axis = mCorner[axisID] - mCorner[0];
		float axisLen = axis.length();
		axis /= axisLen;
		float minLen = 0;
		float maxLen = axisLen;
		for (vec3 &p : points.positions) {
			float projLen = cml::dot(p - mCorner[0], axis);
			minLen = min(minLen, projLen);
			maxLen = max(maxLen, projLen);
		}
		if (index / 3) {
			mCorner[axisID] += (maxLen - axisLen)*axis;
		} else {
			for (int j = 0; j < 4; j++) {
				if (j != axisID) mCorner[j] += minLen*axis;
			}
		}
	}

	float oldArea = mArea;
	if (!computeArea() || oldArea >= mArea) return false;

	return true;
}

bool TTriPrism::refit(TPointSet &points, int index) {

	if (index >= 6) return false;

	// expand the prism to wrap all points

	if (true) {
		int axisID = (index % 3) + 1;
		vec3 axis = mCorner[axisID] - mCorner[0];
		float axisLen = axis.length();
		axis /= axisLen;
		float minLen = 0;
		float maxLen = axisLen;
		for (vec3 &p : points.positions) {
			float projLen = cml::dot(p - mCorner[0], axis);
			minLen = min(minLen, projLen);
			maxLen = max(maxLen, projLen);
		}
		if (index / 3) {
			mCorner[axisID] += (maxLen - axisLen)*axis;
		} else {
			for (int j = 0; j < 4; j++) {
				if (j != axisID) mCorner[j] += minLen*axis;
			}
		}
	}

	float oldArea = mArea;
	if (!computeArea() || oldArea <= mArea) return false;

	return true;
}

bool TCylinder::refit(TPointSet &points, int index) {

	if (index >= 2) return false;

	// expand top/bottom to wrap all points

	vec3 axis = mCenter[1] - mCenter[0];
	float axisLen = axis.length();
	axis /= axisLen;
	float minLen = 0;
	float maxLen = axisLen;
	for (vec3 &p : points.positions) {
		float projLen = cml::dot(p - mCenter[0], axis);
		minLen = min(minLen, projLen);
		maxLen = max(maxLen, projLen);
	}
	if (index) {
		mCenter[1] += (maxLen - axisLen)*axis;
	} else {
		mCenter[0] += minLen*axis;
	}

	float oldArea = mArea;
	if (!computeArea() || oldArea <= mArea) return false;

	return true;
}

bool TTruncCone::refit(TPointSet &points, int index) {

	if (index >= 2) return false;

	// expand top/bottom to wrap all points

	vec3 axis = mCenter[2] - mCenter[0];
	float axisLen = axis.length();
	axis /= axisLen;
	float minLen = (mCenter[1]-mCenter[0]).length();
	float maxLen = axisLen;
	for (vec3 &p : points.positions) {
		float projLen = cml::dot(p - mCenter[0], axis);
		minLen = min(minLen, projLen);
		maxLen = max(maxLen, projLen);
	}
	if (index) {
		mCenter[1] = mCenter[0] + min(minLen, axisLen)*axis;
	} else {
		mCenter[2] = mCenter[0] + max(maxLen, axisLen)*axis;
	}

	float oldArea = mArea;
	if (!computeArea() || oldArea <= mArea) return false;

	return true;
}

float TRecPrism::distance(vec3 point) {

	float minDist = FLT_MAX;
	vec3 basis[3];
	float basisLen[3];
	for (int k = 0; k < 3; k++) {
		basis[k] = mCorner[k + 1] - mCorner[0];
		basisLen[k] = basis[k].length();
		basis[k] /= basisLen[k];
	}
	vec3 &o = mCorner[0];
	for (int k = 0; k < 3; k++) {
		// use one axis as normal (n) and two other axes as basis vectors (d1, d2)
		vec3 &n = basis[k];
		vec3 &d1 = basis[(k + 1) % 3];
		vec3 &d2 = basis[(k + 2) % 3];
		float ln = basisLen[k];
		float l1 = basisLen[(k + 1) % 3];
		float l2 = basisLen[(k + 2) % 3];
		float dist = cml::dot(point - o, n); // distance along n
		vec3 projPoint = point - n*dist;
		float dp1 = cml::dot(projPoint - o, d1); // projected distance along d1
		float dp2 = cml::dot(projPoint - o, d2); // projected distance along d2
		if (dp1 >= 0 && dp1 <= l1 && dp2 >= 0 && dp2 <= l2) { // projected point within surface
			dist = min(fabs(dist), fabs(ln - dist)); // min distance to the two planes
			minDist = min(minDist, dist);
		}
	}

	return minDist;
}

float TTriPrism::distance(vec3 point) {

	vec3 n = mCorner[1] - mCorner[0];
	float ln = n.length();
	n /= ln;
	vec3 c[3];
	c[0] = mCorner[0];
	c[1] = mCorner[2];
	c[2] = mCorner[3];

	float minDist = FLT_MAX;

	float dn = cml::dot(point - c[0], n); // distance to base plane
	if (dn >= 0 && dn <= ln) {
		// check distance to side surfaces
		for (int k = 0; k < 3; k++) {
			// use one base edge as basis vector (e) and the height vector (n) as another basis vector
			vec3 e = c[(k + 1) % 3] - c[k];
			float le = e.length();
			e /= le;
			vec3 pn = cml::cross(e, n); // plane normal (should already be normalized)
			float dist = cml::dot(point - c[k], pn); // distance to plane
			vec3 projPoint = point - pn*dist; // project to side surface
			float de = cml::dot(projPoint - c[k], e); // distance of projected point to edge
			if (de >= 0 && de <= le) { // projected point within surface
				minDist = min(minDist, fabs(dist));
			}
		}
	}

	if (true) {
		// check distance to top/bottom surfaces
		vec3 p = point - n*dn; // project to base surface
		bool inside = true;
		for (int k = 0; k < 3; k++) {
			if (cml::dot(n, cml::cross(p - c[k], p - c[(k + 1) % 3])) < 0) {
				inside = false;
				break;
			}
		}
		if (inside) { // within base surface
			float dist = min(fabs(dn), fabs(ln - dn)); // min distance to the two planes
			minDist = min(minDist, dist);
		}
	}

	return minDist;
}

float TCylinder::distance(vec3 point) {

	vec3 o = mCenter[0];
	vec3 n = mCenter[1] - o;
	float ln = n.length();
	n /= ln;

	float minDist = FLT_MAX;

	float dn = cml::dot(point - o, n); // distance to bottom plane
	vec3 axisPoint = o + n*dn; // project to axis
	float da = (point - axisPoint).length(); // distance to axis

	if (dn >= 0 && dn <= ln) {
		// check distance to side surface
		float dist = fabs(mRadius - da);
		minDist = min(minDist, dist);
	}
	if (da <= mRadius) {
		// check distance to top/bottom surface
		float dist = min(fabs(dn), fabs(ln - dn)); // min distance to the two planes
		minDist = min(minDist, dist);
	}

	return minDist;
}

float TTruncCone::distance(vec3 point) {

	vec3 o = mCenter[0]; // tip point
	vec3 d1 = mCenter[1] - o;
	vec3 d2 = mCenter[2] - o;
	float l1 = d1.length();
	float l2 = d2.length();
	vec3 axis = d2 / l2;

	float minDist = FLT_MAX;

	float dp = cml::dot(point - o, axis); // distance to tip plane
	if (dp < 0) return minDist; // on other side of the cone

	float dt = (point - o).length(); // distance to tip point

	float cosA = cos(mAngle);
	float tanA = tan(mAngle);

	if (dp >= l1 && dp <= l2) {
		// check distance to side surface
		float dist = fabs((o + axis*dp - point).length() - tanA*dp) * cosA;
		minDist = min(minDist, dist);
	}

	if (dp >= dt*cosA) {
		// check distance to top/bottom plane
		float dist = min(fabs(l1 - dp), fabs(l2 - dp)); // min distance to the two planes
		minDist = min(minDist, dist);
	}

	return minDist;
}

bool TRecPrism::inside(vec3 point) {

	for (int k = 0; k < 3; k++) {
		vec3 axis = mCorner[k + 1] - mCorner[0];
		float axisLen = axis.length();
		axis /= axisLen;
		float pointLen = cml::dot(point - mCorner[0], axis);
		if (pointLen < 0 || pointLen > axisLen) return false;
	}

	return true;
}

bool TTriPrism::inside(vec3 point) {

	vec3 normal = mCorner[1] - mCorner[0];
	float normalLen = normal.length();
	normal /= normalLen;
	float pointLen = cml::dot(point - mCorner[0], normal);
	if (pointLen < 0 || pointLen > normalLen) return false; // not within top/bottom surface

	vec3 projPoint = point - normal*pointLen;
	vec3 projVec1 = projPoint - mCorner[0]; // projected vector on bottom triangle
	vec3 projVec2 = projPoint - mCorner[2];
	vec3 projVec3 = projPoint - mCorner[3];
	if (cml::dot(normal, cml::cross(projVec1, projVec2)) < 0) return false; // not within bottom triangle
	if (cml::dot(normal, cml::cross(projVec2, projVec3)) < 0) return false;
	if (cml::dot(normal, cml::cross(projVec3, projVec1)) < 0) return false;

	return true;
}

bool TCylinder::inside(vec3 point) {

	vec3 axis = mCenter[1] - mCenter[0];
	float axisLen = axis.length();
	axis /= axisLen;
	float pointLen = cml::dot(point - mCenter[0], axis);
	if (pointLen < 0 || pointLen > axisLen) return false; // not within top/bottom surface

	vec3 projPoint = point - axis*pointLen;
	if ((projPoint - mCenter[0]).length() > mRadius) return false; // not within cylinder

	return true;
}

bool TTruncCone::inside(vec3 point) {

	vec3 o = mCenter[0]; // tip point
	vec3 d1 = mCenter[1] - o;
	vec3 d2 = mCenter[2] - o;
	float l1 = d1.length();
	float l2 = d2.length();
	vec3 axis = d2 / l2;

	float pointLen = cml::dot(point - mCenter[0], axis);
	if (pointLen < l1 || pointLen > l2) return false; // not within top/bottom surface

	float dp = cml::dot(point - o, axis); // distance to tip plane
	if (dp < 0) return false; // on other side of the cone

	float dt = (point - o).length(); // distance to tip point
	float aa = cml::acos_safe(dp / dt); // angle to axis
	if (aa > mAngle) return false; // not within cone

	return true;
}

bool TRecPrism::tesselate(TTriangleMesh &mesh) {

	auto &p = mCorner; // shorter name

	mesh.positions.clear();
	for (int k = 0; k < 4; k++) {
		mesh.positions.push_back(p[k]);
	}
	for (int k = 0; k < 3; k++) {
		mesh.positions.push_back(p[k + 1] + p[(k + 1) % 3 + 1] - p[0]);
	}
	mesh.positions.push_back(p[1] + p[2] + p[3] - p[0] * 2);
	mesh.amount = (int)mesh.positions.size();

	mesh.indices.clear();
	mesh.indices.push_back(vec3i(0, 2, 1)); mesh.indices.push_back(vec3i(1, 2, 4));
	mesh.indices.push_back(vec3i(0, 3, 2)); mesh.indices.push_back(vec3i(2, 3, 5));
	mesh.indices.push_back(vec3i(0, 1, 3)); mesh.indices.push_back(vec3i(3, 1, 6));
	mesh.indices.push_back(vec3i(1, 4, 6)); mesh.indices.push_back(vec3i(6, 4, 7));
	mesh.indices.push_back(vec3i(2, 5, 4)); mesh.indices.push_back(vec3i(4, 5, 7));
	mesh.indices.push_back(vec3i(3, 6, 5)); mesh.indices.push_back(vec3i(5, 6, 7));

	if (!MeshUtil::recomputeNormals(mesh)) return false;

	return true;
}

bool TTriPrism::tesselate(TTriangleMesh &mesh) {

	auto &p = mCorner; // shorter name

	mesh.positions.clear();
	for (int k = 0; k < 4; k++) {
		mesh.positions.push_back(p[k]);
	}
	mesh.positions.push_back(p[1] + p[2] - p[0]);
	mesh.positions.push_back(p[1] + p[3] - p[0]);
	mesh.amount = (int)mesh.positions.size();

	mesh.indices.push_back(vec3i(0, 2, 1)); mesh.indices.push_back(vec3i(1, 2, 4));
	mesh.indices.push_back(vec3i(0, 1, 3)); mesh.indices.push_back(vec3i(3, 1, 5));
	mesh.indices.push_back(vec3i(2, 3, 4)); mesh.indices.push_back(vec3i(4, 3, 5));
	mesh.indices.push_back(vec3i(0, 3, 2)); mesh.indices.push_back(vec3i(1, 4, 5));

	if (!MeshUtil::recomputeNormals(mesh)) return false;

	return true;
}

bool TCylinder::tesselate(TTriangleMesh &mesh) {

	const int subdiv = 20;

	auto &p = mCenter; // shorter name

	vec3 n = cml::normalize(p[1] - p[0]);
	vec3 d1 = cml::normalize(cml::cross(n, fabs(n[0]) < 0.5f ? vec3(1.0f, 0.0f, 0.0f) : vec3(0.0f, 1.0f, 0.0f)));
	vec3 d2 = cml::cross(n, d1); // local CS: (n, d1, d2)

	mesh.positions.resize(subdiv * 2 + 2);
	for (int k = 0; k < subdiv; k++) {
		float angle = cml::constantsf::two_pi() * k / subdiv;
		vec3 vec = (cos(angle) * d1 + sin(angle)*d2) * mRadius; // swiping vector along circle
		mesh.positions[k * 2] = p[0] + vec;
		mesh.positions[k * 2 + 1] = p[1] + vec;
	}
	mesh.positions[subdiv * 2] = p[0];
	mesh.positions[subdiv * 2 + 1] = p[1];
	mesh.amount = (int)mesh.positions.size();

	mesh.indices.clear();
	int num = subdiv * 2;
	for (int k = 0; k < num; k++) {
		if (k % 2) mesh.indices.push_back(vec3i(k, (k + 1) % num, (k + 2) % num));
		else mesh.indices.push_back(vec3i(k, (k + 2) % num, (k + 1) % num));
	}
	for (int k = 0; k < num; k++) {
		if (k % 2) mesh.indices.push_back(vec3i(num + 1, k, (k + 2) % num));
		else mesh.indices.push_back(vec3i(num, (k + 2) % num, k));
	}

	if (!MeshUtil::recomputeNormals(mesh)) return false;

	return true;
}

bool TTruncCone::tesselate(TTriangleMesh &mesh) {

	const int subdiv = 20;

	auto &p = mCenter; // shorter name

	vec3 n = cml::normalize(p[2] - p[0]);
	vec3 d1 = cml::normalize(cml::cross(n, fabs(n[0]) < 0.5f ? vec3(1.0f, 0.0f, 0.0f) : vec3(0.0f, 1.0f, 0.0f)));
	vec3 d2 = cml::cross(n, d1); // local CS: (n, d1, d2)

	float tangent = tan(mAngle);
	float radius1 = (p[1] - p[0]).length() * tangent;
	float radius2 = (p[2] - p[0]).length() * tangent;

	mesh.positions.resize(subdiv * 2 + 2);
	for (int k = 0; k < subdiv; k++) {
		float angle = cml::constantsf::two_pi() * k / subdiv;
		vec3 vec = (cos(angle) * d1 + sin(angle)*d2); // swiping vector along circle
		mesh.positions[k * 2] = p[1] + vec * radius1;
		mesh.positions[k * 2 + 1] = p[2] + vec * radius2;
	}
	mesh.positions[subdiv * 2] = p[1];
	mesh.positions[subdiv * 2 + 1] = p[2];
	mesh.amount = (int)mesh.positions.size();

	mesh.indices.clear();
	int num = subdiv * 2;
	for (int k = 0; k < num; k++) {
		if (k % 2) mesh.indices.push_back(vec3i(k, (k + 1) % num, (k + 2) % num));
		else mesh.indices.push_back(vec3i(k, (k + 2) % num, (k + 1) % num));
	}
	for (int k = 0; k < num; k++) {
		if (k % 2) mesh.indices.push_back(vec3i(num + 1, k, (k + 2) % num));
		else mesh.indices.push_back(vec3i(num, (k + 2) % num, k));
	}

	if (!MeshUtil::recomputeNormals(mesh)) return false;

	return true;
}
