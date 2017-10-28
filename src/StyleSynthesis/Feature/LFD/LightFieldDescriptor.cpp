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

#include "LightFieldDescriptor.h"

#include <sstream>
#include <omp.h>

#include <Windows.h>
#include <gl/glew.h>
#include <gl/gl.h>

#include "RenderContext.h"
#include "constants.h"

#include "ZernikeDescriptor.h"
#include "FourierDescriptor.h"
#include "EccentricityDescriptor.h"
#include "CircularityDescriptor.h"

#include "PPMWriter.h"

using namespace std;
using namespace LFD;

GLuint texID;
GLuint fboID;
GLuint rbID;

LightFieldDescriptor::LightFieldDescriptor() {

	mpImage = new vector<unsigned char> *[NUM_ANGLE];
	mpCenter = new pair<double, double> *[NUM_ANGLE];
	for (int k = 0; k < NUM_ANGLE; k++) {
		mpImage[k] = new vector<unsigned char>[NUM_CAM];
		mpCenter[k] = new pair<double, double>[NUM_CAM];
	}
}

LightFieldDescriptor::~LightFieldDescriptor() {

	if (mpImage) {
		for (int k = 0; k < NUM_ANGLE; k++) {
			if (mpImage[k]) delete[] mpImage[k];
		}
		delete[] mpImage;
	}

	if (mpCenter) {
		for (int k = 0; k < NUM_ANGLE; k++) {
			if (mpCenter[k]) delete[] mpCenter[k];
		}
		delete[] mpCenter;
	}
}

bool LightFieldDescriptor::init() {

	if (!RenderContext::init()) return false;

	// init flag

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClearDepth(1.0);
	glEnable(GL_DEPTH_TEST);

	// init FBO

	glGenTextures(1, &texID);
	glBindTexture(GL_TEXTURE_2D, texID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	glGenFramebuffers(1, &fboID);
	glBindFramebuffer(GL_FRAMEBUFFER, fboID);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texID, 0);
	glGenRenderbuffers(1, &rbID);
	glBindRenderbuffer(GL_RENDERBUFFER, rbID);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, WINDOW_WIDTH, WINDOW_HEIGHT);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbID);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		cout << "Frame buffer error!" << endl;
		return false;
	}

	// init camera

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-0.5, 0.5, -0.5, 0.5, 0.0, 2.0);
	glViewport(0, 0, (GLsizei)WINDOW_WIDTH, (GLsizei)WINDOW_HEIGHT);

	// init descriptors

	if (!ZernikeDescriptor::init()) return false;

	return true;
}

bool LightFieldDescriptor::finish() {

	glDeleteTextures(1, &texID);
	glDeleteRenderbuffers(1, &rbID);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glDeleteFramebuffers(1, &fboID);

	if (!RenderContext::close()) return false;

	return true;
}

bool LightFieldDescriptor::renderMesh(vector<float> &inVertices, vector<int> &inFaces) {

	if (inFaces.size() % 3 != 0) {
		cout << "Error: invalid size of face vector" << endl;
		return false;
	}
	if (inVertices.size() % 3 != 0) {
		cout << "Error: invalid size of vertex vector" << endl;
		return false;
	}
	int numVertices = (int)inVertices.size() / 3;

	if (!normalize(inVertices, inFaces)) return false;
	for (int angleID = 0; angleID < NUM_ANGLE; angleID++) {
		for (int camID = 0; camID < NUM_CAM; camID++) {
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			gluLookAt(
				CAM_VERTEX[angleID][camID][0],
				CAM_VERTEX[angleID][camID][1],
				CAM_VERTEX[angleID][camID][2],
				0, 0, 0,
				0, 1, 0);
			glPushMatrix();			

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glBegin(GL_TRIANGLES);
			for (int vertexID : inFaces) {
				if (vertexID >= numVertices) {
					cout << "Error: invalid vertex ID in face vector" << endl;
					return false;
				}
				glVertex3f(inVertices[vertexID * 3], inVertices[vertexID * 3 + 1], inVertices[vertexID * 3 + 2]);
			}
			glEnd();
			glPopMatrix();

			mpImage[angleID][camID].resize(WINDOW_HEIGHT*WINDOW_WIDTH);
			glReadPixels(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, &mpImage[angleID][camID].front());

			if (!validateImage(mpImage[angleID][camID])) return false;
			if (!findCenter(mpImage[angleID][camID], mpCenter[angleID][camID])) return false;
		}
	}

	return true;
}

bool LightFieldDescriptor::renderPointCloud(vector<float> &inPoints, float sampleRadius) {

	if (inPoints.size() % 3 != 0) {
		cout << "Error: invalid size of vertex vector" << endl;
		return false;
	}
	int numPoints = (int)inPoints.size() / 3;

	float scale = 1.0f;
	if (!normalize(inPoints, scale)) return false;
	float pixelSize = max(1.0f, 3.0f * (sampleRadius*scale) * min(WINDOW_WIDTH, WINDOW_HEIGHT));

	for (int angleID = 0; angleID < NUM_ANGLE; angleID++) {
		for (int camID = 0; camID < NUM_CAM; camID++) {
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			gluLookAt(
				CAM_VERTEX[angleID][camID][0],
				CAM_VERTEX[angleID][camID][1],
				CAM_VERTEX[angleID][camID][2],
				0, 0, 0,
				0, 1, 0);
			glPushMatrix();

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPointSize(pixelSize);
			glBegin(GL_POINTS);
			for (int pointID = 0; pointID < numPoints; pointID++) {
				glVertex3f(inPoints[pointID * 3], inPoints[pointID * 3 + 1], inPoints[pointID * 3 + 2]);
			}
			glEnd();
			glPopMatrix();

			mpImage[angleID][camID].resize(WINDOW_HEIGHT*WINDOW_WIDTH);
			glReadPixels(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, &mpImage[angleID][camID].front());

			if (!validateImage(mpImage[angleID][camID])) return false;
			if (!findCenter(mpImage[angleID][camID], mpCenter[angleID][camID])) return false;
		}
	}

	return true;
}

bool LightFieldDescriptor::calculate(vector<double> &outDescriptors) {

	// compute descriptors

	vector<double> descriptors[NUM_ANGLE][NUM_CAM];
	for (int angleID = 0; angleID < NUM_ANGLE; angleID++) { // DON'T PARALLELIZE: fftw has global variable

		// radius for Zernike descriptor
		double maxRadius = 0;
		for (int camID = 0; camID < NUM_CAM; camID++) {
			double radius = 0;
			if (!ZernikeDescriptor::radius(mpImage[angleID][camID], mpCenter[angleID][camID], radius)) error("Zernike radius");
			maxRadius = max(maxRadius, radius);
		}

		for (int camID = 0; camID < NUM_CAM; camID++) {

			vector<double> zDescriptor; // Zernike descriptor
			if (!ZernikeDescriptor::calculate(mpImage[angleID][camID], mpCenter[angleID][camID], maxRadius, zDescriptor)) error("Zernike descriptor");

			vector<double> fDescriptor; // Fourier descriptor
			if (!FourierDescriptor::calculate(mpImage[angleID][camID], mpCenter[angleID][camID], fDescriptor)) error("Fourier descriptor");

			double eDescriptor; // Eccentricity descriptor
			if (!EccentricityDescriptor::calculate(mpImage[angleID][camID], mpCenter[angleID][camID], eDescriptor)) error("Eccentricity descriptor");

			double cDescriptor; // Circularity descriptor
			if (!CircularityDescriptor::calculate(mpImage[angleID][camID], cDescriptor)) error("Circularity descriptor");

			vector<double> &desc = descriptors[angleID][camID];
			desc.clear();
			desc.insert(desc.end(), zDescriptor.begin(), zDescriptor.end());
			desc.insert(desc.end(), fDescriptor.begin(), fDescriptor.end());
			desc.push_back(eDescriptor);
			desc.push_back(cDescriptor);
		}
	}

	outDescriptors.clear();
	outDescriptors.reserve(descriptors[0][0].size()*NUM_ANGLE*NUM_CAM);
	for (int angleID = 0; angleID < NUM_ANGLE; angleID++) {
		for (int camID = 0; camID < NUM_CAM; camID++) {
			vector<double> &desc = descriptors[angleID][camID];
			outDescriptors.insert(outDescriptors.end(), desc.begin(), desc.end());
		}
	}

	return true;
}

bool LightFieldDescriptor::visualize(string fileName) {

	vector<unsigned char> fullImage(WINDOW_WIDTH*WINDOW_HEIGHT*NUM_ANGLE*NUM_CAM);
	for (int angleID = 0; angleID < NUM_ANGLE; angleID++) {
		for (int camID = 0; camID < NUM_CAM; camID++) {
			vector<unsigned char> &image = mpImage[angleID][camID];
#pragma omp parallel for
			for (int row = 0; row < WINDOW_HEIGHT; row++) {
				copy(image.begin() + WINDOW_WIDTH*row,
					image.begin() + WINDOW_WIDTH*(row + 1),
					fullImage.begin() + camID*WINDOW_WIDTH
					+ (angleID*WINDOW_HEIGHT + (WINDOW_HEIGHT - row - 1))
					* WINDOW_WIDTH*NUM_CAM);
			}
		}
	}
	if (!PPMWriter::output(fileName, fullImage, WINDOW_WIDTH*NUM_CAM, WINDOW_HEIGHT*NUM_ANGLE)) return false;

	return true;
}

bool LightFieldDescriptor::compare(vector<double> &inDescriptor1, vector<double> &inDescriptor2, vector<double> &outDistance) {

	// constants...

	vector<int> descSizes;
	descSizes.push_back(ZernikeDescriptor::DESCRIPTOR_SIZE);
	descSizes.push_back(FourierDescriptor::DESCRIPTOR_SIZE);
	descSizes.push_back(1); // eccentricity
	descSizes.push_back(1); // circularity

	vector<double> descWeights; // weights for each descriptor for comparison (following original code)
	descWeights.push_back(1.0); // Zernike
	descWeights.push_back(2.0); // Fourier
	descWeights.push_back(1.0); // Eccentricity
	descWeights.push_back(1.0); // Circularity

	int allDescSize = 0;
	for (int sz : descSizes) allDescSize += sz;	
	int totalSize = NUM_ANGLE * NUM_CAM * allDescSize;

	if (totalSize != (int)inDescriptor1.size() || totalSize != (int)inDescriptor2.size()) {
		cout << "Error: incorrect descriptor size" << endl;
		return false;
	}

	// init omp

	struct TLFD {
		int srcAngle;
		int dstAngle;
		int rotID;
	};

#ifdef _OPENMP
	int numThreads = omp_get_max_threads();
#else
	int numThreads = 1;
#endif

	vector<double> ompMinError(numThreads, DBL_MAX);
	vector<TLFD> ompMinRecord(numThreads);

	// find best matched rotation

#pragma omp parallel for num_threads(numThreads)
	for (int angleID = 0; angleID < NUM_ANGLE*NUM_ANGLE; angleID++) {

#ifdef _OPENMP
		int threadID = omp_get_thread_num();
#else
		int threadID = 0;
#endif

		int srcAngle = angleID / NUM_ANGLE;
		int dstAngle = angleID % NUM_ANGLE;

		// compute cost for all pairs of camera positions (not all will be used though...)
		double cost[NUM_CAM][NUM_CAM];

		for (int srcCam = 0; srcCam < NUM_CAM; srcCam++) {
			int srcOffset = allDescSize * (NUM_CAM * srcAngle + srcCam);

			for (int dstCam = 0; dstCam < NUM_CAM; dstCam++) {				
				int dstOffset = allDescSize * (NUM_CAM * dstAngle + dstCam);

				double dist = 0;
				for (int k = 0; k < allDescSize; k++) {
					dist += fabs(inDescriptor1[srcOffset + k] - inDescriptor2[dstOffset + k]); // L1 distance
				}
				cost[srcCam][dstCam] = dist;
			}
		}

		// find minimum cost rotation
		for (int rotID = 0; rotID < NUM_ROTATION; rotID++) {

			double error = 0;
			for (int camID = 0; camID < NUM_CAM; camID++) {
				error += cost[ CAM_ORDER[rotID][camID] ][ CAM_ORDER[0][camID] ];
			}

			if (error < ompMinError[threadID]) {
				ompMinError[threadID] = error;
				ompMinRecord[threadID].srcAngle = srcAngle;
				ompMinRecord[threadID].dstAngle = dstAngle;
				ompMinRecord[threadID].rotID = rotID;				
			}
		}
	}

	// get mininum between threads

	double minError = DBL_MAX;
	int minThreadID = -1;
	for (int threadID = 0; threadID < numThreads; threadID++) {
		if (ompMinError[threadID] < minError) {
			minError = ompMinError[threadID];
			minThreadID = threadID;
		}
	}

	// output distances (seperated for different descriptors)

	outDistance.clear();
	outDistance.resize(descSizes.size(), 0);
	TLFD &minRecord =  ompMinRecord[minThreadID];
	for (int camID = 0; camID < NUM_CAM; camID++) {
		int srcCam = CAM_ORDER[minRecord.rotID][camID];
		int dstCam = CAM_ORDER[0][camID];
		int srcOffset = allDescSize * (NUM_CAM * minRecord.srcAngle + srcCam);
		int dstOffset = allDescSize * (NUM_CAM * minRecord.dstAngle + dstCam);

		int k = 0;
		for (int descID = 0; descID < (int)descSizes.size(); descID++) {
			double dist = 0;
			for (int dim = 0; dim < descSizes[descID]; dim++) {
				dist += fabs(inDescriptor1[srcOffset + k] - inDescriptor2[dstOffset + k]) * descWeights[descID];
				k++;
			}
			outDistance[descID] += dist;
		}
	}

	/* // for debug : output matching result
	cout << "srcAngle = " << minRecord.srcAngle << endl;
	cout << "dstAngle = " << minRecord.dstAngle << endl;
	cout << "rotID = " << minRecord.rotID << endl;

	cout << "src camera pos: ";
	for (double v : CAM_VERTEX[minRecord.srcAngle][CAM_ORDER[minRecord.rotID][0]]) cout << v << " ";
	cout << endl;

	cout << "dst camera pos: ";
	for (double v : CAM_VERTEX[minRecord.dstAngle][CAM_ORDER[0][0]]) cout << v << " ";
	cout << endl;
	*/

	return true;
}

bool LightFieldDescriptor::normalize(vector<float> &vertices, vector<int> &faces) {

	if (faces.size() % 3 != 0) {
		cout << "Error: invalid size of face vector" << endl;
		return false;
	}
	if (vertices.size() % 3 != 0) {
		cout << "Error: invalid size of vertex vector" << endl;
		return false;
	}
	int numVertices = (int)vertices.size() / 3;

	float minBB[3], maxBB[3];
	for (int j = 0; j < 3; j ++ ) {
		minBB[j] = FLT_MAX;
		maxBB[j] = -FLT_MAX;
	}
	for (int vertexID : faces) {
		if (vertexID >= numVertices) {
			cout << "Error: invalid vertex ID in face vector" << endl;
			return false;
		}
		for (int j = 0; j < 3; j++) {
			minBB[j] = min(minBB[j], vertices[vertexID * 3 + j]);
			maxBB[j] = max(maxBB[j], vertices[vertexID * 3 + j]);
		}
	}
	float translate[3];
	float scale = 0;
	for (int j = 0; j < 3; j++) {
		translate[j] = -(minBB[j] + maxBB[j]) / 2;
		scale += (maxBB[j] - minBB[j]) * (maxBB[j] - minBB[j]);
	}
	if (scale>0) scale = 1 / sqrt(scale);
	else scale = 1;

	for (int vertexID = 0; vertexID < numVertices; vertexID++) {
		for (int j = 0; j < 3; j++) {
			vertices[vertexID * 3 + j] = (vertices[vertexID * 3 + j] + translate[j]) * scale;
		}
	}

	return true;
}

bool LightFieldDescriptor::normalize(vector<float> &points, float &scale) {

	if (points.size() % 3 != 0) {
		cout << "Error: invalid size of vertex vector" << endl;
		return false;
	}
	int numPoints = (int)points.size() / 3;

	float minBB[3], maxBB[3];
	for (int j = 0; j < 3; j++) {
		minBB[j] = FLT_MAX;
		maxBB[j] = -FLT_MAX;
	}
	for (int pointID = 0; pointID < numPoints; pointID++) {
		for (int j = 0; j < 3; j++) {
			minBB[j] = min(minBB[j], points[pointID * 3 + j]);
			maxBB[j] = max(maxBB[j], points[pointID * 3 + j]);
		}
	}
	float translate[3];
	for (int j = 0; j < 3; j++) {
		translate[j] = -(minBB[j] + maxBB[j]) / 2;
		scale += (maxBB[j] - minBB[j]) * (maxBB[j] - minBB[j]);
	}
	if (scale>0) scale = 1 / sqrt(scale);
	else scale = 1;

	for (int pointID = 0; pointID < numPoints; pointID++) {
		for (int j = 0; j < 3; j++) {
			points[pointID * 3 + j] = (points[pointID * 3 + j] + translate[j]) * scale;
		}
	}

	return true;
}

bool LightFieldDescriptor::validateImage(vector<unsigned char> &image) {

	int count = 0;
	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			int pos = y*WINDOW_WIDTH + x;
			if (image[pos] < 255) count++;
		}
	}

	if (count == 0) {
		// invisible - make it full black
		for (auto &pixel : image) pixel = 0;
	}

	return true;
}

bool LightFieldDescriptor::findCenter(vector<unsigned char> &inImage, pair<double, double> &outCenter) {

	double minX = DBL_MAX;
	double maxX = -DBL_MAX;
	double minY = DBL_MAX;
	double maxY = -DBL_MAX;
	for (int y = 0; y < WINDOW_HEIGHT; y++) {
		for (int x = 0; x < WINDOW_WIDTH; x++) {
			int pos = y*WINDOW_WIDTH + x;
			if (inImage[pos] < 255) {
				minX = min(minX, x);
				maxX = max(maxX, x);
				minY = min(minY, y);
				maxY = max(maxY, y);
			}
		}
	}

	if (minX <= maxX && minY <= maxY) {
		outCenter.first = (minX + maxX) / 2;
		outCenter.second = (minY + maxY) / 2;
	} else {
		outCenter.first = -1;
		outCenter.second = -1;
	}

	return true;
}