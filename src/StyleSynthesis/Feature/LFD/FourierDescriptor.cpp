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

#include "FourierDescriptor.h"

#include <iostream>
#include <algorithm>

#include "fftw3.h"

#include "constants.h"

#include "PPMWriter.h"

using namespace LFD;

bool FourierDescriptor::calculate(
	vector<unsigned char> &inImage,
	pair<double, double> &inCenter,
	vector<double> &outDescriptor)
{
	vector<int> contourPoints;

	if (!connect(inImage)) return false;
	if (!extractContour(inImage, contourPoints)) return false;
	if (!getDescriptor(contourPoints, inCenter, outDescriptor)) return false;

	return true;
}

bool FourierDescriptor::extractContour(
	vector<unsigned char> &inImage,
	vector<int> &outContour)
{
	// contour chaining

	int startPos = 0;
	while (startPos < WINDOW_HEIGHT*WINDOW_WIDTH) {
		if (inImage[startPos] < 255) break;
		startPos++;
	}
	if (startPos == WINDOW_HEIGHT*WINDOW_WIDTH) { // empty image
		outContour.clear();
		return true;
	}

	const int nextPos[8] = {
		WINDOW_WIDTH - 1, WINDOW_WIDTH,   // down-left, down
		WINDOW_WIDTH + 1, 1,              // down-right, right
		-WINDOW_WIDTH + 1, -WINDOW_WIDTH, // up-right, up
		-WINDOW_WIDTH - 1, -1 };          // up-left, left
	const int nextDir[8] = { 3, 0, 0, 1, 1, 2, 2, 3 }; // 0:down, 1:right, 2:up, 3:left

	int curPos = startPos;
	int curDir = 0;

	// special case : contour may revisit start position but not yet finished
	bool mayRevisit = false;
	if (inImage[curPos + WINDOW_WIDTH - 1] < 255 && inImage[curPos + WINDOW_WIDTH] == 255 &&
		(inImage[curPos + WINDOW_WIDTH + 1] < 255 || inImage[curPos + 1] < 255)) mayRevisit = true;

	outContour.clear();
	vector<unsigned char> outImage(WINDOW_HEIGHT*WINDOW_WIDTH, 255); // for debug
	while (true) {

		outContour.push_back(curPos);
		outImage[curPos] = 0;

		int guess;
		for (guess = 0; guess < 8; guess++) {
			int idx = (curDir * 2 + guess) % 8;
			int newPos = curPos + nextPos[idx];
			if (newPos < 0 || newPos >= WINDOW_HEIGHT*WINDOW_WIDTH) continue;
			if (inImage[newPos] < 255) {
				curPos = newPos;
				curDir = nextDir[idx];
				break;
			}
		}
		if (guess >= 8) return true; // isolated pixel
		if (curPos == startPos && (!mayRevisit || curDir != 1)) break;
	}

	//if (!PPMWriter::output("Style/3.transform/debug/debug-in.ppm", inImage, WINDOW_WIDTH, WINDOW_HEIGHT)) return false;
	//if (!PPMWriter::output("Style/3.transform/debug/debug-out.ppm", outImage, WINDOW_WIDTH, WINDOW_HEIGHT)) return false;
	//system("pause");

	return true;
}

bool FourierDescriptor::getDescriptor(
	vector<int> &inContour,
	pair<double, double> &inCenter,
	vector<double> &outDescriptor)
{
	if (inContour.empty()) { // empty image
		outDescriptor.assign(DESCRIPTOR_SIZE, 0);
		return true;
	}

	// get raw data

	int dim = (int)inContour.size();
	vector<double> cenDist(dim);
	for (int d = 0; d < dim; d++) {
		double dx = inContour[d] % WINDOW_WIDTH - inCenter.first;
		double dy = inContour[d] / WINDOW_WIDTH - inCenter.second;
		 cenDist[d] = sqrt(dx*dx + dy*dy);
	}

	// Discrete Fourier Transform

	double *pData = fftw_alloc_real(dim);
	fftw_complex *pResult = fftw_alloc_complex(dim/2+1);
	fftw_plan plan = fftw_plan_dft_r2c_1d(dim, pData, pResult, FFTW_ESTIMATE);

	for (int j = 0; j < dim; j++) pData[j] = cenDist[j];
	fftw_execute(plan);
	vector<double> spectrum(dim / 2 + 1);
	for (int j = 0; j < dim / 2 + 1; j++) {
		spectrum[j] = sqrt(pResult[j][0] * pResult[j][0] + pResult[j][1] * pResult[j][1]);
	}
	outDescriptor.resize(DESCRIPTOR_SIZE);
	for (int j = 1; j <= DESCRIPTOR_SIZE; j++) {
		if (j < dim / 2 + 1) {
			outDescriptor[j - 1] = spectrum[j] / spectrum[0];
		} else {
			outDescriptor[j - 1] = 0;
		}
	}

	fftw_destroy_plan(plan);
	fftw_free(pData);
	fftw_free(pResult);

	return true;
}

bool FourierDescriptor::isConnected(vector<unsigned char> &image) {

	// establish start point (consistent with "extractContour")
	int startPos = 0;
	while (startPos < WINDOW_HEIGHT*WINDOW_WIDTH) {
		if (image[startPos] < 255) break;
		startPos++;
	}
	if (startPos == WINDOW_HEIGHT*WINDOW_WIDTH) return true; // empty image

	// find connected component
	vector<bool> visited(image.size(), false);
	visited[startPos] = true;
	vector<int> queue(1, startPos);
	int head = 0;
	while (head < (int)queue.size()) {
		int curPos = queue[head];
		for (int y = -1; y < 2; y++) {
			for (int x = -1; x < 2; x++) {
				int newPos = curPos + y*WINDOW_WIDTH + x;
				if (newPos < 0 || newPos >= (int)image.size()) continue;
				if (!visited[newPos] && image[newPos] < 255) {
					visited[newPos] = true;
					queue.push_back(newPos);
				}
			}
		}
		head++;
	}

	// count points
	int componentSize = (int)queue.size();
	int totalSize = 0;
	for (int j = 0; j < WINDOW_HEIGHT*WINDOW_WIDTH; j++) {
		if (image[j] < 255) totalSize++;
	}

	return (double)componentSize > totalSize * 0.95;
}

bool FourierDescriptor::connect(vector<unsigned char> &image) {

	// connect multiple parts

	if (!isConnected(image)) {

		vector<unsigned char> bufferImage(image.begin(), image.end());

		vector<int> mask3X3e;
		vector<int> mask3X3d;
		for (int y = -1; y <= 1; y++) {
			for (int x = -1; x <= 1; x++) {
				mask3X3e.push_back(y*WINDOW_WIDTH + x);
				if (abs(x) + abs(y) < 2) mask3X3d.push_back(y*WINDOW_WIDTH + x);
			}
		}
		if (!erosion(bufferImage, mask3X3e)) return false;

		if (!isConnected(bufferImage)) {
			vector<int> mask5X5e;
			vector<int> mask5X5d;
			for (int y = -2; y <= 2; y++) {
				for (int x = -2; x <= 2; x++) {
					if (abs(x) + abs(y) < 4) mask5X5e.push_back(y*WINDOW_WIDTH + x);
					if (abs(x) + abs(y) < 3) mask5X5d.push_back(y*WINDOW_WIDTH + x);
				}
			}
			if (!erosion(bufferImage, mask5X5e)) return false;

			if (!isConnected(bufferImage)) {
				// fill the bounding box (?!)
				int minX = WINDOW_WIDTH;
				int minY = WINDOW_HEIGHT;
				int maxX = 0;
				int maxY = 0;
				for (int pos = 0; pos < WINDOW_HEIGHT*WINDOW_WIDTH; pos++) {
					if (bufferImage[pos] < 255) {
						minX = min(minX, pos%WINDOW_WIDTH);
						minY = min(minY, pos / WINDOW_WIDTH);
						maxX = max(maxX, pos%WINDOW_WIDTH);
						maxY = max(maxY, pos / WINDOW_WIDTH);
					}
				}
				for (int y = minY; y < maxY; y++) {
					for (int x = minX; x < maxX; x++) {
						bufferImage[y*WINDOW_WIDTH + x] = 0;
					}
				}

			}

			if (!dilation(bufferImage, mask5X5d)) return false;
		}

		if (!dilation(bufferImage, mask3X3d)) return false;

		image.swap(bufferImage);
	}

	return true;
}

bool FourierDescriptor::erosion(vector<unsigned char> &image, vector<int> &mask) {

	vector<unsigned char> buffer(image.begin(), image.end());
	image.assign(image.size(), 255);
	
#pragma omp parallel for shared(image, mask, buffer)
	for (int pos = 0; pos < WINDOW_HEIGHT*WINDOW_WIDTH; pos++) {
		if (buffer[pos] < 255) {
			for (int offset : mask) {
				int newPos = pos + offset;
				if(newPos>=0 && newPos<(int)image.size()) image[newPos] = 0;
			}
		}
	}

	return true;
}

bool FourierDescriptor::dilation(vector<unsigned char> &image, vector<int> &mask) {

	vector<unsigned char> buffer(image.begin(), image.end());
	image.assign(image.size(), 0);

#pragma omp parallel for shared(image, mask, buffer)
	for (int pos = 0; pos < WINDOW_HEIGHT*WINDOW_WIDTH; pos++) {
		if (buffer[pos] == 255) {
			for (int offset : mask) {
				int newPos = pos + offset;
				if (newPos>=0 && newPos<(int)image.size()) image[newPos] = 255;
			}
		}
	}

	return true;
}