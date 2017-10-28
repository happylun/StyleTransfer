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

#include <vector>
#include <string>

#include "Library/CMLHelper.h"

#include "Data/StyleSynthesisTypes.h"

using namespace std;

namespace StyleSynthesis {

	class SampleSimplePoissonDisk {

	private:

		struct TPoissonSample {
			vec3 position;
			vec3 normal;
			int faceID;
		};

		static const double RELATIVE_RADIUS;

	public:

		SampleSimplePoissonDisk(TTriangleMesh *mesh);
		~SampleSimplePoissonDisk();

	public:

		bool runSampling(int numSamples);
		bool runSampling(double sampleRadius);
		bool exportSample(TSampleSet &samples);

	private:

		bool calculateCDF();
		bool buildGrids();
		bool generateSamples();

		void clearUp();

	private:

		bool sampleOnTriangle(TPoissonSample &sample, int faceID);
		bool checkGrid(TPoissonSample &sample, int &gridPos);

	private:

		TTriangleMesh *mpMesh;

	private:

		vector<double> mTriangleAreaCDF; // CDF normalized to [0..1] : # of faces
		double mTotalArea; // mesh total area
		int mDesiredNumSamples; // desired number of samples
		double mSampleRadius; // sample points' minimum distance
		double mGridRadius; // grid's side length (may be larger than mSampleRadius)
		vec3 mBBMin; // mesh bounding box corner
		vec3 mBBMax; // mesh bounding box corner
		vector<vector<int>> mGrids; // sampleID : # of samples in grid cell : (gridXSize*gridYSize*gridZSize)
		vec3i mGridSize; // dimensions of grids
		vector<TPoissonSample> mSamples; // sample : # of samples
	};
}