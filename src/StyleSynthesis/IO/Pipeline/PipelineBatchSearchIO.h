#pragma once

#include <string>

#include "IO/BaseIO.h"

using namespace std;

namespace StyleSynthesis {

	class ContextPartGraph;

	class PipelineBatchSearchIO : public BaseIO {

	public:

		static bool process();

	private:

		static bool runModelPairs(string sourceName, string targetName, string pairName);
		static bool organizeResults(string sourceName, string targetName, string pairName);

		static bool error(string s);
	};

}