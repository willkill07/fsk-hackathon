#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>

#include <openacc.h>

#include "params.h"

int main (int argc, char* argv[]) {

  const float DELTA = 0.2f;
	std::cerr << "Loading Data" << std::endl;
  loadData (argv[1]);
  sim.resize (sizes.size() * sizes.size());
	std::cerr << "Compute Kernel " << std::endl;
  computeSimilarity(subtrees.data(), subtrees.size(), offsets.data(), sizes.data(), offsets.size(), sim.data(), sim.size(), DELTA);
	std::cout << sim[0] << std::endl;
	#ifdef OUTPUT
  std::ofstream ofs (argv[2]);
  std::copy (sim.begin(), sim.end(), std::ostream_iterator<float> (ofs, "\n"));
  ofs.close();
	#endif
  return 0;
}

void loadData (const std::string & fileName) {
  int cumulativeIndex = 0;
  std::string fileToRead;
  std::ifstream ifs (fileName.c_str());
  while (std::getline (ifs, fileToRead)) {
    int currentIndex = 0;
    std::string line;
		std::ifstream bin (fileToRead.c_str());
    while (std::getline (bin, line)) {
      ++currentIndex;
      Subtree st;
      std::istringstream iss (line);
      std::copy (std::istream_iterator <float> (iss),
                 std::istream_iterator <float> (),
                 st.fv);
      subtrees.push_back (st);
    }
    offsets.push_back (cumulativeIndex);
    cumulativeIndex += currentIndex;
    sizes.push_back (currentIndex);
		if (sizes.size() >= LIMIT)
			break;
  }
}

#pragma acc routine seq
float simFunc (const Subtree * restrict s1, const Subtree * restrict s2) {
  float normdiff = 0.0f;
  for (int k = 0; k < FV_SIZE; ++k) {
    normdiff += fabsf (s1->fv[k] - s2->fv[k]) / fmaxf (1.0f, fmaxf (s1->fv[k], s2->fv[k]));
  }
  return normdiff;
}

void computeSimilarity (const Subtree * restrict data, const int subtreeCount,
                        const int * restrict offsets, const int * restrict sizes, const int binaryCount,
                        float * restrict sim, const int simCount, const float delta) {

  #pragma acc parallel loop collapse(2) copyin (data[0:subtreeCount], offsets[0:binaryCount], sizes[0:binaryCount]), copyout (sim[0:simCount]) gang
	for (int i1 = 0; i1 < binaryCount; ++i1) {
		for (int i2 = 0; i2 < binaryCount; ++i2) {
			float globalSim = 0.0f;
      #pragma acc loop collapse(2) vector
			for (int j1 = 0; j1 < sizes[i1]; ++j1) {
				for (int j2 = 0; j2 < sizes[i2]; ++j2) {
					const Subtree* s1 = data + offsets [i1] + j1;
					const Subtree* s2 = data + offsets [i2] + j2;
					float localSim = simFunc (s1, s2);
					if (localSim > delta) {
						globalSim += 1.0f;
					}
				}
			}
			int i = (binaryCount * (binaryCount - 1) / 2) - (binaryCount - i1) * (binaryCount - i1 - 1) / 2 + i2 - i1 - 1;
			sim [i] = globalSim / fminf (sizes[i1], sizes[i2]);
		}
	}
	std::cerr << "> Fin" << std::endl;
}
