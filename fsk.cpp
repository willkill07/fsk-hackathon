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
  sim.resize ((STOP_L1 - START_L1) * (STOP_L2 - START_L2));
	std::cerr << "Compute Kernel " << std::endl;
  computeSimilarity (subtrees.data(), offsets.data(), sizes.data(), START_L1, STOP_L1, START_L2, STOP_L2, sim.data(), DELTA);
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
	offsets.push_back (cumulativeIndex);
	sizes.push_back (0);
}

#pragma acc routine seq
float simFunc (const Subtree * restrict s1, const Subtree * restrict s2) {
  float normdiff = 0.0f;
  for (int k = 0; k < FV_SIZE; ++k) {
    normdiff += fabsf (s1->fv[k] - s2->fv[k]) / fmaxf (1.0f, fmaxf (s1->fv[k], s2->fv[k]));
  }
  return normdiff;
}

void computeSimilarity (const Subtree * restrict data, const int * restrict offsets, const int * restrict sizes,
												const int LOOP1_START, const int LOOP1_END, const int LOOP2_START, const int LOOP2_END,
                        float * restrict sim, const float delta) {

	const Subtree* data1 = data + offsets[LOOP1_START];
	int binSize1 = offsets[LOOP1_END] - offsets[LOOP1_START];

	const Subtree* data2 = data + offsets[LOOP2_START];
	int binSize2 = offsets[LOOP2_END] - offsets[LOOP2_START];

	const int* offsets1 = offsets + LOOP1_START;
	const int* sizes1 = sizes + LOOP1_START;
	int size1 = LOOP1_END - LOOP1_START;

	const int* offsets2 = offsets + LOOP2_START;
	const int* sizes2 = sizes + LOOP2_START;
	int size2 = LOOP2_END - LOOP2_START;


  #if USE_DATA
		{
			int begin = std::min(offsets[LOOP1_START], offsets[LOOP2_START]);
			int end = std::max(offsets[LOOP1_END], offsets[LOOP2_END]);
			int size = end - begin;
      #pragma acc enter data copyin(data[begin:size])
		}
		{
			int begin = std::min(LOOP1_START, LOOP2_START);
			int end = std::max(LOOP1_END, LOOP2_END);
			int size = end - begin;
      #pragma acc enter data copyin(offsets[begin:size],sizes[begin:size])
		}
    #pragma acc parallel loop collapse(2) \
			present (data1,data2,offsets1,offsets2,sizes1,sizes2),	\
			copyout (sim[0:size1*size2])
	#else
  	#pragma acc parallel loop collapse(2) \
		  copyin (data1[0:binSize1],data2[0:binSize2],offsets1[0:size1],offsets2[0:size2],sizes1[0:size1],sizes2[0:size2]), \
		  copyout (sim[0:size1*size2])
	#endif
	for (int i1 = 0; i1 < size1; ++i1) {
		for (int i2 = 0; i2 < size2; ++i2) {
			float globalSim = 0.0f;
      #pragma acc loop collapse(2) vector
			for (int j1 = 0; j1 < sizes1 [i1]; ++j1) {
				for (int j2 = 0; j2 < sizes2 [i2]; ++j2) {
					const Subtree* s1 = data1 + offsets1 [i1] + j1;
					const Subtree* s2 = data2 + offsets2 [i2] + j2;
					float localSim = simFunc (s1, s2);
					if (localSim > delta) {
						globalSim += 1.0f;
					}
				}
			}
			sim [size1 * i2 + i1] = globalSim / fminf (sizes1[i1], sizes2[i2]);
		}
	}

	#if USE_DATA
    #pragma acc exit data delete(data,offsets,sizes)
	#endif
	std::cerr << "> Fin" << std::endl;
}
