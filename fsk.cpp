#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdio>

#include <openacc.h>

#include "params.h"
#include "timer.h"

int main (int argc, char* argv[]) {
  const float DELTA = 0.2f;
	int params[4];
	std::transform (argv + 3, argv + 7, params, ::atoi);
  std::cerr << "Loading Data" << std::endl;
  loadData (argv[1]);
  sim.resize ((params[1] - params[0]) * (params[3] - params[2]));
  std::cerr << "Compute Kernel " << std::endl;
	timer start = getTime();
  computeSimilarity (subtrees.data(), offsets.data(), sizes.data(), params[0], params[1], params[2], params[3], sim.data(), DELTA);
	timer end = getTime();
	printf ("Time (us): %12.2f\n", usTime(start, end));
	std::copy (sim.begin(), sim.begin() + 5, std::ostream_iterator<float> (std::cout, " "));
	std::cout << std::endl;
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
  }
  offsets.push_back (cumulativeIndex);
  sizes.push_back (0);
}

#pragma acc routine
float simFunc (const Subtree * const restrict s1, const Subtree * const restrict s2) {
  float normdiff = 0.0f;
  for (int k = 0; k < FV_SIZE; ++k) {
    normdiff += fabsf (s1->fv[k] - s2->fv[k]) / fmaxf (1.0f, fminf (s1->fv[k], s2->fv[k]));
  }
  return normdiff / FV_SIZE;
}

void computeSimilarity (const Subtree * const restrict data, const int * const restrict offsets, const int * const restrict sizes,
                        const int LOOP1_START, const int LOOP1_END, const int LOOP2_START, const int LOOP2_END,
                        float * restrict sim, const float delta) {

  const Subtree* const data1 = data + offsets[LOOP1_START];
  int binSize1 = offsets[LOOP1_END] - offsets[LOOP1_START];

  const Subtree* const data2 = data + offsets[LOOP2_START];
  int binSize2 = offsets[LOOP2_END] - offsets[LOOP2_START];

  const int* const offsets1 = offsets + LOOP1_START;
  const int* const sizes1 = sizes + LOOP1_START;
  int size1 = LOOP1_END - LOOP1_START;

  const int* const offsets2 = offsets + LOOP2_START;
  const int* const sizes2 = sizes + LOOP2_START;
  int size2 = LOOP2_END - LOOP2_START;

  #pragma acc parallel loop gang collapse(2), \
    copyin (data1[0:binSize1], data2[0:binSize2], offsets1[0:size1], \
            offsets2[0:size2], sizes1[0:size1], sizes2[0:size2]), \
    copyout (sim[0:size1*size2])
  for (int i1 = 0; i1 < size1; ++i1) {
    for (int i2 = 0; i2 < size2; ++i2) {
      float globalSim = 0.0f;
      #pragma acc loop vector
      for (int j1 = 0; j1 < sizes1 [i1]; ++j1) {
        float localSim12 = 1.0f;
        for (int j2 = 0; j2 < sizes2 [i2]; ++j2) {
          localSim12 = fminf (localSim12, simFunc (data1 + offsets1 [i1] + j1, data2 + offsets1 [i2] + j2));
        }
        globalSim += (localSim12 < delta);
      }
      #pragma acc loop vector
      for (int j2 = 0; j2 < sizes2 [i2]; ++j2) {
        float localSim21 = 1.0f;
        for (int j1 = 0; j1 < sizes1 [i1]; ++j1) {
          localSim21 = fminf (localSim21, simFunc (data1 + offsets1 [i1] + j1, data2 + offsets1 [i2] + j2));
        }
        globalSim += (localSim21 < delta);
      }
      sim [size1 * i2 + i1] = 1.0f - globalSim / (sizes1[i1] + sizes2[i2]);
    }
  }
}
