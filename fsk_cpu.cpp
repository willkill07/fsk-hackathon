#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <chrono>

#include <cmath>

#define restrict __restrict__

#include "params.h"

int main (int argc, char* argv[]) {

  const float DELTA = 0.2f;
  std::cerr << "Loading Data" << std::endl;
  loadData (argv[1]);
  sim.resize ((STOP_L1 - START_L1) * (STOP_L2 - START_L2));
  std::cerr << "Compute Kernel " << std::endl;
  auto begin = std::chrono::high_resolution_clock::now();
  computeSimilarity (subtrees.data(), offsets.data(), sizes.data(), START_L1, STOP_L1, START_L2, STOP_L2, sim.data(), DELTA);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "TIME: " << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count() << " us" << std::endl;
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

float simFunc (const Subtree * restrict s1, const Subtree * restrict s2) {
  float normdiff = 0.0f;
  #pragma omp simd
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
  
  #pragma omp parallel for collapse(2) schedule(dynamic,1)
  for (int i1 = 0; i1 < size1; ++i1) {
    for (int i2 = 0; i2 < size2; ++i2) {
      float globalSim = 0.0f;
      const Subtree* base1 = data1 + offsets1 [i1];
      const Subtree* base2 = data2 + offsets2 [i2];
      for (int j1 = 0; j1 < sizes1 [i1]; ++j1) {
        float localSim = 1.0f;
        for (int j2 = 0; j2 < sizes2 [i2]; ++j2) {
          localSim = fminf (localSim, simFunc (base1 + j1, base2 + j2));
        }
        globalSim += (localSim < delta);
      }
      for (int j2 = 0; j2 < sizes2 [i2]; ++j2) {
        float localSim = 1.0f;
        for (int j1 = 0; j1 < sizes1 [i1]; ++j1) {
          localSim = fminf (localSim, simFunc (base1 + j1, base2 + j2));
        }
        globalSim += (localSim < delta);
      }
      sim [size1 * i2 + i1] = 1.0f - globalSim / (sizes1[i1] + sizes2[i2]);
    }
  }
  std::cerr << "> Fin" << std::endl;
}
