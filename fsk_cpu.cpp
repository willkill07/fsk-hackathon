#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdio>

#define restrict __restrict__

#include "params.h"
#include "timer.h"

int main (int argc, char* argv[]) {
  const float DELTA = 0.35f;
  int params[4];
  timer start, end;
  std::transform (argv + 3, argv + 7, params, ::atoi);
  std::cerr << "Loading Data";
  start = getTime();
  loadData (argv[1]);
  end = getTime();
  fprintf (stderr, "Time (us): %12.2f\n", usTime(start, end));
  sim.resize ((params[1] - params[0]) * (params[3] - params[2]));
  std::cerr << "Compute Kernel ";
  start = getTime();
  computeSimilarity (subtrees.data(), offsets.data(), sizes.data(), params[0], params[1], params[2], params[3], sim.data(), DELTA);
  end = getTime();
  fprintf (stderr, "Time (us): %12.2f\n", usTime(start, end));
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
  const Subtree* data2 = data + offsets[LOOP2_START];

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
}
