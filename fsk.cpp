#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#define FV_SIZE 54

struct Subtree {
  float fv[FV_SIZE];
};

std::vector <Subtree> subtrees;
std::vector <int> offsets;
std::vector <int> sizes;
std::vector <float> sim;

void loadData (const std::string &);
void computeSimilarity (const Subtree * restrict, const int * restrict, const int * restrict, const int, const float, float * restrict);
void index (int n, int k, int * restrict, int * restrict);
float simFunc (const Subtree * restrict, const Subtree * restrict);

int main (int argc, char* argv[]) {

  const float DELTA = 0.7f;
  loadData (argv[1]);

  int subtreeCount = subtrees.size();
  Subtree * subtreeData = subtrees.data();

  int binCount = sizes.size();
  int * offsetData = offsets.data();
  int * sizeData = sizes.data();

  int simCount = binCount * (binCount - 1) / 2;
  sim.reserve (simCount);
  float * simData = sim.data();

  #pragma acc data copyin (subtreeData[0:subtreeCount], offsetData[0:binCount], sizeData[0:binCount]), copyout (simData[0:simCount])
  computeSimilarity(subtreeData, offsetData, sizeData, binCount, DELTA, simData);

  std::ofstream ofs (argv[2]);
  std::copy (sim.begin(), sim.end(), std::ostream_iterator<float> (ofs, "\n"));
  ofs.close();

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
}

#pragma acc routine seq
void index (int n, int k, int * restrict i, int * restrict j) {
  int _i = n - 2 - floor (sqrt (4 * n * (n - 1) - 8 * k - 7) / 2.0 - 0.5);
  *j = k + _i + 1 + ((n - _i) * ((n - _i) - 1) - n * (n - 1)) / 2;
  *i = _i;
}

#pragma acc routine seq
float simFunc (const Subtree * restrict s1, const Subtree * restrict s2) {
  float normdiff = 0.0f;
  for (int k = 0; k < FV_SIZE; ++k) {
    normdiff += fabsf (s1->fv[k] - s2->fv[k]) / fmaxf (s1->fv[k], s2->fv[k]);
  }
  return normdiff;
}

void computeSimilarity (const Subtree * restrict data,
                        const int * restrict offsets,
                        const int * restrict sizes,
                        const int binaryCount, const float delta,
                        float * restrict sim) {

  #pragma acc parallel loop collapse(2) gang
  for (int i1 = 0; i1 < binaryCount; ++i1) {
    for (int i2 = i1 + 1; i2 < binaryCount; ++i2) {
      float globalSim = 0.0f;
      #pragma acc loop collapse(2) vec
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
      sim [i] = globalSim / fminf (size1, size2);
    }
  }
}
