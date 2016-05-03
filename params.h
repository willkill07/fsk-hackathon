#define LIMIT 2000

#define FV_SIZE 54

#define START_L1 0
#define STOP_L1 1000
#define START_L2 1000
#define STOP_L2 2000

struct Subtree {
  float fv[FV_SIZE];
};

std::vector <Subtree> subtrees;
std::vector <int> offsets;
std::vector <int> sizes;
std::vector <float> sim;

void loadData (const std::string &);
void computeSimilarity (const Subtree * restrict data, const int * restrict offsets, const int * restrict sizes,
												const int LOOP1_START, const int LOOP1_END, const int LOOP2_START, const int LOOP2_END,
                        float * restrict sim, const float delta);
float simFunc (const Subtree * restrict, const Subtree * restrict);
