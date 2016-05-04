#define LIMIT 2000

#define FV_SIZE 54

struct Subtree {
  float fv[FV_SIZE];
};

std::vector <Subtree> subtrees;
std::vector <int> offsets;
std::vector <int> sizes;
std::vector <float> sim;

void loadData (const std::string &);
void computeSimilarity (const Subtree * const restrict, const int * const restrict, const int * const restrict,
			const int, const int, const int, const int, float * const restrict, const float);
float simFunc (const Subtree * const restrict, const Subtree * const restrict);
