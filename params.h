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
void computeSimilarity (const Subtree * restrict, const int,
                        const int * restrict, const int * restrict, const int,
												float * restrict, const int, const float);
float simFunc (const Subtree * restrict, const Subtree * restrict);
