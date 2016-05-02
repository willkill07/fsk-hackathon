const int FV_SIZE = 54;

struct Subtree {
	float fv[FV_SIZE];
};

float computeSubtreeKernel (const Subtree* begin1, const Subtree* begin2, int size1, int size2) {
	float globalSim = 0.0f;
	for (int i1 = 0; i1 < size1; ++i1) {
		for (int i2 = 0; i2 < size2; ++i2) {
			Subtree c1 = begin1 [i1];
			Subtree c2 = begin2 [i2];
			float normdiff = 0.0f;
			for (int i = 0; i < FV_SIZE; ++i) {
				normdiff += fabs (c1.fv[i] - c2.fv[i]) / fmax (c1.fv[i], c2.fv[i]);
			}
			if (normdiff > delta) {
				globalSim += normdiff;
			}
		}
	}
	return globalSim / (size1 * size2);
}
