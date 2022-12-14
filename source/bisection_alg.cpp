#include "header.hpp"

int BisecAlg (int n, double* A, double* x, double eps, double leftb,
              double rightb) {
  if (!IsSymmetrical(A, n))
    return -1;

  if (rightb < leftb)
    return -2;

  PrintMat(A, n, n);
  Tridiagonal(A, n);
	PrintMat(A, n, n);

  rightb += std::numeric_limits<double>::epsilon();
	leftb -= std::numeric_limits<double>::epsilon();
	x[0] = (double)n_(n, A, rightb) - n_(n, A, leftb);
	if (std::fabs(x[0]) < std::numeric_limits<double>::epsilon())
		return 0;

	int beforeLeftBorderCount = n_(n, A, leftb);
	double curLeft = leftb;
	double curRight = rightb;
  double curMid;

  int doneCount = 0, count;
	while (doneCount < (int)x[0]) {
		while (curRight - curLeft > eps) {
			curMid = 0.5 * (curLeft + curRight);

			if (n_(n, A, curMid) < doneCount + 1 + beforeLeftBorderCount)
				curLeft = curMid;
			else
				curRight = curMid;
		}

		curMid = 0.5 * (curLeft + curRight);
		count = n_(n, A, curRight) - n_(n, A, curLeft);
		for (int j = 0; j < count; j++) {
        //printf("%d  %d  %d  %d  %f  %f  %f\n", doneCount, j, count, doneCount + j + 1, curLeft, curMid, curRight);
        x[doneCount + j + 1] = curMid;
    }

		doneCount += count;
		curLeft = curMid;
		curRight = rightb;
	}

	return 0;
}

bool IsSymmetrical(double* A, int n) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      if (std::fabs(A[i * n + j] - A[j * n + i]) >
       std::numeric_limits<double>::epsilon())
        return false;
    }
  }
  return true;
}

int n_(int n, double* A, double lambda)
{
	int res;
	double x, y;

	x = A[0 * n + 0] - lambda;
	y = 1.0;
	if (x < 0.)
    res = 1;
  else
    res = 0;

  double u, v, tmp, a_k, b_k1, gamma;
	for (int i = 1; i < n; i++)
	{
		a_k = A[i * n + i] - lambda;
		b_k1 = A[i * n + i - 1];

		tmp = std::fabs(b_k1 * b_k1 * y);
		if (std::fabs(x) > tmp)
			tmp = std::fabs(x);

		if (tmp < std::numeric_limits<double>::epsilon())
			tmp = 10 * std::numeric_limits<double>::epsilon();

		gamma = (1 / std::numeric_limits<double>::epsilon()) / tmp;
		u = gamma * (a_k * x - b_k1 * b_k1 * y);
		v = gamma * x;
		if (u * x < 0.0)
			res++;
		x = u;
		y = v;
	}

	return res;
}

void Tridiagonal(double* A, int n)
{
	double x, y, r, s;
	double A_ii, A_ij, A_ji, A_jj;
  double cos, sin;

	for (int i = 1; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			x = A[i * n + i - 1];
			y = A[j * n + i - 1];
			if (std::fabs(y) < std::numeric_limits<double>::epsilon())
				continue;

			r = std::sqrt(x * x + y * y);
			if (r < std::numeric_limits<double>::epsilon())
				continue;

			cos = x / r;
			sin = -y / r;
			A[i * n + i - 1] = A[(i - 1) * n + i] = r;
			A[j * n + i - 1] = A[(i - 1) * n + j] = 0.0;

			for (int k = i + 1; k < n; k++)
			{
				if (k == j)
					continue;
				x = A[i * n + k];
				y = A[j * n + k];
				A[k * n + i] = A[i * n + k] = x * cos - y * sin;
				A[k * n + j] = A[j * n + k] = x * sin + y * cos;
			}
			x = A[i * n + i];
			y = A[j * n + j];
			r = A[i * n + j];
			s = A[j * n + i];

			A_ii = x * cos - s * sin;
			A_ji = x * sin + s * cos;
			A_ij = r * cos - y * sin;
			A_jj = r * sin + y * cos;

			A[i * n + i] = A_ii * cos - A_ij * sin;
			A[j * n + i] = A_ii * sin + A_ij * cos;
			A[i * n + j] = A[j * n + i];
			A[j * n + j] = A_ji * sin + A_jj * cos;
		}
	}
}
