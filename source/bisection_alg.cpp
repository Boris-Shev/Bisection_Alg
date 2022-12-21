#include "header.hpp"

//TODO: Разобраться как алгоритм обрабатывает невырожденные матрицы
//TODO: Убрать глобальную переменную
int Iter = 0;

int BisecAlg (int n, double* A, double* x, double eps, double leftb,
              double rightb) {
  if (!IsSymmetrical(A, n))
    return -1;

  if (rightb < leftb)
    return -2;

  //PrintMat(A, n, n);
  Tridiagonal(A, n);
	//PrintMat(A, n, n);


  // int k = 20;
  // for (int i = 0; i <= k; i++) {
  //   printf("%.2f:  %d\n", leftb + i * (rightb - leftb)/k, n_(n, A, leftb + i * (rightb - leftb)/k + eps));
  // }
  // return -1;
  //printf("%.2f:  %d\n", 1., n_(n, A, rightb));
  //return -1;

  rightb += eps;
	leftb -= eps;
	x[0] = (double)n_(n, A, rightb) - n_(n, A, leftb);
	if (std::fabs(x[0]) < std::numeric_limits<double>::epsilon())
		return 0;
////////////////////////////////////////////////////////////////////
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

    //printf("test1: %d\n", doneCount);
		curMid = 0.5 * (curLeft + curRight);
    //printf("curMid: %.20lf\n", curMid);
    count = n_(n, A, curRight) - n_(n, A, curLeft);
		for (int j = 0; j < count; j++) {
        //printf("%d  %d  %d  %d  %f  %f  %f\n", doneCount, j, count, doneCount + j + 1, curLeft, curMid, curRight);
        //printf("%d:  %d\n", j, doneCount + j + 1);
        x[doneCount + j + 1] = curMid;
    }

		doneCount += count;
		curLeft = curMid;
		curRight = rightb;
	}

  printf("Iter = %d\n", Iter);
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
  // REMIND: Проверять что х и у не равны одновременно 0
  // REMIND: а = alphaInv * (a - lambda);    b = alphaInv * b;
  // TODO: Рассматривать случаи когда внедиагональные элементы 0
  Iter ++;
  lambda += 1e-10;

  double alphaInv = std::fabs(A[0 * n + 0] - lambda);
  for (int i = 1; i < n; i++) {
    if (std::fabs(A[i * n + i]) > alphaInv)
      alphaInv = std::fabs(A[i * n + i] - lambda);

    if (std::fabs(A[i * n + i - 1]) > alphaInv)
      alphaInv = std::fabs(A[i * n + i - 1]);
  }
  if (std::fabs(alphaInv) < std::numeric_limits<double>::epsilon()) {
    if (lambda < 0)
      return 0;
    else
      return n;
  }
  alphaInv = 1. / (4 * alphaInv);

  int count = 0;
  double w = std::numeric_limits<double>::epsilon();
  double M = 1. / (8. * w);
  double x = alphaInv * (A[0] - lambda);
  double y = 1;
  if (A[0] < lambda)
    count++;

  double a, b;
  double u, v;
  double t, q, l, m;
  for (int i = 1; i < n; i++) {

    if (std::fabs(x) < w) {
      if (std::fabs(y) < w) {
        if (lambda > 0)
          count += n - i - 1;
        break;
      }
      x = - Sgn(y);
      y = 0;
      count++;
      continue;
    }

    a = alphaInv * (A[i * n + i] - lambda);
    b = alphaInv * A[i * n + i - 1];

    if (std::fabs(b) < w) {
      if (std::fabs(a) < w) {
        y = x;
        x = 0;
        continue;
      }
      if (a < 0)
        count++;
      y = x;
      x = - Sgn(x);
      continue;
    }

    t = std::fmax(std::fabs(b * (b * y)), std::fabs(x));
    v = (t < 1) * ((x / t) * M) + !(t < 1) * ((M / t) * x);
    q = a * v;

    m = (M * b) * b;
    l = (((m < 1) && (t < 1)) || (!(m < 1) && !(t < 1))) * ((m / t) * y) +
        (((m < 1) && !(t < 1)) || (!(m < 1) && (t < 1))) * ((y / t) * m);
    u = q - l;

    if (std::fabs(u) < w){
      x = 0;
      y = v; // maybe problems
      continue;
    }

    if ((x > 0 && u < 0) || (x < 0 && u > 0))
      count++;
    x = u; // maybe problems
    y = v; // maybe problems
  }

	return count;
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

int Sgn(double x) {
  return (x > 0) - (x < 0);
}
