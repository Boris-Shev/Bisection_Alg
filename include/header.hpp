#ifndef _HEADER_HPP_
#define _HEADER_HPP_

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <limits>
#include <type_traits>
#include <string>

int TestInitArg (int argc, char* argv[], int* n, int* m, int* k, double* eps,
                 double* leftb, double* rightb);

int InMat (int size, int formula, double* matrx, char* file);

double HelperInMat (int formula, int size, int i, int j);

void PrintMat (double* matrx, int numRow, int numCol, int limiter);

void PrintMat (double* matrx, int numRow, int numCol);

int BisecAlg (int n, double* matrx, double* x, double eps, double leftb,
              double rightb);

bool IsSymmetrical(double* matrx, int n);

int n_(int n, double* a, double lambda);

void Tridiagonal(double* matrx, int n);

double Residual (double* matrx, int n, double* x, int mode);

double Inaccuracy (double* x, int size);


#endif // _HEADER_HPP_
