#include "header.hpp"

int TestInitArg (int argc, char* argv[], int* n, int* m, int* k, double* eps) {
  if (argc < 5)
    return -1;

  for (int i = 1; i < 5; i++) {
    std::string s(argv[i]);
    std::string::const_iterator it = s.begin();
    if (i != 3) {
      while (it != s.end() && std::isdigit(*it)) ++it;
      if (!(!s.empty() && it == s.end()))
        return -2;
    } else {
  M:  while (it != s.end() && std::isdigit(*it)) ++it;
      if (*it == '.') {
        it++; goto M;
      }
      if (*it == 'e') {
        if (*(it + 1) == '-') {
          if (it + 2 == s.end() || it == s.begin())
            return -2;
          it += 2; goto M;
        } else {
          if (it + 1 == s.end() || it == s.begin())
            return -2;
          it += 1; goto M;
        }
      }
      if (!(!s.empty() && it == s.end()))
        return -2;
      }
  }
  *n = std::stoi(argv[1]);
  *m = std::stoi(argv[2]);
  *eps = std::stod(argv[3]);
  *k = std::stoi(argv[4]);
  return 0;
}

int InMat (int size, int formula, double* matrx, char* file) {
  if (formula == 0) {
    std::ifstream fin(file);
    if(!fin.is_open()){ // Ошибка открытия файла
      return -1;
    }
    int i;
    for (i = 0; fin >> matrx[i] && i < size*size; i++) {}
    if (fin.eof() && i == size*size)
      {}
    else if (fin.eof() && i != size*size) // Недостаточное количество элементов
      return -4;
    else if (fin.fail()) // Неверный формат данных
      return -2;
    else if (fin.bad()) // Ошибка ввода-вывода при чтении
      return -3;
  } else {
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        matrx[i*size + j] = HelperInMat(formula, size, i, j);
  }
  return 0;
}

double HelperInMat (int formula, int size, int i, int j) {
  switch (formula) {
    case 1:
      return size - std::max(i + 1, j + 1) + 1;
    case 2:
      return 2 * (i == j) - 1 * (std::abs(i - j) == 1);
    case 3:
      return 1 * (j == i && j + 1 < size) + (i + 1) * (j + 1 == size) +
              + (j + 1) * ((i + 1) == size) - size * (i == j && i == size - 1);
    case 4:
      return 1. / (i + j + 1);
    default:
      return 0;
  }
  return 0;
}

void PrintMat (double* matrx, int numRow, int numCol, int limiter) {
  if (limiter <= numRow) {
    for (int i = 0; i < limiter; i++) {
      for (int j = 0; j < numCol; j++) {
        printf(" %10.3e", matrx[i*numCol + j]);
      }
      printf("\n");
    }
  } else {
    for (int i = 0; i < numRow; i++) {
      for (int j = 0; j < numCol; j++) {
        printf(" %10.3e", matrx[i*numCol + j]);
      }
      printf("\n");
    }
  }
}

double Residual (double* A, int n, double* b, double* x) {
  double norm_b = 0;
  double norm_num = 0;
  double sum;

  for (int i = 0; i < n; i++) {
    norm_b += b[i] * b[i];
  }

  for (int i = 0; i < n; i++) {
    sum = 0;
    for (int k = 0; k < n; k++) {
      sum += A[i*n + k] * x[k];
    }
    norm_num += (sum - b[i]) * (sum - b[i]);
  }
  return std::sqrt(std::fabs(norm_num)) / std::sqrt(std::fabs(norm_b));
}

double Inaccuracy (double* x, int size) {
  double sum = 0;
  for (int i = 0; i < size; i += 2) {
    sum += (x[i] - 1) * (x[i] - 1);
  }
  for (int i = 1; i < size; i += 2) {
    sum += x[i] * x[i];
  }
  return std::sqrt(std::fabs(sum));
}
