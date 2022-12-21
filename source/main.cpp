#include "header.hpp"

int main(int argc, char* argv[]) {
  // 1: размер, 2: ограничитель печати матрицы, 3: точность вычислений
  // 4: левый конец интервала, 5: правый конец интервала
  // 6: номер формулы заполнения матрицы (см double HelperInMat(...))
  // 7:(если 6: = 0) имя файла c матрицой
  int n, m, k, err;
  double eps, leftb, rightb;
  err = TestInitArg(argc, argv, &n, &m, &k, &eps, &leftb, &rightb);
  switch (err) {
    case -1:
      printf("Недостаточное количество аргументов\n");
      return -1;

    case -2:
      printf("Некорректные аргументы\n");
      return -2;

    default:
      break;
  }

  if (n <= 0) {
    printf("Некорректный размер матрицы\n");
    return -3;
  }

  double* matrx = new double[n*n];
  if (k == 0)
    err = InMat(n, k, matrx, argv[7]);
  else
    err = InMat(n, k, matrx, nullptr);
  switch (err) {
    case -1:
      std::cout << "Ошибка открытия файла" << std::endl;
      delete[] matrx;
      return -4;

    case -2:
      std::cout << "Неверный формат данных" << std::endl;
      delete[] matrx;
      return -5;

    case -3:
      std::cout << "Ошибка ввода-вывода при чтении" << std::endl;
      delete[] matrx;
      return -6;

    case -4:
      std::cout << "Недостаточное количество элементов" << std::endl;
      delete[] matrx;
      return -7;

    default:
      break;
  }

  //PrintMat(matrx, n, n);
///////////// Конец проверок

  double* x = new double[n + 1];
  err = 0;
  double time = (double)clock();
  err = BisecAlg (n, matrx, x, eps, leftb, rightb);
  time = (double)(clock() - time) / CLOCKS_PER_SEC;

  switch (err) {
    case -1:
      printf("Матрица несимметрична\n");
      delete[] matrx;
      delete[] x;
      return -8;

    case -2:
      printf("Некорректно заданный отрезок\n");
      delete[] matrx;
      delete[] x;
      return -9;

    default:
      break;
  }

  printf("Количество собственных значений: %d\n", (int)x[0]);
  for (int i = 0, limiter = std::min(m, int(x[0])); i < limiter; i++)
    printf("%.20lf\n", x[i + 1]);

  printf("\nВремя алгоритма: %.3lf\n", time);
  if((int)x[0] == n) {
    printf("Невязка в 1-ом инварианте: %lf\n", Residual(matrx, n, x, 1));
    printf("Невязка во 2-ом инварианте: %lf\n", Residual(matrx, n, x, 2));
  }

  delete[] matrx;
  delete[] x;
  return 0;
}
