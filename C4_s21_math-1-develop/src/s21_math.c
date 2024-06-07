#include "s21_math.h"

int s21_abs(int x) {
  // Если x меньше 0, вернуть x умноженное на -1. Иначе, если x больше или
  // равно 0, вернуть x
  return x < 0 ? x * (-1) : x;
}

long double s21_fabs(double x) {
  // Если x меньше 0, вернуть x умноженное на -1. Иначе, если x больше или
  // равно 0, вернуть x
  return x < 0 ? x * (-1) : x;
}

long double s21_fmod(double x, double y) {
  long double remainder = 0.0;
  if (y == 0.0 || s21_is_nan(x) || s21_is_nan(y) || s21_is_inf(x)) {
    // Обработка случаев с NaN, inf и деления на ноль
    remainder = S21_NAN;
  } else if (s21_is_inf(y)) {
    remainder = x;
  } else {
    // Вычисляем целую часть от деления x на y
    double quotient = x / y;
    long long intPart = (long long)quotient;  // Приводим к целому числу
    // Находим остаток от деления
    remainder = x - (y * intPart);
  }

  return remainder;
}

long double s21_ceil(double x) {
  // Получаем целую часть числа с помощью приведения типов
  // Получаем целую часть числа
  long double result = (long double)((long long)x);
  // Проверяем, является ли значение x специальным или максимальным значением
  // double
  if (x == -S21_INF || x == S21_INF || s21_is_nan(x) || x == DBL_MAX)
    result = x;
  // Если число положительное и меньше исходного значения, увеличиваем его на 1
  if (x > 0 && x > result) {
    result += 1.0L;
    // Увеличиваем на 1, если число положительное и оно было меньше исходного
    // значения
  }
  // Возвращаем округленное значение
  return result;
}

long double s21_floor(double x) {
  // Получаем целую часть числа
  long double result = (long double)((long long)x);
  // Если x равен бесконечности, минус бесконечности или NaN, вернуть x
  if (x == -S21_INF || x == S21_INF || s21_is_nan(x)) result = x;
  // Если число отрицат. и оно было больше исходного значения, уменьшаем на 1
  if (x < 0 && x < result) {
    result -= 1.0L;
  }
  // Вернуть результат
  return result;
}

long double s21_sqrt(double x) {
  long double result = 0;
  // Если x меньше 0 или является NaN, вернуть NaN
  if (x < 0 || s21_is_nan(x)) {
    result = S21_NAN;
  }
  // Если x равен бесконечности или меньше или равен эпсилону корня, вернуть x
  else if (s21_is_inf(x) || x <= S21_SQRT_EPS) {
    result = x;
  }
  // Иначе вычислить корень из x с помощью функции pow
  else {
    result = s21_pow(x, 0.5);
  }
  // Вернуть результат
  return result;
}

long double s21_exp(double x) {
  long double result = 1.0;
  // Если x равен бесконечности, вернуть бесконечность
  if (x == S21_INF) {
    result = S21_INF;
  }
  // Если x является NaN, вернуть NaN
  else if (x != x) {
    result = S21_NAN;
  }
  // Если x равен минус бесконечности, вернуть 0
  else if (x == -S21_INF) {
    result = 0;
  }
  // Иначе вычислить экспоненту с помощью ряда Тейлора
  else {
    int n = 1;
    long double a = x, sn = 1;
    // Если x отрицательный, изменить знак и вычислять экспоненту положительного
    // числа
    if (x < 0.0) a = -x;
    // Вычислить экспоненту с помощью ряда Тейлора
    for (int i = 0; i < 1600; i++) {
      sn *= a / n++;
      result += sn;
    }
    // Если x отрицательный, вычислить обратное значение
    if (x < 0.0) result = 1 / result;
  }
  // Вернуть результат
  return result;
}

long double s21_log(double x) {
  long double result = 0;
  // Если x меньше 0, равен -бесконечности или NaN, результатом будет NaN
  if (x < 0 || x == -S21_INF || s21_is_nan(x)) {
    result = S21_NAN;
  }
  // Если x равен 0, результатом будет -бесконечность
  else if (x == 0) {
    result = -S21_INF;
  }
  // Если x равен бесконечности, результатом будет бесконечность
  else if (x == S21_INF) {
    result = S21_INF;
  }
  // Если x равен 1, результатом будет 0
  else if (x == 1) {
    result = 0;
  } else {
    // Инициализация переменных и вычисление целой части логарифма и остатка
    double N = 0.0, P = x, A = 0;
    while (P >= S21_E) {
      P /= S21_E;
      N++;  // Увеличиваем целую часть логарифма на 1
    }
    N += (P / S21_E);  // Добавляем остаток к целой части логарифма
    P = x;      // Восстанавливаем P
    int j = 0;  // Инициализируем счётчик итераций
    // Используем метод Ньютона для нахождения логарифма
    do {
      double L, R;
      A = N;  // Сохраняем предыдущее значение N
      L = (P / (s21_exp(N - 1.0)));  // Вычисляем левую границу
      R = ((N - 1.0) * S21_E);  // Вычисляем правую границу
      N = ((L + R) / S21_E);  // Вычисляем новое приближение логарифма
      j++;  // Увеличиваем счётчик итераций
    } while (N != A && j < 10000);
    // Цикл продолжается, пока значение N не перестанет
    // изменяться или пока число итераций не достигнет 10000
    // Результатом является значение N
    result = N;
  }
  return result;
}

long double s21_pow(double base, double exp) {
  long double result = 0;
  int tmp = 0;  // Инициализируем временную переменную
  if (base == 1 || exp == 0 || (base == -1 && (s21_is_inf(exp)))) {
    // Если base равен 1, или exp равен 0, или base равен -1 и exp является
    // бесконечностью
    result = 1;  // Результат равен 1
  } else if ((s21_fabs(base) < 1 && exp == -S21_INF) ||
             // Если base по модулю меньше 1 и exp равен минус бесконечности ИЛИ
             (s21_fabs(base) > 1 && exp == S21_INF) ||
             // Если base по модулю больше 1 и exp равен бесконечности ИЛИ
             (base == -S21_INF &&
              ((exp > 0 && (int)exp % 2 == 0) || exp == S21_INF)) ||
             // Если base равен минус бесконечности и exp чётное и больше 0 ИЛИ
             // exp равен бесконечности ИЛИ
             (base == S21_INF && exp > 0)) {
    // Если base равен бесконечности и exp больше 0
    result = S21_INF;  // Результат равен бесконечности
  } else if ((base < 0 && (int)exp - exp != 0) ||
             // Если base меньше 0 и exp не является целым числом ИЛИ
             (base != base && exp != exp)) {
    // Если base или exp равны NaN
    result = S21_NAN;  // Результат равен NaN
  } else if ((base == 0 && exp > 0) ||
             // Если base равен 0 и exp больше 0 ИЛИ
             (s21_fabs(base) > 1 && exp == -S21_INF) ||
             // Если base по модулю больше 1 и exp равен минус бесконечности ИЛИ
             (s21_fabs(base) < 1 && exp == S21_INF) ||
             // Если base по модулю меньше 1 и exp равен бесконечности ИЛИ
             (base == S21_NEGZERO && (exp > 0 && (int)exp % 2 == 0)) ||
             // Если base равен отрицательному нулю и exp чётное и больше 0 ИЛИ
             (base == -S21_INF && (exp < 0 && (int)exp % 2 == 0)) ||
             // Если base равен минус бесконечности и exp отрицательное и чётное
             (base == S21_INF && exp < 0)) {
    // Если base равен бесконечности и exp отрицательное
    result = 0;  // Результат равен 0
  } else if ((base == S21_NEGZERO && (exp > 0 && (int)exp % 2 != 0)) ||
             // Если base равен минус нулю и exp нечётное и больше 0 ИЛИ
             (base == -S21_INF && (exp < 0 && (int)exp % 2 != 0))) {
    // Если base равен минус бесконечности и exp отрицательное и нечётное
    result = S21_NEGZERO;  // Результат равен отрицательному нулю
  } else if (base == -S21_INF && (exp > 0 && (int)exp % 2 != 0)) {
    // Если base равен минус бесконечности и exp нечётное и больше 0
    result = -S21_INF;  // Результат равен минус бесконечности
  } else {
    if (base < 0) {  // Если base меньше 0
      if ((int)exp % 2 != 0) {  // Если exp не является целым числом
        tmp = 1;  // Устанавливаем значение временной переменной в 1
      }
      base = s21_fabs(base);  // Берем модуль base
    }
    result = tmp ? -s21_exp(exp * s21_log(base)) : s21_exp(exp * s21_log(base));
    // Если значение tmp равно 1, то  результат равен -exp * log(base), иначе
    // результат равен exp * log(base)
  }
  return result;  // Возвращаем результат
}

long double s21_acos(double x) {
  long double rezult = 0;
  if (x < -1.0 || x > 1.0) {
    rezult = S21_NAN;
  } else if (x == 1.0) {
    rezult = 0.0;
  } else if (x == -1.0) {
    rezult = S21_PI;
  } else {
    rezult = S21_PI_2 - s21_asin(x);
  }
  return rezult;
}

long double s21_asin(double x) {
  long double rezult = 0;
  if (x < -1.0 || x > 1.0) {
    rezult = S21_NAN;
  } else if (x == 1.0) {
    rezult = S21_PI_2;
  } else if (x == -1.0) {
    rezult = -S21_PI_2;
  } else {
    rezult = s21_atan(x / s21_sqrt(1.0 - x * x));
  }
  return rezult;
}

long double s21_atan(double x) {
  long double result = 0.0;
  if (x != 0.0) {
    // Если значение x превышает 1.0 или меньше -1.0, то вычисляем арктангенс
    // для x = 1/x
    if (x > 1.0) {
      result = S21_PI_2 - s21_atan(1.0 / x);
    } else if (x < -1.0) {
      result = -S21_PI_2 - s21_atan(1.0 / x);
    } else {
      long double power = x;
      long double term = x;
      int sign = 1;
      int n = 1;
      while (n < 1000000) {
        result += sign * term;
        power *= x * x;
        term = power / (2 * n + 1);
        sign = -sign;
        n++;
      }
    }
  }
  return result;
}

long double s21_sin(double x) {
  // Приводим угол к диапазону от 0 до 2π радиан
  double two_pi = 2.0 * S21_PI;
  x = s21_fmod(x, two_pi);
  // Если угол отрицательный, приводим его к положительному эквиваленту
  if (x < 0) x += two_pi;
  long double result = x;  // Инициализируем результат значением для n=1
  double term = x;  // Текущий член ряда
  int n = 2;  // Начинаем с n=2, так как n=0 и n=1 уже учли в результате
  while (s21_fabs(term) >= 1e-6) {
    term = -term * x * x / (n * (n + 1));
    result += term;
    n += 2;
  }
  return result;
}

long double s21_cos(double x) {
  // Приводим угол к диапазону от 0 до 2π радиан
  double two_pi = 2.0 * S21_PI;
  x = s21_fmod(x, two_pi);
  // Если угол отрицательный, приводим его к положительному эквиваленту
  if (x < 0) x += two_pi;
  long double result = 1.0;
  double term = 1.0;
  int n = 1;
  while (s21_fabs(term) >= 1e-6) {
    term = -term * x * x / (2 * n * (2 * n - 1));
    result += term;
    n++;
  }
  return result;
}

long double s21_tan(double x) {
  long double result = 0.0;
  result = s21_sin(x) / s21_cos(x);
  return result;
}
