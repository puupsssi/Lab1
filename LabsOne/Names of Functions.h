#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip> // Для setw
using namespace std;
/*Для тестовой и первой основной задачи :*/

/*double u_test - это функция, которая возвращает решение тестовой задачи, где
double x - текущая точка, в которой мы хотим посчитать точное решение тестовой задачи,
double u0 - начальное условие 
*/
double u_test(double x, double u0);

/*double test_f - это функция, которая возвращает значение функции f(x, u) в правой части дифференциального уравнения
из тестовой задачи, где 
double x - текущая точка, в которой мы хотим посчитать решение,
double v - текущая точка численной траектории, на основе которой мы хотим посчитать численное решение 
*/
double test_f(double x, double v);

/*double f_first_task - это функция, которая возвращает значение функции f(x, u) в правой части дифференциального уравнения
из первой основной задачи, где
double x - текущая точка, в которой мы хотим посчитать решение,
double v - текущая точка численной траектории, на основе которой мы хотим посчитать численное решение 
*/
double f_first_task(double x, double v);

/* pair<double, double> step_of_the_method_for_equation - функция, совершающая один шаг методом рунге кутта 4го порядка, 
возвращает следующую точку численной траектории (X_n+1, V_n+1)
double(*f)(double, double) - функция правой части дифференциального уравнения, на основе которой считается следующая точка численной траектории
double xn, double vn - текущая точка численной траектории
double h - шаг - то, на какую величину изменится x 
*/
pair<double, double> step_of_the_method_for_equation(double(*f)(double, double), double xn, double vn, double h);

/* double calculate_S - функция, которая вычисляет (и возвращает) величину S для контроля локальной погрешности
double vn_plus_1 - текущая точка, полученная методом с обычным шагом,
double v2n - текущая точка, полученная методом с половинным шагом
double epsilon - параметр контроля локальной погрешности
*/
double calculate_S(double vn_plus_1, double v2n, double epsilon);

/* int check_the_point - функция, контролирующая локальную погрешность, 
возвращает 0, если точка "плохая" (то есть уменьшаем шаг и возвращаемся в предыдущую точку),
возвращает 1, если точка хорошая и шаг подходящий
возвращает 2, если точка хорошая и погрешность такая маленькая, что можно увеличить шаг
double vn_plus_1 - текущая точка, полученная методом с обычным шагом,
double v2n - текущая точка, полученная методом с половинным шагом
double epsilon - параметр контроля локальной погрешности
double* h - текущий шаг (неявно возвращает новый шаг)
*/
int check_the_point(double vn_plus_1, double v2n, double* h, double epsilon);

/* vector<vector<double>> runge_kutta_4th_order - основная функция, соединяющая все воедино и считающая траекторию
возвращает явно:
vector<vector<double>> - сама численная и другие величины, требуемые для вывода, а именно:
если вычисляем траекторию с контролем локальной погрешности:
{ xn - текущая точка траектории,
  vn - текущее значение численной траектории, 
  v2n - текущее значение численной траектории, вычисленное с шагом в два раза меньше,
  (vn - v2n) - разница значиний с обычным и половинным шагом,
  S - величина для контроля локальной погрешности,
  hn - текущий шаг }

если вычисляем траекторию без контроля локальной погрешности:
{ xn- текущая точка траектории,
  vn - текущее значение численной траектории, 
  hn - текущий шаг }

  параметры:
  double(*f)(double, double) - функция правой части дифференциального уравнения, на основе которой считается следующая точка численной траектории
  double x0, double v0 - начальные условия
  double h - начальный шаг
  int n_steps - количество шагов
  double epsilon - параметр контроля локальной погрешности 
  int need_epsilon - нужно ли контролировать локальную погрешность 
  double right_boarder - левая граница отрезка
  double epsilon_boarder - в какой окрестности от границы нужно закончить счет
  vector<pair<double, double>>* changes_step - изменения шага (неявно возвращает пары изменения шага)
  vector<pair<double, double>>* test_task - неявно возвращает для тестовой задачи:
  пары (точное значение функции u(x) тестовой задачи, разница значений точного и численного решений в данной точке)
*/
vector<vector<double>> runge_kutta_4th_order(double(*f)(double, double), double x0, double v0,
    double h, int n_steps, double epsilon, int need_epsilon, double right_boarder, double epsilon_boarder,
    vector<pair<double, double>>* changes_step, vector<pair<double, double>>* test_task);


//Для второй основной задачи

/*double f_second_task - это функция, которая возвращает значение функции f(x, u) в правой части дифференциального уравнения
из первой основной задачи, где
double x - текущая точка, в которой мы хотим посчитать решение,
double u_1, double u_2 - текущие точки численной траектории, на основе которых мы хотим посчитать численное решение
double a, double b - параметры системы
*/
pair<double, double> f_second_task(double u_1, double u_2, double x, double a, double b);

/* tuple<double, double, double> step_of_the_method_for_the_system - функция, совершающая один шаг методом рунге кутта 4го порядка,
возвращает следующую точку численной траектории (X_n+1, V1_n+1, V2_n+1)
double(*f)(double, double) - функция правой части дифференциального уравнения, на основе которой считается следующая точка численной траектории
double xn, double v_1n, double v_2n - текущая точка численной траектории
double h - шаг - то, на какую величину изменится x
double a, double b - параметры системы
*/
tuple<double, double, double> step_of_the_method_for_the_system(pair<double, double>(*f)(double, double, double, double, double), double xn, double v_1n, double v_2n, double h, double a, double b);

/* double calculate_S_for_system - функция, которая вычисляет (и возвращает) величину S для контроля локальной погрешности
double v_1n_plus_1,double v_2n_plus_1 - текущая точка, полученная методом с обычным шагом,
double v_1n_plus_1,double v_2n_plus_1 - текущая точка, полученная методом с половинным шагом
double epsilon - параметр контроля локальной погрешности
*/
double calculate_S_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double epsilon);

/* int check_the_point_for_system - функция, контролирующая локальную погрешность,
возвращает 0, если точка "плохая" (то есть уменьшаем шаг и возвращаемся в предыдущую точку),
возвращает 1, если точка хорошая и шаг подходящий
возвращает 2, если точка хорошая и погрешность такая маленькая, что можно увеличить шаг
double v_1n_plus_1,double v_2n_plus_1 - текущая точка, полученная методом с обычным шагом,
double v_1n_plus_1,double v_2n_plus_1 - текущая точка, полученная методом с половинным шагом
double epsilon - параметр контроля локальной погрешности
double* h - текущий шаг (неявно возвращает новый шаг)
*/
int check_the_point_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double* h, double epsilon);


/* vector<vector<double>> runge_kutta_4th_order_for_system - основная функция, соединяющая все воедино и считающая траекторию
возвращает явно:
vector<vector<double>> - сама численная и другие величины, требуемые для вывода, а именно:
если вычисляем траекторию с контролем локальной погрешности:
{ xn - текущая точка траектории,
  (v_1n,v_2n) -текущее значение численной траектории,
  (v_12n,v_22n) - текущее значение численной траектории, вычисленное с шагом в два раза меньше,
  (v_1n - v_12n),(v_2n - v_22n)- разница значиний с обычным и половинным шагом,
  S - величина для контроля локальной погрешности,
  hn - текущий шаг }

если вычисляем траекторию без контроля локальной погрешности:
{ xn- текущая точка траектории,
  (v_1n, v_2n) - текущее значение численной траектории,
  hn - текущий шаг }

  параметры:
  pair<double, double>(*f)(double, double, double, double, double) - функция правой части дифференциального уравнения, на основе которой считается следующая точка численной траектории
  double x0, double v_10, double v_20 - начальные условия
  double a, double b - параметры системы
  double h - начальный шаг
  int n_steps - количество шагов
  double epsilon - параметр контроля локальной погрешности
  int need_epsilon - нужно ли контролировать локальную погрешность
  double right_boarder - левая граница отрезка
  double epsilon_boarder - в какой окрестности от границы нужно закончить счет
  vector<pair<double, double>>* changes_step - изменения шага (неявно возвращает пары изменения шага)
  vector<pair<double, double>>* test_task - неявно возвращает для тестовой задачи:
  пары (точное значение функции u(x) тестовой задачи, разница значений точного и численного решений в данной точке)
*/
vector<vector<double>> runge_kutta_4th_order_for_system(pair<double, double>(*f)(double, double, double, double, double), double x0, double v_10, double v_20, double a, double b,
    double h, int n_steps, double epsilon, int need_epsilon, double right_boarder, double epsilon_boarder,
    vector< pair<double, double>>* changes_step);
