#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <iomanip> // ��� setw
using namespace std;
/*��� �������� � ������ �������� ������ :*/

/*double u_test - ��� �������, ������� ���������� ������� �������� ������, ���
double x - ������� �����, � ������� �� ����� ��������� ������ ������� �������� ������,
double u0 - ��������� ������� 
*/
double u_test(double x, double u0);

/*double test_f - ��� �������, ������� ���������� �������� ������� f(x, u) � ������ ����� ����������������� ���������
�� �������� ������, ��� 
double x - ������� �����, � ������� �� ����� ��������� �������,
double v - ������� ����� ��������� ����������, �� ������ ������� �� ����� ��������� ��������� ������� 
*/
double test_f(double x, double v);

/*double f_first_task - ��� �������, ������� ���������� �������� ������� f(x, u) � ������ ����� ����������������� ���������
�� ������ �������� ������, ���
double x - ������� �����, � ������� �� ����� ��������� �������,
double v - ������� ����� ��������� ����������, �� ������ ������� �� ����� ��������� ��������� ������� 
*/
double f_first_task(double x, double v);

/* pair<double, double> step_of_the_method_for_equation - �������, ����������� ���� ��� ������� ����� ����� 4�� �������, 
���������� ��������� ����� ��������� ���������� (X_n+1, V_n+1)
double(*f)(double, double) - ������� ������ ����� ����������������� ���������, �� ������ ������� ��������� ��������� ����� ��������� ����������
double xn, double vn - ������� ����� ��������� ����������
double h - ��� - ��, �� ����� �������� ��������� x 
*/
pair<double, double> step_of_the_method_for_equation(double(*f)(double, double), double xn, double vn, double h);

/* double calculate_S - �������, ������� ��������� (� ����������) �������� S ��� �������� ��������� �����������
double vn_plus_1 - ������� �����, ���������� ������� � ������� �����,
double v2n - ������� �����, ���������� ������� � ���������� �����
double epsilon - �������� �������� ��������� �����������
*/
double calculate_S(double vn_plus_1, double v2n, double epsilon);

/* int check_the_point - �������, �������������� ��������� �����������, 
���������� 0, ���� ����� "������" (�� ���� ��������� ��� � ������������ � ���������� �����),
���������� 1, ���� ����� ������� � ��� ����������
���������� 2, ���� ����� ������� � ����������� ����� ���������, ��� ����� ��������� ���
double vn_plus_1 - ������� �����, ���������� ������� � ������� �����,
double v2n - ������� �����, ���������� ������� � ���������� �����
double epsilon - �������� �������� ��������� �����������
double* h - ������� ��� (������ ���������� ����� ���)
*/
int check_the_point(double vn_plus_1, double v2n, double* h, double epsilon);

/* vector<vector<double>> runge_kutta_4th_order - �������� �������, ����������� ��� ������� � ��������� ����������
���������� ����:
vector<vector<double>> - ���� ��������� � ������ ��������, ��������� ��� ������, � ������:
���� ��������� ���������� � ��������� ��������� �����������:
{ xn - ������� ����� ����������,
  vn - ������� �������� ��������� ����������, 
  v2n - ������� �������� ��������� ����������, ����������� � ����� � ��� ���� ������,
  (vn - v2n) - ������� �������� � ������� � ���������� �����,
  S - �������� ��� �������� ��������� �����������,
  hn - ������� ��� }

���� ��������� ���������� ��� �������� ��������� �����������:
{ xn- ������� ����� ����������,
  vn - ������� �������� ��������� ����������, 
  hn - ������� ��� }

  ���������:
  double(*f)(double, double) - ������� ������ ����� ����������������� ���������, �� ������ ������� ��������� ��������� ����� ��������� ����������
  double x0, double v0 - ��������� �������
  double h - ��������� ���
  int n_steps - ���������� �����
  double epsilon - �������� �������� ��������� ����������� 
  int need_epsilon - ����� �� �������������� ��������� ����������� 
  double right_boarder - ����� ������� �������
  double epsilon_boarder - � ����� ����������� �� ������� ����� ��������� ����
  vector<pair<double, double>>* changes_step - ��������� ���� (������ ���������� ���� ��������� ����)
  vector<pair<double, double>>* test_task - ������ ���������� ��� �������� ������:
  ���� (������ �������� ������� u(x) �������� ������, ������� �������� ������� � ���������� ������� � ������ �����)
*/
vector<vector<double>> runge_kutta_4th_order(double(*f)(double, double), double x0, double v0,
    double h, int n_steps, double epsilon, int need_epsilon, double right_boarder, double epsilon_boarder,
    vector<pair<double, double>>* changes_step, vector<pair<double, double>>* test_task);


//��� ������ �������� ������

/*double f_second_task - ��� �������, ������� ���������� �������� ������� f(x, u) � ������ ����� ����������������� ���������
�� ������ �������� ������, ���
double x - ������� �����, � ������� �� ����� ��������� �������,
double u_1, double u_2 - ������� ����� ��������� ����������, �� ������ ������� �� ����� ��������� ��������� �������
double a, double b - ��������� �������
*/
pair<double, double> f_second_task(double u_1, double u_2, double x, double a, double b);

/* tuple<double, double, double> step_of_the_method_for_the_system - �������, ����������� ���� ��� ������� ����� ����� 4�� �������,
���������� ��������� ����� ��������� ���������� (X_n+1, V1_n+1, V2_n+1)
double(*f)(double, double) - ������� ������ ����� ����������������� ���������, �� ������ ������� ��������� ��������� ����� ��������� ����������
double xn, double v_1n, double v_2n - ������� ����� ��������� ����������
double h - ��� - ��, �� ����� �������� ��������� x
double a, double b - ��������� �������
*/
tuple<double, double, double> step_of_the_method_for_the_system(pair<double, double>(*f)(double, double, double, double, double), double xn, double v_1n, double v_2n, double h, double a, double b);

/* double calculate_S_for_system - �������, ������� ��������� (� ����������) �������� S ��� �������� ��������� �����������
double v_1n_plus_1,double v_2n_plus_1 - ������� �����, ���������� ������� � ������� �����,
double v_1n_plus_1,double v_2n_plus_1 - ������� �����, ���������� ������� � ���������� �����
double epsilon - �������� �������� ��������� �����������
*/
double calculate_S_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double epsilon);

/* int check_the_point_for_system - �������, �������������� ��������� �����������,
���������� 0, ���� ����� "������" (�� ���� ��������� ��� � ������������ � ���������� �����),
���������� 1, ���� ����� ������� � ��� ����������
���������� 2, ���� ����� ������� � ����������� ����� ���������, ��� ����� ��������� ���
double v_1n_plus_1,double v_2n_plus_1 - ������� �����, ���������� ������� � ������� �����,
double v_1n_plus_1,double v_2n_plus_1 - ������� �����, ���������� ������� � ���������� �����
double epsilon - �������� �������� ��������� �����������
double* h - ������� ��� (������ ���������� ����� ���)
*/
int check_the_point_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double* h, double epsilon);


/* vector<vector<double>> runge_kutta_4th_order_for_system - �������� �������, ����������� ��� ������� � ��������� ����������
���������� ����:
vector<vector<double>> - ���� ��������� � ������ ��������, ��������� ��� ������, � ������:
���� ��������� ���������� � ��������� ��������� �����������:
{ xn - ������� ����� ����������,
  (v_1n,v_2n) -������� �������� ��������� ����������,
  (v_12n,v_22n) - ������� �������� ��������� ����������, ����������� � ����� � ��� ���� ������,
  (v_1n - v_12n),(v_2n - v_22n)- ������� �������� � ������� � ���������� �����,
  S - �������� ��� �������� ��������� �����������,
  hn - ������� ��� }

���� ��������� ���������� ��� �������� ��������� �����������:
{ xn- ������� ����� ����������,
  (v_1n, v_2n) - ������� �������� ��������� ����������,
  hn - ������� ��� }

  ���������:
  pair<double, double>(*f)(double, double, double, double, double) - ������� ������ ����� ����������������� ���������, �� ������ ������� ��������� ��������� ����� ��������� ����������
  double x0, double v_10, double v_20 - ��������� �������
  double a, double b - ��������� �������
  double h - ��������� ���
  int n_steps - ���������� �����
  double epsilon - �������� �������� ��������� �����������
  int need_epsilon - ����� �� �������������� ��������� �����������
  double right_boarder - ����� ������� �������
  double epsilon_boarder - � ����� ����������� �� ������� ����� ��������� ����
  vector<pair<double, double>>* changes_step - ��������� ���� (������ ���������� ���� ��������� ����)
  vector<pair<double, double>>* test_task - ������ ���������� ��� �������� ������:
  ���� (������ �������� ������� u(x) �������� ������, ������� �������� ������� � ���������� ������� � ������ �����)
*/
vector<vector<double>> runge_kutta_4th_order_for_system(pair<double, double>(*f)(double, double, double, double, double), double x0, double v_10, double v_20, double a, double b,
    double h, int n_steps, double epsilon, int need_epsilon, double right_boarder, double epsilon_boarder,
    vector< pair<double, double>>* changes_step);
