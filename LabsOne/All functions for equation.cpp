#include"Names of Functions.h"
double u_test(double x, double u0) {
    return u0 * exp(x);
}

double test_f(double x, double v) {
    return v;
}

double f_first_task(double x, double v) {
    return (x / (pow(x, 2) + 1)) * pow(v, 2) + v - pow(v, 3) * sin(10 * x);
}

pair<double, double> step_of_the_method_for_equation(double(*f)(double, double), double xn, double vn, double h) {
    // ���������� ������������� k1, k2, k3, k4 - ��� ������ ����� ����� 4�� �������
    double k1 = f(xn, vn);
    double k2 = f(xn + h / 2.0, vn + h / 2.0 * k1);
    double k3 = f(xn + h / 2.0, vn + h / 2.0 * k2);
    double k4 = f(xn + h, vn + h * k3);

    // ���������� �������� ���������� x � v
    xn = xn + h;  // ��������� x �� ���� ���
    vn = vn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);  // ��������� v
    return { xn, vn };
}
double calculate_S(double vn_plus_1, double v2n, double epsilon) {
    //�������� ��� �������� ��������� �����������
    double S = abs(v2n - vn_plus_1) / 15.0; //2^4 - 1 = 15;
    return S;
}
int check_the_point(double vn_plus_1, double v2n, double* h, double epsilon) {
    double S = calculate_S(vn_plus_1, v2n, epsilon);
    if ((epsilon / 32.0) <= S && S <= epsilon) { //������� �����
        return 1; //������, ��� ����� ������� � �� ���������� ����
    }
    if (S < (epsilon / 32.0)) {
        *h = *h * 2;
        return 2; //���������� � ���������� �������, ������ ��� �������� ���
    }
    if (S > epsilon) {//����� ������, ������ ��� � �������� ��� �� ����� (xn,vn)
        *h = *h / 2.0;
        return 0;
    }
    return 0;
}
// ������� ��� ������ �����-����� 4-�� �������
vector<vector<double>> runge_kutta_4th_order(double(*f)(double, double), double x0, double v0,//������� �������, ������� �������� �� ������ �����
    double h, int n_steps, double epsilon, int need_epsilon, double right_border, double epsilon_border,
    vector<pair<double, double>>* changes_step, vector<pair<double, double>>* test_task) {
    double xn = x0; //������� �����
    double vn = v0; //������� ��������
    double hn = h;  //������� ���
    double v2n = v0;//���������� ��� �����, ����������� � ���������� �����
    int c1 = 0, c2 = 0, is_task_test = (f == test_f); //c1 - ������� ������� ���� �������, c2 - ������� �������� ����, 
                                                        //is_task_test ������ �� �� ������ �������� ������ 
    //������ ���� ��� �����, ����� ��� ��� �� ���������� ������, � ������ ���������������� ��������
    pair<double, double> half_step_point;//������������� ����� - �������� ����
    pair<double, double> new_point_with_half_step;//�������� �����, ����������� � ���������� �����
    pair<double, double>new_point;//����� ����� ������� ����������� �����

    // ������ ��� �������� �����������
    vector<vector<double>> numerical_solution;
    test_task->push_back({ u_test(xn, v0), abs(u_test(xn, v0) - vn) });
    if (need_epsilon) { // � ��������� ��������� �����������
        vector<double> new_raw = { xn, vn, NAN, NAN, NAN,hn };//�������� ������ ������ �����
        numerical_solution.push_back(new_raw);
        for (int i = 1; i < n_steps; i++) {
            c1 = 0, c2 = 0;
            while (1) {//���� �� ������ ������� �����
                new_point = step_of_the_method_for_equation(f, xn, vn, hn); //��c������ ����� ����� � ������� �����

                half_step_point = step_of_the_method_for_equation(f, xn, vn, hn / 2.0);//������� ����� � ���������� �����
                new_point_with_half_step = step_of_the_method_for_equation(f, half_step_point.first, half_step_point.second, hn / 2.0);
                v2n = new_point_with_half_step.second;

                int olp = check_the_point(new_point.second, new_point_with_half_step.second, &hn, epsilon); // 1 - ����� �������, ��� ��� ��,
                //2 - ����� �������, ��� � ��� ���� ������,
                //0 - ����� ������, ��� � ��� ���� ������
                if (olp == 2) {//������� �������� ����
                    c2 += 1;
                    hn /= 2;
                }
                else if (olp == 0) {//������� ������� ���� �� ���
                    c1 += 1;
                }
                if (olp) {//���� ����� �������
                    xn = new_point.first;
                    vn = new_point.second;//��������� �����
                    double h_to_turn = hn;
                    if (olp == 2) {
                        h_to_turn /= 2;
                    }
                    vector<double> new_raw = { xn, vn, v2n, (vn - v2n), calculate_S(vn,v2n,epsilon),h_to_turn };//��� ��������� �� ���� ���
                    numerical_solution.push_back(new_raw);
                    if (is_task_test) {//��� �������� ������ ��� ������������� �������
                        test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                    }
                    break;
                }
            }
            changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//�������� ������� ��������� ���� ��� ����� ���� 
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)   // ���� �� ��� ��������� � �������-����������� ������ �������, 
                // �� ��� ���� ��������� ����� � �� ��������� ����  
            {
                break;
            }
            if (xn + hn > right_border + epsilon_border) { //���� ��������� ��� ������� ��� �� �����������, �� ������� ��������� ����� �� ������� � ��������� ����
                i++;
                hn = right_border - xn;
                while (1) {//���� �� ������ ������� �����
                    new_point = step_of_the_method_for_equation(f, xn, vn, hn); //��c������ ����� ����� � ������� �����

                    half_step_point = step_of_the_method_for_equation(f, xn, vn, hn / 2.0);//������� ����� � ���������� �����
                    new_point_with_half_step = step_of_the_method_for_equation(f, half_step_point.first, half_step_point.second, hn / 2.0);
                    v2n = new_point_with_half_step.second;

                    int olp = check_the_point(new_point.second, new_point_with_half_step.second, &hn, epsilon); // 1 - ����� �������, ��� ��� ��,
                    //2 - ����� �������, ��� � ��� ���� ������,
                    //0 - ����� ������, ��� � ��� ���� ������
                    if (olp == 2) {//������� �������� ����
                        c2 += 1;
                    }
                    else if (olp == 0) {//������� ������� ���� �� ���
                        c1 += 1;
                    }
                    if (olp) {//���� ����� �������
                        xn = new_point.first;
                        vn = new_point.second;//��������� �����
                        vector<double> new_raw = { xn, vn, v2n, (vn - v2n), calculate_S(vn,v2n,epsilon),hn };//��� ��������� �� ���� ���
                        numerical_solution.push_back(new_raw);
                        if (is_task_test) {//��� �������� ������ ��� ������������� �������
                            test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                        }
                        break;
                    }
                }
                changes_step->push_back({ (*changes_step)[i - 1].first,(*changes_step)[i - 1].second});
            }
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border) { // ���� �� ��� ��������� � �������-����������� ������ �������, 
                                                                                              //�� ��� ���� ��������� ����� � �� ��������� ���� 
                break;
            }
        }
    }
    else { //��� �������� ��������� �����������
        vector<double> new_raw = { xn, vn, hn };//�������� ������ ������ �����
        numerical_solution.push_back(new_raw);
        for (int i = 0; xn < right_border && i < n_steps; i++) {
            new_point = step_of_the_method_for_equation(f, xn, vn, hn);
            xn = new_point.first;
            vn = new_point.second;
            if (is_task_test) {
                new_raw = { xn, vn, hn };//��� ��������� �� ���� ���
                numerical_solution.push_back(new_raw);
                test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
            }
            else {
                vector<double> new_raw = { xn, vn, hn };//��� ��������� �� ���� ���
                numerical_solution.push_back(new_raw);
            }
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border) { // ���� �� ��� ��������� � �������-����������� ������ �������, 
                                                                                              //�� ��� ���� ��������� ����� � �� ��������� ����
                break;
            }
            else if (xn + hn > right_border + epsilon_border) { //���� ��������� ��� ������� ��� �� �����������, �� ������� ��������� ����� �� ������� � ��������� ����
                hn = right_border - xn;
                new_point = step_of_the_method_for_equation(f, xn, vn, hn);
                xn = new_point.first;
                vn = new_point.second;
                if (is_task_test) {
                    new_raw = { xn, vn, hn };//��� ��������� �� ���� ���
                    numerical_solution.push_back(new_raw);
                    test_task->push_back({ u_test(xn,v0),abs(u_test(xn,v0) - vn) });
                }
                else {
                    vector<double> new_raw = { xn, vn, hn };//��� ��������� �� ���� ���
                    numerical_solution.push_back(new_raw);
                }
                break;
            }
        }
    }
    return numerical_solution;
}