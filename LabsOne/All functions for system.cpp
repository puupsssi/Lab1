#include"Names of Functions.h"
// ������� f(x, v) ��� ������� ���������


pair<double, double> f_second_task(double x, double u_1, double u_2,  double a, double b) {
    double du1_dx = u_2;
    double du2_dx = -a*u_2-b*sin(u_1);

    return { du1_dx, du2_dx };
}


tuple<double, double, double> step_of_the_method_for_the_system(pair<double, double>(*f)(double, double, double, double, double), double xn, double v_1n, double v_2n, double h, double a, double b) {
    // ���������� ������������� k1, k2, k3, k4 ��� u1 � u2 - ��� ������ ����� ����� 4�� �������
    double k1_u1 = f(xn, v_1n, v_2n, a, b).first;
    double k1_u2 = f(xn, v_1n, v_2n, a, b).second;

    double k2_u1 = f(xn + h / 2.0, v_1n + h / 2.0 * k1_u1, v_2n + h / 2.0 * k1_u2, a, b).first;
    double k2_u2 = f(xn + h / 2.0, v_1n + h / 2.0 * k1_u1, v_2n + h / 2.0 * k1_u2, a, b).second;

    double k3_u1 = f(xn + h / 2.0, v_1n + h / 2.0 * k2_u1, v_2n + h / 2.0 * k2_u2, a, b).first;
    double k3_u2 = f(xn + h / 2.0, v_1n + h / 2.0 * k2_u1, v_2n + h / 2.0 * k2_u2, a, b).second;

    double k4_u1 = f(xn + h, v_1n + h * k3_u1, v_2n + h * k3_u2, a, b).first;
    double k4_u2 = f(xn + h, v_1n + h * k3_u1, v_2n + h * k3_u2, a, b).second;

    // ���������� �������� ���������� u1 � u2
    xn = xn + h;  // ��������� x �� ���� ���
    v_1n += (h / 6.0) * (k1_u1 + 2 * k2_u1 + 2 * k3_u1 + k4_u1);
    v_2n += (h / 6.0) * (k1_u2 + 2 * k2_u2 + 2 * k3_u2 + k4_u2);

    return { xn, v_1n, v_2n }; // ���������� ����������� ��������
}


double calculate_S_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double epsilon) {
    double S_1 = abs(v_12n - v_1n_plus_1) / 15.0; //2^4 - 1 = 15;
    double S_2 = abs(v_22n - v_2n_plus_1) / 15.0;
    double S = sqrt(pow(S_1, 2) + pow(S_2, 2));//����������� �������� ��� �������� ��������� �����������
    return S;
}


int check_the_point_for_system(double v_1n_plus_1, double v_12n, double v_2n_plus_1, double v_22n, double* h, double epsilon) {
    double S = calculate_S_for_system(v_1n_plus_1, v_12n, v_2n_plus_1, v_22n, epsilon);
    if ((epsilon / 32.0) <= S && S < epsilon) { //������� �����
        return 1; //������, ��� ����� ������� � �� ���������� ����
    }
    if (S < (epsilon / 32.0)) {
        *h = *h * 2;
        return 2; //���������� � ���������� �������, ������ ��� �������� ���
    }
    if (S >= epsilon) {//����� ������, ������ ��� � �������� ��� �� ����� (xn,v_1,v_2)
        *h = *h / 2.0;
        return 0;
    }
    return 0;
}


vector<vector<double>> runge_kutta_4th_order_for_system(pair<double, double>(*f)(double, double, double, double, double), double x0, double v_10, double v_20, double a, double b,//������� �������, ������� �������� �� ������ �����
    double h, int n_steps, double epsilon, int need_epsilon, double right_border, double epsilon_border,
    vector< pair<double, double>>* changes_step) {
    double xn = x0;//������� �����
    double v_1n = v_10;//������� �������� ������ ����������
    double v_2n = v_20;//������� �������� ������ ����������
    double hn = h; //������� ���
    double v_12n, v_22n;//���������� ��� �����, ����������� � ���������� �����
    int c1 = 0, c2 = 0; // c1 - ������� ������� ���� �������, c2 - ������� �������� ����
    //������ ����� ��� �����, ����� ��� ��� �� ���������� ������, � ������ ���������������� ��������
    tuple<double, double, double> half_step_point;//������������� ����� - �������� ����
    tuple<double, double, double> new_point_with_half_step;//�������� �����, ����������� � ���������� �����
    tuple<double, double, double> new_point;//����� ����� ������� ����������� �����

    // ������ ��� �������� �����������
    vector<vector<double>> numerical_solution;

    if (need_epsilon) { // � ��������� ��������� �����������
        vector<double> new_raw = { xn, v_1n,v_2n, NAN, NAN, NAN,NAN,NAN,hn };//�������� ������ ������ �����
        numerical_solution.push_back(new_raw);
        for (int i = 1; i < n_steps + 1; i++) {
            while (1) {//���� �� ������ ������� �����
                c1 = 0, c2 = 0;
                new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b); //��c������ ����� ����� � ������� �����

                half_step_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn / 2.0, a, b);//������� ����� � ���������� �����
                new_point_with_half_step = step_of_the_method_for_the_system(f, get<0>(half_step_point), get<1>(half_step_point), get<2>(half_step_point), hn / 2.0, a, b);
                v_12n = get<1>(new_point_with_half_step);
                v_22n = get<2>(new_point_with_half_step);

                int olp = check_the_point_for_system(get<1>(new_point), get<1>(new_point_with_half_step),
                                                        get<2>(new_point), get<2>(new_point_with_half_step), &hn, epsilon);
                // 1 - ����� �������, ��� ��� ��,
                //2 - ����� �������, ��� � ��� ���� ������,
                //0 - ����� ������, ��� � ��� ���� ������
                if (olp == 2) {//������� �������� ����
                    c2 += 1;
                }
                else if (olp == 0) {//������� ������� ���� �� ���
                    c1 += 1;
                }

                if (olp) {//���� ����� �������
                    xn = get<0>(new_point);
                    v_1n = get<1>(new_point);//��������� �����
                    v_2n = get<2>(new_point);//��������� �����
                    double h_to_turn = hn;
                    if (olp == 2) {
                        h_to_turn /= 2;
                    }
                    vector<double> new_raw = { xn, v_1n,v_2n, v_12n,v_22n, (v_1n - v_12n),(v_2n - v_22n),32*calculate_S_for_system(v_1n,v_1n,v_12n,v_22n,epsilon),h_to_turn };//��� ��������� �� ���� ���
                    numerical_solution.push_back(new_raw);
                    break;
                }
            }
            changes_step->push_back({ (*changes_step)[i - 1].first + c1,(*changes_step)[i - 1].second + c2 });//�������� ������� ��������� ���� ��� ����� ����
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)   // ���� �� ��� ��������� � �������-����������� ������ �������, 
                // �� ��� ���� ��������� ����� � �� ��������� ����  
            {
                break;
            }
            if (xn + hn > right_border + epsilon_border) //���� ��������� ��� ������� ��� �� �����������, �� ������� ��������� ����� �� ������� � ��������� ���� 
                                                        //(���� ������ ��� �� ����������� ���������, ��� �� ���� �� ������� � ������� ����������� - � ���� ������ ������ ���������� �������)
            {
                i++;
                hn = right_border - xn;
                while (1) {
                    c1 = 0, c2 = 0;
                    new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b); //��c������ ����� ����� � ������� �����

                    half_step_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn / 2.0, a, b);//������� ����� � ���������� �����
                    new_point_with_half_step = step_of_the_method_for_the_system(f, get<0>(half_step_point), get<1>(half_step_point), get<2>(half_step_point), hn / 2.0, a, b);
                    v_12n = get<1>(new_point_with_half_step);
                    v_22n = get<2>(new_point_with_half_step);

                    int olp = check_the_point_for_system(get<1>(new_point), get<1>(new_point_with_half_step),
                        get<2>(new_point), get<2>(new_point_with_half_step), &hn, epsilon);
                    // 1 - ����� �������, ��� ��� ��,
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
                        xn = get<0>(new_point);
                        v_1n = get<1>(new_point);//��������� �����
                        v_2n = get<2>(new_point);//��������� �����
                        vector<double> new_raw = { xn, v_1n,v_2n, v_12n,v_22n, (v_1n - v_12n),(v_2n - v_22n),32*calculate_S_for_system(v_1n,v_1n,v_12n,v_22n,epsilon),hn };//��� ��������� �� ���� ���
                        numerical_solution.push_back(new_raw);
                        break;
                    }
                }
                changes_step->push_back({ (*changes_step)[i - 1].first,(*changes_step)[i - 1].second });//�������� ������� ��������� ���� ��� ����� ����
            }
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)   // ���� �� ��� ��������� � �������-����������� ������ �������, 
                                                                                              // �� ��� ���� ��������� ����� � �� ��������� ����  
            {
                break;
            }
        }
    }
    else { //��� �������� ��������� �����������
        vector<double> new_raw = { xn, v_1n,v_2n, hn };//�������� ������ ������ �����
        numerical_solution.push_back(new_raw);
        for (int i = 0;xn < right_border && i < n_steps; i++) {
            //�������� ����� �� �������
            new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b);
            xn = get<0>(new_point);
            v_1n = get<1>(new_point);//��������� �����
            v_2n = get<2>(new_point);//��������� �����
            vector<double> new_raw = { xn, v_1n,v_2n, hn };//��� ��������� �� ���� ���
            numerical_solution.push_back(new_raw);
            if (right_border - epsilon_border <= xn && xn <= right_border + epsilon_border)   // ���� �� ��� ��������� � �������-����������� ������ �������, 
                                                                                              //�� ��� ���� ��������� ����� � �� ��������� ����  
            {
                break;
            }
            else if (xn + hn > right_border + epsilon_border) //���� ��������� ��� ������� ��� �� �����������, �� ������� ��������� ����� �� ������� � ��������� ����
            {
                hn = right_border - xn;
                new_point = step_of_the_method_for_the_system(f, xn, v_1n, v_2n, hn, a, b);
                xn = get<0>(new_point);
                v_1n = get<1>(new_point);//��������� �����
                v_2n = get<2>(new_point);//��������� �����
                vector<double> new_raw = { xn, v_1n,v_2n, hn };//��� ��������� �� ���� ���
                numerical_solution.push_back(new_raw);
                break;
            }
        }
    }
    return numerical_solution;
}