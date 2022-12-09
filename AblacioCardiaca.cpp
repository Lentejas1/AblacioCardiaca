#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;

double Tc = 273.15 + 36.5;
double V = 40;
double cond = 0.472, k = 0.56;
double P = cond * pow(V, 2) / 2;
double dz = 1. / 100;
double alpha = 0.56 / (3683 * 1081);
double tn = pow(0.02,2) / alpha;
double cosa= 0.5*dz;
int N = 101;
double dt_exp(float ratio){
    return ratio * pow(dz, 2);
}

double dt_imp(float ratio){
    return ratio * dz;
}

double f(double x, double t){
    double sum = 0;
    for (int n = 0; n < 801; n++){
        if (n % 2 != 0) {
            sum += 4 / (n * M_PI) * (1 - exp(-1 * t * pow(n * M_PI, 2))) / pow(n * M_PI, 2) * sin(
                    n * M_PI * x);
        }
    }
    return sum + Tc * k / P;
}

double analytic_sol(){
    double T_f, T_opt;
    T_f = (36.5 + 273.15) * k / P;
    T_opt = (50 + 273.15) * k / P;
    double t = 0;
    while (T_f < T_opt){
        t += cosa;
        T_f = f(0.375,t);
    }
    return t*tn;
}

void explicit_method(double ratio, int M){
    double T[M][N];
    for (int j = 0; j < N; j++){
        T[0][j] = Tc * k / P;
    }
    for (int i = 0; i < M; i++){
        T[i][0] = Tc * k / P;
        T[i][N-1] = Tc * k / P;
    }


    for (int i = 0; i + 1 < M; i++) {
        for (int j = 1; j < N -1; j++) {
            T[i+1][j] =
                    dt_exp(ratio) / (pow(dz, 2)) * (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) + dt_exp(ratio) + T[i][j];
        }
    }


    FILE* txt;

    txt = fopen("Ablacio_Explicit051.txt", "w");

    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            fprintf(txt, "%lf,", T[i][j]*P/k);
        }
        fprintf(txt, "%lf\n", T[i][N-1]*P/k);
    }

    fclose(txt);

}

double gamma(double ratio){
    return ratio / dz;
}

void implicit_method(float ratio, int M){
    double T[M][N];

    for (int j = 0; j < N; j++){
        T[0][j] = Tc * k / P;
    }
    for (int i = 0; i < M; i++){
        T[i][0] = Tc * k / P;
        T[i][N-1] = Tc * k / P;
    }


    for (int i = 1; i < M; i++) {
        int j;
        T[i][j] = T[i-1][j];

        for (int g = 1; g < 2000; g++){
            for (int j = 1; j+1<N; j++){
                T[i][j] = ((T[i][j+1]+T[i][j-1])*gamma(ratio) + dt_imp(ratio) + T[i-1][j]) / (1 + 2 * gamma(ratio));
            }
        }
    }

    FILE* txt;
    txt = fopen("Ablacio_Implicit100.txt", "w");

    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            fprintf(txt, "%lf,", T[i][j]*P/k);
        }
        fprintf(txt, "%lf\n", T[i][N-1]*P/k);
    }

    fclose(txt);
}



int main(){
    //explicit_method(0.25, 1001);
    //explicit_method(0.49, 511);
    explicit_method(0.51, 492);
    //implicit_method(0.5, 6);
    //implicit_method(1, 3);
    return 0;
}