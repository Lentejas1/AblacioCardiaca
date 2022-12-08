#include <iostream>
#include <cstdio>
#include <cmath>
#define N 100
#define M 2000
using namespace std;

double Tc = 273.15 + 36.5;
double V = 40;
double cond = 0.472, k = 0.56;
double P = cond * pow(V, 2) / 2;
double dz = 1. / N;
double alpha = 0.56 / (3683 * 1081);
double tn = pow(0.02,2) / alpha;

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
        t += 0.0001;
        T_f = f(0.375,t);
    }
    return t*tn;
}

void explicit_method(float ratio){
    double T[M][N];
    for (int j = 0; j < N; j++){
        T[0][j] = Tc * k / P;
    }
    for (int i = 0; i < M; i++){
        T[i][0] = Tc * k / P;
        T[i][N-1] = Tc * k / P;
    }


    for (int i = 0; i + 1 < M; i++) {
        for (int j = 1; j < N; j++) {
            T[i+1][j] =
                    dt_exp(ratio) / (pow(dz, 2)) * (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) + dt_exp(ratio) + T[i][j];
        }
    }


    FILE* txt;

    txt = fopen("Ablacio_Explicit.txt", "w");

    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            fprintf(txt, "%lf,", T[i][j]*P/k);
        }
        fprintf(txt, "%lf\n", T[i][N-1]*P/k);
    }

    fclose(txt);

}

void implicit_method(float ratio){
    double T[M][N];

    for (int j = 0; j < N; j++){
        T[0][j] = Tc * k / P;
    }
    for (int i = 0; i < M; i++){
        T[i][0] = Tc * k / P;
        T[i][N-1] = Tc * k / P;
    }

    double T_gs[N];
    for (int i = 0; i + 1 < M; i++) {
        for (int g = 1; g < 400; g++){
            for (int j = 0; j<N; j++){
                T_gs[j] = T[i][j];
            }
            for (int j = 1; j < N; j++) {
                T[i][j] = ((T_gs[j+1]+T_gs[j-1])*ratio + dt_imp(ratio) + T[i - 1][j]) / (1 + 2 * ratio);
            }
        }

    }


    FILE* txt;

    txt = fopen("Ablacio_Implicit.txt", "w");

    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            fprintf(txt, "%lf,", T[i][j]*P/k);
        }
        fprintf(txt, "%lf\n", T[i][N-1]*P/k);
    }

    fclose(txt);
}



int main(){
    explicit_method(0.51);
    //implicit_method(0.49);
    cout << analytic_sol();
    return 0;
}