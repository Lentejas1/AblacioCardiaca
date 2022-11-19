#include <iostream>
#include <cstdio>
#include <cmath>
#define N 100
#define M 2000
using namespace std;

double Tc = 273.15 + 36.5;
double V = 40;
double cond = 0.472, k = 0.56, cv = 3683, rho = 1081;
double P = cond * pow(V, 2) / 4;
double dz = 1. / N, ta = 0.025;


double dt(float ratio){
    return ratio * pow(dz, 2);
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
            T[i+1][j] = dt(ratio) / (pow(dz, 2)) * (T[i][j+1] - 2 * T[i][j] + T[i][j-1]) + dt(ratio) + T[i][j];
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
    double gamma = dt(ratio) / pow(dz,2);

    for (int j = 0; j < N; j++){
        T[0][j] = Tc * k / P;
    }
    for (int i = 0; i < M; i++){
        T[i][0] = Tc * k / P;
        T[i][N-1] = Tc * k / P;
    }

    double T_gs[N];
    for (int i = 0; i + 1 < M; i++) {
        for (int g; g < 400; g++){
            for (int j = 1; j<N; j++){
                T_gs[j] = T[i][j-1];
            }
            for (int j = 1; j < N; j++) {
                T[i][j] = ((T_gs[j+1]+T_gs[j-1])*gamma + dt(ratio) + T[i-1][j])/(1+2*gamma);
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
        explicit_method(0.49);
        implicit_method(0.49);
        return 0;
}
