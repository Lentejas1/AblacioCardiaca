#include <iostream>
#include <cstdio>
#include <cmath>
#define N 100
#define M 2000
using namespace std;

double T[M][N];
double Tc = 273.15 + 36.5;
double V = 40;
double cond = 0.472, k = 0.56, cv = 3683, rho = 1081;
double P = cond * pow(V, 2) / 4;
double dz = 1. / N, dt = 0.49 * pow(dz, 2), ta = 0.025;

int main(){
        for (int j = 0; j < N; j++){
            T[0][j] = Tc * k / P;
            }
        for (int i = 0; i < M; i++){
            T[i][0] = Tc * k / P;
            T[i][N-1] = Tc * k / P;
        }


    for (int i = 0; i + 1 < M; i++) {
        for (int j = 1; j < N; j++) {
            T[i+1][j] = dt / (pow(dz, 2)) * (T[i][j+1] - 2 * T[i][j] + T[i][j-1]) + dt + T[i][j];
        }
    }


    FILE* txt;

    txt = fopen("Ablacio.txt", "w");

    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            fprintf(txt, "%lf,", T[i][j]*P/k);
        }
        fprintf(txt, "%lf\n", T[i][N-1]*P/k);
    }

    fclose(txt);
}
