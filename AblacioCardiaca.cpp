#include <iostream>
#include <cstdio>
#include <cmath>
#define N 100

using namespace std;

double T[N][N];
double T_0 = 36.5;
double V = 80/sqrt(2);
double cond = 0.472, z = 2, S = 0.1, k = 0.56, cv = 3683, rho = 1081, P = S*cond/z* pow(V, 2);
double l = sqrt(k/P);
double dz = z / N, dt = 0.51* pow(dz, 2), ta = 0.025;

int main(){
    for (int i = 0; i < N; i++){
        T[i][0] = 1;
        T[0][i] = 1;
        T[N-1][0] = 1;
    }


    for (int i = 1; i + 1 < N; i++) {
        for (int j = 1; j < N; j++) {
            T[i+1][j] = dt / (pow(dz, 2)) * (T[i][j+1] - 2 * T[i][j] + T[i][j-1]) + dt + T[i][j];
        }
    }


    FILE* txt;

    txt = fopen("Ablacio.txt", "w");

    for (int j = 0; j < N; j++){
        for (int i = 0; i < N - 1; i++){
            fprintf(txt, "%lf,", T[i][j]*T_0);
        }
        fprintf(txt, "%lf\n", T[N-1][j]);
    }

    fclose(txt);
}
