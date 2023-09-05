// Simulate Lennard-Jones fluid of Argon using NPT ansambl

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define N 100    // broj čestica
#define Nw 10    // broj šetača
#define Nk 100   // broj koraka
#define Nb 10    // broj blokova
#define Nbskip 0 // broj blokova koje preskačemo radi stabilizacije

float lennardJones_reduced(float x1, float x2, float y1, float y2, float z1, float z2, float sigma)
{
  float Ulj_red;
  float x = sigma / sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2)); // x = sigma/r
  x = pow(x, 6);                                                                     // x = (sigma/r)**6
  if (x <= 1)                                                                        // ignoriraj parove koji su preblizu
    Ulj_red = (pow(x, 2) - x);
  else
    Ulj_red = 0;
  return Ulj_red;
}

float rxLJforce_reduced(float x1, float x2, float y1, float y2, float z1, float z2, float sigma)
{
  float product;
  float x = sigma / sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2)); // x = sigma/r
  x = pow(x, 6);                                                                     // x = (sigma/r)**6
  if (x <= 1)
    product = 24 * (2 * pow(x, 2) - x);
  else
    product = 0;
  return product;
}

float Sum_rf_reduced(float x[Nw + 1][N + 1], float y[Nw + 1][N + 1], float z[Nw + 1][N + 1], float sigma, int i)
{
  float result = 0;
  for (int a = 1; a <= N; a++)
  {
    for (int b = a + 1; b <= N; b++)
    {
      result += rxLJforce_reduced(x[i][a], x[i][b], y[i][a], y[i][b], z[i][a], z[i][b], sigma);
    }
  }
  return result;
}

float Sum_Ulj_reduced(float x[Nw + 1][N + 1], float y[Nw + 1][N + 1], float z[Nw + 1][N + 1], float sigma, int i)
{
  float potEtemp = 0;
  for (int a = 1; a <= N; a++)
  {
    for (int b = a + 1; b <= N; b++)
    {
      potEtemp += lennardJones_reduced(x[i][a], x[i][b], y[i][a], y[i][b], z[i][a], z[i][b], sigma);
    }
  }
  return potEtemp;
}

float Sum_v2(float v[Nw + 1][N + 1], int i)
{
  float sum = 0;
  for (int a = 1; a <= N; a++)
  {
    sum += v[i][a] * v[i][a];
  }
  return sum;
}

int main(void)
{
#pragma region // KONSTANTE & VARIJABLE
  int i, j, k, ib;
  long idum = -1234;
  // konstante:
  float sigma = 1;                                         // 3.4 * 10^(-10) [m]
  float epsilon = 1;                                       // 1.65 * 10^(-21) [J]
  float mass = 1;                                          // 6.69 * 10^(-26) [kg]
  float rho = 1.0481 * pow(10, -3) * mass / pow(sigma, 3); // 1.784 [kg/m^3]
  float L0 = 10 * sigma * cbrt(N / 1.0481);                // L0 = cbrt(N*m/rho)
  float v_avg = 2.8654 * sqrt(epsilon / mass);             // 450 [m/s]
  float k_B = 1;                                           // 1.380649 * pow(10, -23) [m^2*kg/(s^2*K)] // ili računam preko neke T pa dobijem neš na 10^-2?
  // veličine:
  float x[Nw + 1][N + 1]; // broj šetača, broj čestica
  float y[Nw + 1][N + 1]; // broj šetača, broj čestica
  float z[Nw + 1][N + 1]; // broj šetača, broj čestica
  float v[Nw + 1][N + 1]; // brzine čestica
  float L[Nw + 1];        // stranica
  float V[Nw + 1];        // volumen
  float U[Nw + 1];        // unutarnja energija (Upot + Ukin)
  float Uk[Nw + 1];       // kinetička energija
  float Up[Nw + 1];       // potencijalna energija
  float T[Nw + 1];        // temperatura
  float press[Nw + 1];    // tlak
  // promjene
  float x0, y0, z0, v0, vmin = 0.9 * v_avg, vmax = 1.1 * v_avg; // prosječna brzina čestice u plinu je oko 450 m/s
  float dL, dLMax = L0 / 100;                                   // metara
  float dx, dy, dz, dv;                                         // promjene koordinata i brzina
  float dxyzMax = L0 / 100, dvMax = v_avg / 100;                // maksimalne promjene
  // pomoćne varijable
  float delta_V, delta_H, delta_U, delta_W;
  float x_old[N + 1], y_old[N + 1], z_old[N + 1], v_old[N + 1]; // za reversat promjenu koordinata ako se korak odbacuje
  float V_old, L_old, U_old, Up_old, Uk_old, T_old, press_old;  // za reversat promjenu u slučaju odbacivanja
  float p, ran;                                                 // vjerojatnost prihvaćanja promjene, random broj pri određivanju odbacivanja/prihvaćanja koraka
  int N_change_volume = 5;
  float L_mean, V_mean, U_mean, p_mean, T_mean; // vrijednosti u koraku usrednjene po svim šetačima
  int accepted = 0, rejected = 0;
  float ratio;
  float block_avg_L, block_avg_V, block_avg_U, block_avg_T, block_avg_p; // usrednjivanje po bloku
#pragma endregion

  FILE *data, *block_data;
  data = fopen("data.txt", "w");
  block_data = fopen("block_data.txt", "w");

  // inicijalizacija čestica
  for (i = 1; i <= Nw; i++) // po šetačima
  {
    for (j = 1; j <= N; j++) // po česticama
    {
      // dobro je da su sve čestice na početku udaljene barem za σ.
      x0 = ran1(&idum) * L0; // ran1(&idum) = random broj iz[0, 1]
      y0 = ran1(&idum) * L0; // ran1(&idum) = random broj iz[0, 1]
      z0 = ran1(&idum) * L0; // ran1(&idum) = random broj iz[0, 1]
      v0 = ran1(&idum) * (vmax - vmin) + vmin;
      x[i][j] = x0; // broj šetača, broj čestica
      y[i][j] = y0; // broj šetača, broj čestica
      z[i][j] = z0; // broj šetača, broj čestica
      v[i][j] = v0;
    }
    // računamo unutarnju energiju početne konfiguracije
    L[i] = L0;
    V[i] = L0 * L0 * L0;
    Uk[i] = 0.5 * mass * Sum_v2(v, i);                        // Ukin = Sum(i=1,...,N) 1/2*mv_i^2
    Up[i] = 4 * epsilon * Sum_Ulj_reduced(x, y, z, sigma, i); // Upot = Sum(i=1,...,N)Sum(j=i+1,...,N) 4*epsilon*((sigma/r_ij)^12-(sigma/r_ij)^6)
    U[i] = Uk[i] + Up[i];                                     // U_ukupna = Ukin + Upot
    T[i] = 0.6666 * Uk[i] / k_B / N;                          // T = 2/3 Ukin /(kB*N)
    press[i] = 1 / V[i] * (N * k_B * T[i] - epsilon / (3 * k_B * T[i]) * Sum_rf_reduced(x, y, z, sigma, i));
    printf("T=%f; p=%f; Upot=%f; Ukin=%f; Uuk=%f; V=%f\n", T[i], press[i], Up[i], Uk[i], U[i], V[i]);
  }

  for (ib = 1; ib <= Nb; ib++) // po bloku
  {
    block_avg_L = 0;
    block_avg_V = 0;
    block_avg_U = 0;
    block_avg_T = 0;
    block_avg_p = 0;
    for (i = 1; i <= Nk; i++) // po koracima
    {
      for (j = 1; j <= Nw; j++) // po šetačima - mikrostanja - šetači po 3N dimenzionalnom faznom prostoru
      {
        L_old = L[j];
        V_old = V[j];
        Uk_old = Uk[j];
        Up_old = Up[i];
        U_old = U[j];
        T_old = T[j];
        press_old = press[j];
        for (k = 1; k <= N; k++) // po česticama
        {
          dx = (ran1(&idum) * 2 - 1) * dxyzMax;
          dy = (ran1(&idum) * 2 - 1) * dxyzMax;
          dz = (ran1(&idum) * 2 - 1) * dxyzMax;
          dv = (ran1(&idum) * 2 - 1) * dvMax;
          x_old[k] = x[j][k];
          y_old[k] = y[j][k];
          z_old[k] = z[j][k];
          v_old[k] = v[j][k];
          x[j][k] += dx;
          y[j][k] += dy;
          z[j][k] += dz;
          v[j][k] += dv;
          // brzina ne smije biti negativna (jer odražava apsolutnu vrijednost)
          if (v[j][k] < 0)
            v[j][k] = -v[j][k];
          // rubni uvjeti
          if (x[j][k] < 0)
            x[j][k] = -x[j][k];
          if (x[j][k] > L[j])
            x[j][k] = L[j] - (x[j][k] - L[j]);
          if (y[j][k] < 0)
            y[j][k] = -y[j][k];
          if (y[j][k] > L[j])
            y[j][k] = L[j] - (y[j][k] - L[j]);
          if (z[j][k] < 0)
            z[j][k] = -z[j][k];
          if (z[j][k] > L[j])
            z[j][k] = L[j] - (z[j][k] - L[j]);
        }

        if (i % N_change_volume == 0) // SVAKI N-ti korak promijeni i volumen...
        {
          dL = (ran1(&idum) * 2 - 1) * dLMax;
          L[j] = L[j] + dL;
          if (L[j] < 0) // duljina stranice ne može biti negativna
            L[j] = -L[j];
          V[j] = L[j] * L[j] * L[j];
          for (k = 1; k <= N; k++) // po česticama
          {
            x[j][k] = L[j] / L_old * x[j][k]; // skaliranje svih x koordinata faktorom L'/L
            y[j][k] = L[j] / L_old * y[j][k]; // skaliranje svih y koordinata faktorom L'/L
            z[j][k] = L[j] / L_old * z[j][k]; // skaliranje svih z koordinata faktorom L'/L
          }
        }
        // ponovno računanje termodinamičkih vrijednosti
        Uk[j] = 0.5 * mass * Sum_v2(v, j);                        // Ukin = Sum(i=1,...,N) 1/2*mv_i^2
        Up[j] = 4 * epsilon * Sum_Ulj_reduced(x, y, z, sigma, j); // Upot = Sum(i=1,...,N)Sum(j=i+1,...,N) 4*epsilon*((sigma/r_ij)^12-(sigma/r_ij)^6)
        U[j] = Uk[j] + Up[j];                                     // U_ukupna = Ukin + Upot
        T[j] = 0.6666 * Uk[j] / k_B / N;                          // T = 2/3 Ukin /(kB*N)
        press[j] = 1 / V[j] * (N * k_B * T[j] - epsilon / (3 * k_B * T[j]) * Sum_rf_reduced(x, y, z, sigma, j));
        // delte
        delta_V = V[j] - V_old;
        delta_U = U[j] - U_old;
        delta_H = delta_U + press[j] * delta_V - k_B * T[j] * N * log(V[j] / V_old);
        delta_W = 1 / (k_B * T[j]) * delta_H;
#pragma region // Metropolis algoritam
        if (delta_W > 0)
        {
          p = exp(-delta_W); // vjerojatnost prihvaćanja promjene
          ran = ran1(&idum);
          if (ran < p)
          {
            accepted++;
            ratio = (float)accepted / (float)(accepted + rejected);
            // printf("%f\n", ratio);
          }
          else // odbacivanje promjena
          {
            rejected++;
            ratio = (float)accepted / (float)(accepted + rejected);
            // printf("%f\n", ratio);
            // resetiranje vrijednosti
            L[i] = L_old;
            V[j] = V_old;
            Uk[j] = Uk_old;
            Up[j] = Up_old;
            U[j] = U_old;
            T[j] = T_old;
            press[j] = press_old;
            for (k = 1; k <= N; k++) // po česticama
            {
              x[j][k] = x_old[k];
              y[j][k] = y_old[k];
              z[j][k] = z_old[k];
              v[j][k] = v_old[k];
            }
          }
        }
        else
        {
          accepted++;
          ratio = (float)accepted / (float)(accepted + rejected);
          // printf("%f\n", ratio);
        }
#pragma endregion
        // kraj petlje šetača
      }
      // // Maksimalnu duljinu koraka podesavamo kako bi prihvacanje bilo oko 50%
      if (ratio > 0.5)
      {
        dxyzMax = dxyzMax * 1.05;
        dvMax = dvMax * 1.05;
        if (i % N_change_volume == 0)
          dLMax = dLMax * 1.05;
      }
      if (ratio < 0.5)
      {
        dxyzMax = dxyzMax * 0.95;
        dvMax = dvMax * 0.95;
        if (i % N_change_volume == 0)
          dLMax = dLMax * 0.95;
      }
      // za svaki korak zapisujemo parametre, usrednjene po šetačima
      if (ib >= Nbskip) // ako je završena stabilizacija
      {
#pragma region // CALCULATING THE MEANS
        L_mean = 0;
        V_mean = 0;
        U_mean = 0;
        T_mean = 0;
        p_mean = 0;
        for (j = 1; j <= Nw; j++)
        {
          L_mean += L[j];
          V_mean += V[j];
          U_mean += U[j];
          T_mean += T[j];
          p_mean += press[j];
        }
        L_mean = L_mean / Nw;
        V_mean = V_mean / Nw;
        U_mean = U_mean / Nw;
        T_mean = T_mean / Nw;
        p_mean = p_mean / Nw;
        fprintf(data, "%d\t%f\t%f\t%f\t%f\t%f\n", ib * Nk + i, L[1], V_mean, U_mean, T_mean, p_mean);
        // za svaki korak dodaj srednju vrijednost šetača blok average sumi
        block_avg_L += L_mean;
        block_avg_V += V_mean;
        block_avg_U += U_mean;
        block_avg_T += T_mean;
        block_avg_p += p_mean;
#pragma endregion
      }
      // kraj korak levela
    }
    if (ib >= Nbskip) // ako je završena stabilizacija
    {
      block_avg_L = block_avg_L / Nk;
      block_avg_V = block_avg_V / Nk;
      block_avg_U = block_avg_U / Nk;
      block_avg_T = block_avg_T / Nk;
      block_avg_p = block_avg_p / Nk;
      fprintf(block_data, "%d\t%f\t%f\t%f\t%f\t%f\n", ib * Nk, block_avg_L, block_avg_V, block_avg_U, block_avg_T, block_avg_p);
      // printf("%d\t%f\t%f\t%f\t%f\t%f\n", ib, block_avg_L, block_avg_V, block_avg_U, block_avg_T, block_avg_p);
    }
  }
  printf("Final ratio: %f\n", ratio);

  fclose(data);
  fclose(block_data);
  return 0;
}
