// Solves the 1-dimensional Euler equations (inviscid Navier-Stokes equations), using MacCormack's method (2nd order)
// for the case of a Sod shock-tube

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/******************************************************************************/
// CONSTANTES Y VARIABLES GLOBALES

// Parámetros constantes de la simulación
const int    NX = 400;         // Tamaño de la malla
const double XL = 0.0;         // Coordenada física del extremo izquierdo
const double XR = 1.0;         // Coordenada física del extremo derecho
const double TFIN = 0.2;       // Tiempo final de integración
const double CFL = 0.9;        // Parametro de Courant
const double dtprint = 0.01;   // Intervalo para escribir a disco

const double gam=1.4;

const int ieq=3;

const double eta = 0.1;

// Para tubo de choque: coord de la separación entre estados iniciales
const double X0 = 0.5;

// Para ecuación de advección: velocidad de ondas ('a' en las notas)
const double A = 1.0;

// Constantes derivadas de las anteriores
const double DX = (XR-XL)/NX;      // Espaciamiento de la malla

// Variables globales
double U[ieq][NX+2];       // Variables conservadas actuales
double UP[ieq][NX+2];      // Variables conservadas "avanzadas"
double F[ieq][NX+2];       // Flujos físicos
double P[ieq][NX+2];
double UT[ieq][NX+2];
double dt;            // Paso de tiempo
double t;          // Tiempo actual
int it;               // Iteración actual
clock_t start;        // Tiempo de inicio
double tprint;        // Tiempo para el siguiente output
int itprint;          // Número de salida


/******************************************************************************/

// Impone las condiciones iniciales
void initflow(double U[ieq][NX+2]) {

  // Inicializar los valores de U en todo el dominio
  // Nótese que también llenamos las celdas fantasma
double x;
const double rhol=1.0;
const double vl=0.0;
const double pl=1.0;
const double rhor=0.125;
const double vr=0.0;
const double pr=0.1;

  for (int i=0; i <= NX+1; i++){
    x=XL+i*DX;
    if (x<=X0){
      U[0][i]=rhol;
      U[1][i]=rhol*vl;
      U[2][i]=0.5*rhol*vl*vl+pl/(gam-1);
    } else {
      U[0][i]=rhor;
      U[1][i]=rhor*vr;
      U[2][i]=0.5*rhor*vr*vr+pr/(gam-1);
    }
  }

  // Inicializar otras variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;

}

/******************************************************************************/

// Escribe a disco el estado de la simulación
void output(double U[ieq][NX+2]) {

  // Generar el nombre del archivo de salida
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // Abrir el archivo
  fstream fout(fname, ios::out);

  // Escribir los valores de U al archivo
  double x;
  for (int i=0; i <= NX; i++) {
    x = XL + i*DX;
    fout << x << " " << P[0][i] << " " << P[1][i] << " " << P[2][i] << " " << endl;
  }

  // Cerrar archivo
  fout.close();

  printf("Se escribió salida %s\n", fname);

  // Avanzar variables de output
  itprint = itprint + 1;
  tprint = itprint * dtprint;

}

/******************************************************************************/

// Aplicar condiciones de frontera a celdas fantasma
// El arreglo pasado es al que aplicaremos las BCs
void boundary(double U[ieq][NX+2]) {

  for(int iieq=0; iieq<=2;iieq++){
    U[iieq][0]=U[iieq][1];
    U[iieq][NX+1]=U[iieq][NX];
  }
}

/******************************************************************************/

// Calcular los flujos físicos F !tambien se incluyen las celdas fantasma
void primitivas(double U[ieq][NX+2], double P[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {
    //F[i]=A*U[i];
    P[0][i]=U[0][i];
    P[1][i]=U[1][i]/U[0][i];
    P[2][i]=(gam-1)*(U[2][i]-pow(U[1][i],2)/(2*U[0][i]));

  }

}


// Calcular los flujos físicos F !tambien se incluyen las celdas fantasma
void fluxes(double P[ieq][NX+2], double F[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {

    F[0][i]=P[0][i]*P[1][i];
    F[1][i]=P[0][i]*pow(P[1][i],2)+P[2][i];
    F[2][i]=P[1][i]*(0.5*P[0][i]*pow(P[1][i],2)+P[2][i]*gam/(gam-1));
  }

}

/******************************************************************************/

// Calcula el paso de tiempo resultante de la condición CFL
double timestep(double P[ieq][NX+2]) {

  double dt;

  // Para advección, max_u es simplemente A
  //double max_speed = abs(A);

  // Para otros casos, debe calcular el valor máximo de |velocidad|
  double max_speed = 0.0;
  double cs;
  double k;
  for (int i = 1; i <= NX; i++) {
    cs = sqrt(gam*P[2][i]/P[0][i]);
    k = abs(P[1][i])+cs;
    if (k > max_speed) max_speed = k;

}

  dt = CFL * DX / max_speed;
printf("dt = %f, t = %f, it = %i \n",dt, t, it);
  return dt;

}

/******************************************************************************/

// Aplica el método de Lax para obtener las UP a partir de las U
// Supone que los flujos F ya fueron actualizados
void mac(double P[ieq][NX+2], double UT[ieq][NX+2], double U[ieq][NX+2], double F[ieq][NX+2], double UP[ieq][NX+2]) {

  for (int i=1; i<=NX; i++){
    for (int iieq=0; iieq<=2; iieq++){
      UT[iieq][i]=U[iieq][i]-dt/DX*(F[iieq][i+1]-F[iieq][i]);
      //UP[iieq][i]=(U[iieq][i+1]+U[iieq][i-1])/2.0-dt/(2*DX)*(F[iieq][i+1]-F[iieq][i-1]);
    }
  }

  boundary(UT);

  primitivas(UT,P);

  fluxes(P,F);

  for (int i=1; i<=NX; i++){
    for (int iieq=0; iieq<=2; iieq++){
      UP[iieq][i]=(U[iieq][i]+UT[iieq][i])/2.0-dt/(2.0*DX)*(F[iieq][i]-F[iieq][i-1]);
      //UP[iieq][i]=(U[iieq][i+1]+U[iieq][i-1])/2.0-dt/(2*DX)*(F[iieq][i+1]-F[iieq][i-1]);
    }
  }

}

/******************************************************************************/

// Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
void stepviscoso(double U[ieq][NX+2], double UP[ieq][NX+2]) {

  for (int i = 1; i <= NX; i++) {
    for (int iieq=0; iieq<=2; iieq++){
      //U[iieq][i]=UP[iieq][i];
      if (((UP[iieq][i+1]-UP[iieq][i])*(UP[iieq][i]-UP[iieq][i-1]))<0) {
        U[iieq][i]=UP[iieq][i]+eta*(UP[iieq][i+1]+UP[iieq][i-1]-2.0*UP[iieq][i]);
      }
      else {
        U[iieq][i]=UP[iieq][i];
      }

    }
  }

  t = t + dt;
  it = it + 1;

}

/******************************************************************************/

int main() {

  // Condición inicial e inicializaciones
  initflow(U);

  primitivas(U,P);

  // Escribir condición inicial a disco
  output(U);

  // Tiempo de inicio de la simulación
  start = clock();
  while (t <= TFIN) {

    primitivas(U,P); // actualiza primitivas

    // Actualizar el paso de tiempo
    dt = timestep(P);

    // Actualizar flujos físicos
    fluxes(P, F);

    // Aplicar método de Lax para actualizar las UP
    mac(P,UT, U, F, UP);

    // Aplicar condiciones de frontera a las UP
    boundary(UP);

    // Avanzar en el tiempo
    stepviscoso(U, UP);

    // Escribir a disco
    if (t >= tprint) {
      primitivas(U,P);
      output(U);
    }

  }

// Terminar
cout << "\nSe calcularon " << it << " iteraciones en "
     << (double)(clock() - start)/CLOCKS_PER_SEC << "s.\n\n";
}
