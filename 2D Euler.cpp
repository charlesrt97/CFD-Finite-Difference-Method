// Solves the 2-dimensional Euler equations (inviscid Navier-Stokes equations), using MacCormack's method (2nd order)
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
const int    NX = 500;         // Tamaño de la malla en x
const int    NY = 500;


const double X1 = 0.0;         // Coordenada física del extremo izquierdo
const double X2 = 2.0;         // Coordenada física del extremo derecho

const double Y1 = 0.0;
const double Y2 = 2.0;

const double TFIN = 1.0;       // Tiempo final de integración
const double CFL = 0.2;        // Parametro de Courant
const double dtprint = 0.1;   // Intervalo para escribir a disco

const double gam=1.4;

const double eta = 0.1;

const int ieq=4;   // numero de ecuaciones

// Para tubo de choque: coord de la separación entre estados iniciales
const double R = 0.4;

// Para ecuación de advección: velocidad de ondas ('a' en las notas)
const double A = 1.0;

// Constantes derivadas de las anteriores
const double DX = (X2-X1)/NX;      // Espaciamiento de la malla en x
const double DY = (Y2-Y1)/NY;

// Variables globales
double U[ieq][NX+2][NY+2];       // Variables conservadas actuales
double G[ieq][NX+2][NY+2];       // Variables conservadas actuales
double UP[ieq][NX+2][NY+2];      // Variables conservadas "avanzadas"
double GP[ieq][NX+2][NY+2];       // Variables conservadas "avanzadas"
double F[ieq][NX+2][NY+2];       // Flujos físicos
double P[ieq][NX+2][NY+2];
double UT[ieq][NX+2][NY+2];
double c_x[NX+2];
double c_y[NY+2];
double dt;            // Paso de tiempo
double t;          // Tiempo actual
int it;               // Iteración actual
clock_t start;        // Tiempo de inicio
double tprint;        // Tiempo para el siguiente output
int itprint;          // Número de salida


/******************************************************************************/

// Impone las condiciones iniciales
void initflow(double U[ieq][NX+2][NY+2]) {

  // Inicializar los valores de U en todo el dominio
  // Nótese que también llenamos las celdas fantasma
double x;
double y;
double r;
const double xc=1;
const double yc=1;

const double rhoi=1.0;
const double ui=0.0;
const double vi=0.0;
const double pii=1.0;

const double rhoo=0.125;
const double uo=0.0;
const double vo=0.0;
const double po=0.1;

  for (int i=0; i <= NX+1; i++){
    for (int j=0; j <= NY+1; j++){
      x=X1+i*DX;
      y=Y1+j*DY;
      r = sqrt(pow(x-xc,2)+pow(y-yc,2));
      if (r <= R){
        U[0][i][j]=rhoi;
        U[1][i][j]=rhoi*ui;
        U[2][i][j]=rhoi*vi;
        U[3][i][j]=0.5*rhoi*(ui*ui+vi*vi)+pii/(gam-1);
      } else {
        U[0][i][j]=rhoo;
        U[1][i][j]=rhoo*uo;
        U[2][i][j]=rhoo*vo;
        U[3][i][j]=0.5*rhoo*(uo*uo+vo*vo)+po/(gam-1);
      }
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
void output(double U[ieq][NX+2][NY+2]) {

  // Generar el nombre del archivo de salida
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // Abrir el archivo
  fstream fout(fname, ios::out);

  // Escribir los valores de U al archivo
  double x;
for (int iieq=0; iieq <= ieq-1; iieq++){
    for (int i=1; i <= NX; i++){
      for (int j=1; j <= NY; j++){
        fout << P[iieq][i][j] << " ";
      }
      fout << endl;
    }
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
void boundary(double U[ieq][NX+2][NY+2]) {

  for (int iieq = 0; iieq <= ieq-1; iieq++){
    for (int i=0; i <= NX+1; i++){

      U[iieq][i][NY+1]=U[iieq][i][NY]; // lado arriba
      U[iieq][i][0]=U[iieq][i][1];     // lado abajo
    }
  }

  for (int iieq = 0; iieq <= ieq-1; iieq++){
    for (int j=0; j <= NY+1; j++){
      U[iieq][NX+1][j]=U[iieq][NX][j]; // lado derecho
      U[iieq][0][j]=U[iieq][1][j];     // lado izquierdo
    }
  }
}
/******************************************************************************/

// Calcular las primitivas
void primitivas(double U[ieq][NX+2][NY+2], double P[ieq][NX+2][NY+2]) {

  for (int i=0; i<=NX+1; i++) {
    for (int j=0; j<=NY+1; j++){
      P[0][i][j]=U[0][i][j];
      P[1][i][j]=U[1][i][j]/U[0][i][j];
      P[2][i][j]=U[2][i][j]/U[0][i][j];
      P[3][i][j]=(gam-1)*(U[3][i][j]-pow(U[1][i][j],2)/(2*U[0][i][j])-pow(U[2][i][j],2)/(2*U[0][i][j]));
    }
  }
}

// Calcular los flujos físicos F !tambien se incluyen las celdas fantasma
void fluxes(double P[ieq][NX+2][NY+2], double F[ieq][NX+2][NY+2], double G[ieq][NX+2][NY+2]) {

  for (int i=0; i<=NX+1;i++) {
    for (int j=0; j<=NY+1;j++){

      F[0][i][j]=P[0][i][j]*P[1][i][j];
      F[1][i][j]=P[0][i][j]*pow(P[1][i][j],2)+P[3][i][j];
      F[2][i][j]=P[0][i][j]*P[1][i][j]*P[2][i][j];
      F[3][i][j]=P[1][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2))+P[3][i][j]*gam/(gam-1));

      G[0][i][j]=P[0][i][j]*P[2][i][j];
      G[1][i][j]=P[0][i][j]*P[1][i][j]*P[2][i][j];
      G[2][i][j]=P[0][i][j]*pow(P[2][i][j],2)+P[3][i][j];
      G[3][i][j]=P[2][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2))+P[3][i][j]*gam/(gam-1));

    }
  }
}


/******************************************************************************/

// Calcula el paso de tiempo resultante de la condición CFL
double timestep(double P[ieq][NX+2][NY+2]) {

  double dt;
  double dtx;
  double dty;
  double min;

  double max_speed_x = 0.0;
  double max_speed_y = 0.0;
  double csx;
  double csy;
  double kx;
  double ky;

  for (int i=1; i<=NX; i++){
    for (int j=1; j<=NY; j++){
      csx = sqrt(gam*P[3][i][j]/P[0][i][j]);
      kx = abs(P[1][i][j]) + csx;
      if (kx > max_speed_x) max_speed_x = kx;
    }
  }

  dtx = DX / max_speed_x;

  for (int i=1; i <= NX; i++){
    for (int j=1; j <= NY; j++){
      csy = sqrt(gam*P[3][i][j]/P[0][i][j]);
      ky = abs(P[2][i][j]) + csy;
      if (ky > max_speed_y) max_speed_y = ky;
    }
  }

  dty = DY / max_speed_y;

  if (dtx > dty){
    min = dty;
  } else {
    min = dtx;
  }

  dt = CFL * min / sqrt(2);

  return dt;

}

/******************************************************************************/

// Aplica el método de Lax para obtener las UP a partir de las U
// Supone que los flujos F ya fueron actualizados
void mac(double U[ieq][NX+2][NY+2], double F[ieq][NX+2][NY+2], double UP[ieq][NX+2][NY+2], double G[ieq][NX+2][NY+2],double UT[ieq][NX+2][NY+2],double P[ieq][NX+2][NY+2]) {

for (int iieq=0; iieq<=ieq-1; iieq++){
  for (int i=0; i<=NX+1; i++){
    for (int j=0; j<=NY+1; j++){

      UT[iieq][i][j]=U[iieq][i][j]-dt/DX*(F[iieq][i+1][j]-F[iieq][i][j])-dt/DY*(G[iieq][i][j+1]-G[iieq][i][j]);

      }
    }
  }

  boundary(UT);

  primitivas(UT,P);

  fluxes(P,F,G);

  //fluxes(P,G);


  for (int iieq=0; iieq<=ieq-1; iieq++){
    for (int i=0; i<=NX+1; i++){
      for (int j=0; j <= NY+1; j++){

        UP[iieq][i][j]=(U[iieq][i][j]+UT[iieq][i][j])/2.0-dt/(2.0*DX)*(F[iieq][i][j]-F[iieq][i-1][j])-dt/(2.0*DY)*(G[iieq][i][j]-G[iieq][i][j-1]);

        //UP[iieq][i][j]=(U[iieq][i+1][j]+U[iieq][i-1][j]+U[iieq][i][j+1]+U[iieq][i][j-1])/4 - dt/(2*DX)*(F[iieq][i+1][j]-F[iieq][i-1][j]) - dt/(2*DY)*(G[iieq][i][j+1]-G[iieq][i][j-1]);

      }
    }
  }
}

/******************************************************************************/

// Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
void stepviscoso(double U[ieq][NX+2][NY+2], double UP[ieq][NX+2][NY+2]) {

for (int iieq=0; iieq<=ieq-1; iieq++){
  for (int i = 1; i <= NX; i++) {
    for (int j=1; j<= NY;j++){
      //U[iieq][i]=UP[iieq][i];
      if ( ((U[iieq][i+1][j]-U[iieq][i][j])*(U[iieq][i][j]-U[iieq][i-1][j])<0) ||  ((U[iieq][i][j+1]-U[iieq][i][j])*(U[iieq][i][j]-U[iieq][i][j-1])<0)) {
        U[iieq][i][j]=UP[iieq][i][j]+eta*(UP[iieq][i+1][j]+UP[iieq][i-1][j]+UP[iieq][i][j+1]+UP[iieq][i][j-1]-4*UP[iieq][i][j]);
      }
      else {
        U[iieq][i][j]=UP[iieq][i][j];
      }

    }
  }
  }

  t = t + dt;
  it = it + 1;

}


//  for (int iieq=0; iieq<=ieq-1; iieq++){
//    for (int i = 1; i <= NX; i++) {
  //    for (int j = 1; j <= NY; j++){
  //      U[iieq][i][j]=UP[iieq][i][j];
  //    }
  //  }
  //}

  //t = t + dt;
  //it = it + 1;

//}

/******************************************************************************/

int main() {

  // Condición inicial e inicializaciones
  initflow(U);

  primitivas(U,P);

  // Escribir condición inicial a disco
  output(P);

  // Tiempo de inicio de la simulación
  start = clock();
  while (t <= TFIN) {

    primitivas(U,P);

    // Actualizar el paso de tiempo
    dt = timestep(U);

    // Actualizar flujos físicos
    fluxes(P, F, G);

    // Aplicar método de Lax para actualizar las UP
    mac(U, F, UP, G, UT,P);

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
