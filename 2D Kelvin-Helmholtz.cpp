// 2-dimensional Kelvin-Helmholtz instability, using MacCormack's method (2nd order)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

/******************************************************************************/

// constant parameters used for the simulation
const int    NX = 500;         // mesh size in the x-direction
const int    NY = 500;


const double X1 = 0.0;         // left physical coordinate
const double X2 = 1.0;         // right physical coordinate

const double Y1 = 0.0;
const double Y2 = 1.0;

const double TFIN = 10.0;       // integration time
//const double CFL = 0.5;        
const double CFL = 0.9;        // courant parameter
const double dtprint = 0.01;   // time interval to write to disk

const double gam=1.4;

//const double eta = 0.005;
//const double eta = 0.005;
const double eta = 0.005;

const int ieq=4;   // number of equations

// for the schock-tube: coordinate of separation between initial states
//const double R = 0.4;

const double A = 0.01;

//const double y0 = 0.5;
const double yy0 = 0.5;

const double DX = (X2-X1)/NX;      // Espaciamiento de la malla en x
const double DY = (Y2-Y1)/NY;

// global variables
double U[ieq][NX+2][NY+2];       // Variables conservadas actuales
double G[ieq][NX+2][NY+2];       // Variables conservadas actuales
double UP[ieq][NX+2][NY+2];      // Variables conservadas "avanzadas"
double GP[ieq][NX+2][NY+2];       // Variables conservadas "avanzadas"
double F[ieq][NX+2][NY+2];       // Flujos físicos
double P[ieq][NX+2][NY+2];
double UT[ieq][NX+2][NY+2];
//double c_x[NX+2];
//double c_y[NY+2];
double dt;            // time step
double t;          // current time
int it;               // current iteration
clock_t start;        // initial time
double tprint;        // time for the following output
int itprint;          // output number


/******************************************************************************/

// sets initial conditions
void initflow(double U[ieq][NX+2][NY+2]) {

// initializes U in all the domain
// includes ghost cellsfantasma
  
const double rho_arriba = 1.0;
const double u_arriba = -0.5;
const double v_arriba = 0.0;
const double p_arriba = 2.5;

const double rho_abajo = 2.0;
const double u_abajo = 0.5;
const double v_abajo = 0.0;
const double p_abajo = 2.5;

double x;
double y;
double rho;
double u;
double v;
double pres;

//srand(time(NULL));
srand(0);

for (int i = 0; i <= NX+1; i++){
  for (int j = 0; j <= NY+1; j++){
    x=X1+i*DX;
    y=Y1+j*DY;
    if (y <= yy0){
      rho = rho_abajo;
      u = u_abajo;
      v = v_abajo;
      pres = p_abajo;
    } else {
      rho = rho_arriba;
      u = u_arriba;
      v = v_arriba;
      pres = p_arriba;
    }
    u=u+A*(((double)rand()/RAND_MAX)-0.5);
    v=v+A*(((double)rand()/RAND_MAX)-0.5);

    U[0][i][j] = rho;
    U[1][i][j] = rho*u;
    U[2][i][j] = rho*v;
    U[3][i][j] = 0.5*rho*(u*u+v*v)+pres/(gam-1);

    }
  }
  // Initializes other variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;
}


  /* //--------------------------------------------------------
double x;
double y;
double r;

double ui2;
double vi2;
double uo2;
double vo2;

const double rhoi=2.0;
const double ui=0.5;
const double vi=0.0;
const double pii=2.5;

const double rhoo=1;
const double uo=-0.5;
const double vo=0.0;
const double po=2.5;

//srand(0);

srand(time(NULL)); //con time.h

  for (int i=0; i <= NX+1; i++){
    for (int j=0; j <= NY+1; j++){
      x=X1+i*DX;
      y=Y1+j*DY;
      ui2=ui+A*(((double)rand()/RAND_MAX)-0.5);
      vi2=vi+A*(((double)rand()/RAND_MAX)-0.5);
      uo2=uo+A*(((double)rand()/RAND_MAX)-0.5);
      vo2=vo+A*(((double)rand()/RAND_MAX)-0.5);
  //    printf("%f %f \n",ui2);
      if (y <= yy0){
        U[0][i][j]=rhoi;
        U[1][i][j]=rhoi*(ui2);
        U[2][i][j]=rhoi*(vi2);
        U[3][i][j]=0.5*rhoi*(ui2*ui2+vi2*vi2)+pii/(gam-1);
      } else {
        U[0][i][j]=rhoo;
        U[1][i][j]=rhoo*(uo2);
        U[2][i][j]=rhoo*(vo2);
        U[3][i][j]=0.5*rhoo*(uo2*uo2+vo2*vo2)+po/(gam-1);
      }
    }
  }

  // Inicializar otras variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;

*/ //--------------------------------------------------------






/******************************************************************************/

// writes to disk the state of the simulation
void output(double P[ieq][NX+2][NY+2]) {

  // generates output file's name
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // opens the file
  fstream fout(fname, ios::out);

  // writes U values to disk
  //double x;
for (int iieq=0; iieq <= ieq-1; iieq++){
    for (int i=1; i <= NX; i++){
      for (int j=1; j <= NY; j++){
        fout << P[iieq][i][j] << " ";
      }
      fout << endl;
    }
  }

  // closes the file
  fout.close();

  printf("Output %s\n", fname);

  itprint = itprint + 1;
  tprint = itprint * dtprint;

}

/******************************************************************************/

// applies boundary conditions to ghost cells
void boundary(double U[ieq][NX+2][NY+2]) {


    for (int i=1; i <= NX; i++){
      for (int iieq = 0; iieq <= ieq-1; iieq++){
      U[iieq][i][0]=U[iieq][i][1];     
      U[iieq][i][NY+1]=U[iieq][i][NY]; 
    }
  }


    for (int j=1; j <= NY; j++){
      for (int iieq = 0; iieq <= ieq-1; iieq++){
      U[iieq][0][j]=U[iieq][NX][j];     
      U[iieq][NX+1][j]=U[iieq][1][j];
    }
  }
}
/******************************************************************************/

// computes primitives, including ghost cells
void primitives(double U[ieq][NX+2][NY+2], double P[ieq][NX+2][NY+2]) {

  for (int i = 0; i <= NX+1; i++) {
    for (int j = 0; j <= NY+1; j++){
      P[0][i][j]=U[0][i][j];
      P[1][i][j]=U[1][i][j]/U[0][i][j];
      P[2][i][j]=U[2][i][j]/U[0][i][j];
      //P[3][i][j]=(gam-1)*(U[3][i][j]-pow(U[1][i][j],2)/(2*U[0][i][j])-pow(U[2][i][j],2)/(2*U[0][i][j]));
      P[3][i][j]=(gam-1)*(U[3][i][j]-0.5*(U[1][i][j]*U[1][i][j]+U[2][i][j]*U[2][i][j])/U[0][i][j]);
    }
  }
}

// computes physical fluxes, including ghost cells
void fluxes(double P[ieq][NX+2][NY+2], double F[ieq][NX+2][NY+2], double G[ieq][NX+2][NY+2]) {

  for (int i = 0; i <= NX+1; i++) {
    for (int j = 0; j <= NY+1; j++){

      F[0][i][j]=P[0][i][j]*P[1][i][j];
      //F[1][i][j]=P[0][i][j]*pow(P[1][i][j],2)+P[3][i][j];
      F[1][i][j]=P[0][i][j]*P[1][i][j]*P[1][i][j]+P[3][i][j];
      F[2][i][j]=P[0][i][j]*P[1][i][j]*P[2][i][j];
      //F[3][i][j]=P[1][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2))+P[3][i][j]*gam/(gam-1));
      F[3][i][j]=P[1][i][j]*(0.5*P[0][i][j]*(P[1][i][j]*P[1][i][j]+P[2][i][j]*P[2][i][j])+P[3][i][j]*gam/(gam-1));

      G[0][i][j]=P[0][i][j]*P[2][i][j];
      G[1][i][j]=P[0][i][j]*P[1][i][j]*P[2][i][j];
      //G[2][i][j]=P[0][i][j]*pow(P[2][i][j],2)+P[3][i][j];
      G[2][i][j]=P[0][i][j]*P[2][i][j]*P[2][i][j]+P[3][i][j];
      //G[3][i][j]=P[2][i][j]*(0.5*P[0][i][j]*(pow(P[1][i][j],2)+pow(P[2][i][j],2))+P[3][i][j]*gam/(gam-1));
      G[3][i][j]=P[2][i][j]*(0.5*P[0][i][j]*(P[1][i][j]*P[1][i][j]+P[2][i][j]*P[2][i][j])+P[3][i][j]*gam/(gam-1));
    }
  }
}


/******************************************************************************/

// computes new time step resulting from the CFL condition
double timestep(double P[ieq][NX+2][NY+2]) {

  double dt;
  double dtx;
  double dty;
  //double min;

  double max_speed_x = 0.0;
  double max_speed_y = 0.0;
  double cs;
  double kx;
  double ky;

  for (int i=1; i<=NX; i++){
    for (int j=1; j<=NY; j++){
      cs = sqrt(gam*P[3][i][j]/P[0][i][j]);
      kx = abs(P[1][i][j]) + cs;
      ky = abs(P[2][i][j]) + cs;
      if (kx > max_speed_x) max_speed_x = kx;
      if (ky > max_speed_y) max_speed_y = ky;
    }
  }

  dtx = CFL * DX / max_speed_x;

  //for (int i=1; i <= NX; i++){
  //  for (int j=1; j <= NY; j++){
  //    csy = sqrt(gam*P[3][i][j]/P[0][i][j]);
  //    ky = abs(P[2][i][j]) + csy;
  //    if (ky > max_speed_y) max_speed_y = ky;
  //  }
  //}

  dty = CFL * DY / max_speed_y;

  //if (dtx >= dty){
  //  min = dty;
  //} else {
  //  min = dtx;
  //}

  dt = min(dtx,dty)/ sqrt(2);

  printf("dt = %f, t = %f, it = %i \n",dt, t, it);

  return dt;

}

/******************************************************************************/

// applies the MacCormack's method to obtain the numerical intercell fluxes
void mac(double U[ieq][NX+2][NY+2], double F[ieq][NX+2][NY+2], double UP[ieq][NX+2][NY+2], double G[ieq][NX+2][NY+2],double UT[ieq][NX+2][NY+2],double P[ieq][NX+2][NY+2]) {

  //primitives(U,P);

  //fluxes(P,F,G);

for (int iieq=0; iieq<=ieq-1; iieq++){
  for (int i=1; i<=NX; i++){
    for (int j=1; j<=NY; j++){

      UT[iieq][i][j]=U[iieq][i][j]-dt/DX*(F[iieq][i+1][j]-F[iieq][i][j])-dt/DY*(G[iieq][i][j+1]-G[iieq][i][j]);

      }
    }
  }

  boundary(UT);

  primitives(UT,P);

  fluxes(P,F,G);


for (int iieq=0; iieq<=ieq-1; iieq++){
    for (int i=1; i<=NX; i++){
      for (int j=1; j <= NY; j++){

        UP[iieq][i][j]=(U[iieq][i][j]+UT[iieq][i][j])/2.0-dt/(2.0*DX)*(F[iieq][i][j]-F[iieq][i-1][j])-dt/(2.0*DY)*(G[iieq][i][j]-G[iieq][i][j-1]);

        //UP[iieq][i][j]=(U[iieq][i+1][j]+U[iieq][i-1][j]+U[iieq][i][j+1]+U[iieq][i][j-1])/4 - dt/(2*DX)*(F[iieq][i+1][j]-F[iieq][i-1][j]) - dt/(2*DY)*(G[iieq][i][j+1]-G[iieq][i][j-1]);

      }
    }
  }
}

/******************************************************************************/

// this represents one complete time step
void stepviscoso(double U[ieq][NX+2][NY+2], double UP[ieq][NX+2][NY+2]) {


  for (int iieq=0; iieq<=ieq-1; iieq++){
    for (int i = 0; i <= NX+1; i++) {
      for (int j = 0; j <= NY+1; j++){
        if ((i == 0) || (i == NX+1) || (j == 0) || (j == NY+1)){
          U[iieq][i][j]=UP[iieq][i][j];
        }
        else if ( ((U[iieq][i+1][j]-U[iieq][i][j])*(U[iieq][i][j]-U[iieq][i-1][j]) < 0) ||  ((U[iieq][i][j+1]-U[iieq][i][j])*(U[iieq][i][j]-U[iieq][i][j-1])<0) ) {
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

/******************************************************************************/

int main() {

  // initial conditions and initializes variables
  initflow(U);

  primitives(U,P);

  // writes initial conditions to disk
  output(P);

  // simulation's initial time
  start = clock();
  while (t <= TFIN) {

    primitives(U,P);

    // updates time step
    dt = timestep(P);

    // updates phyisical fluxes
    fluxes(P, F, G);

    // applies MacCormack's method
    mac(U, F, UP, G, UT, P);

    // applies boundary conditions to the UP-variables
    boundary(UP);
    
    // this represents one complete time step
    stepviscoso(U, UP);

    // writes to disk
    if (t >= tprint) {
      primitives(U,P);
      output(U);
    }

  }

// end
cout << "\n Number of iterations: " << it << ". Time: "
     << (double)(clock() - start)/CLOCKS_PER_SEC << "s.\n\n";
}
