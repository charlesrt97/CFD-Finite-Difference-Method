// Solves the 2-dimensional Euler equations (inviscid Navier-Stokes equations), using MacCormack's method (2nd order)
// for the case of a cylindrical explosion

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/******************************************************************************/

// constant parameters used for the simulation
const int    NX = 500;         // mesh size in the x-direction
const int    NY = 500;


const double X1 = 0.0;         // left physical coordinate
const double X2 = 2.0;         // right physical coordinate

const double Y1 = 0.0;
const double Y2 = 2.0;

const double TFIN = 1.0;       // integration time
const double CFL = 0.2;        // Courant parameter
const double dtprint = 0.1;   // time interval to write to disk

const double gam=1.4;

const double eta = 0.1;

const int ieq=4;   // number of equations

// for the schock-tube: coordinate of separation between initial states
const double R = 0.4;

// for advection: wave's velocity
const double A = 1.0;

const double DX = (X2-X1)/NX;      // Espaciamiento de la malla en x
const double DY = (Y2-Y1)/NY;

// global variables
double U[ieq][NX+2][NY+2];       // "current" conservative variables
double G[ieq][NX+2][NY+2];       // "current" conservative variables
double UP[ieq][NX+2][NY+2];      // "advanced" conservative variables
double GP[ieq][NX+2][NY+2];       // "advanced" conservative variables
double F[ieq][NX+2][NY+2];       // phyisical fluxes
double P[ieq][NX+2][NY+2];
double UT[ieq][NX+2][NY+2];
double c_x[NX+2];
double c_y[NY+2];
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
// includes ghost cells
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

  // Initializes other variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;

}

/******************************************************************************/

// writes to disk the state of the simulation
void output(double U[ieq][NX+2][NY+2]) {

  // generates output file's name
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // opens the file
  fstream fout(fname, ios::out);

  // writes U values to disk
  double x;
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

  for (int iieq = 0; iieq <= ieq-1; iieq++){
    for (int i=0; i <= NX+1; i++){

      U[iieq][i][NY+1]=U[iieq][i][NY]; 
      U[iieq][i][0]=U[iieq][i][1];    
    }
  }

  for (int iieq = 0; iieq <= ieq-1; iieq++){
    for (int j=0; j <= NY+1; j++){
      U[iieq][NX+1][j]=U[iieq][NX][j];
      U[iieq][0][j]=U[iieq][1][j];   
    }
  }
}
/******************************************************************************/

// computes primitives, including ghost cells
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

// computes physical fluxes, including ghost cells
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

// computes new time step resulting from the CFL condition
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

// applies the MacCormack's method to obtain the numerical intercell fluxes
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

// this represents one complete time step
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

  // initial conditions and initializes variables
  initflow(U);

  primitivas(U,P);

  // writes initial conditions to disk
  output(P);

  // simulation's initial time
  start = clock();
  while (t <= TFIN) {

    primitivas(U,P);

    // updates time step
    dt = timestep(U);

    // updates phyisical fluxes
    fluxes(P, F, G);

    // applies MacCormack's method
    mac(U, F, UP, G, UT,P);

    // applies boundary conditions to the UP-variables
    boundary(UP);

    // this represents one complete time step
    stepviscoso(U, UP);

    // writes to disk
    if (t >= tprint) {
      primitivas(U,P);
      output(U);
    }

  }

// end
cout << "\n Number of iterations: " << it << ". Time: "
     << (double)(clock() - start)/CLOCKS_PER_SEC << "s.\n\n";
}
