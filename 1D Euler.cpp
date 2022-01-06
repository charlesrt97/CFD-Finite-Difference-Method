// Solves the 1-dimensional Euler equations (inviscid Navier-Stokes equations), using MacCormack's method (2nd order)
// for the case of a Sod shock-tube

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/******************************************************************************/
// constant paramters and global variables

// constant parameters used in the simulation
const int    NX = 400;         // mesh size
const double XL = 0.0;         // left physicial coordinate
const double XR = 1.0;         // right physical coordinate
const double TFIN = 0.2;       // time integration
const double CFL = 0.9;        // Courant number
const double dtprint = 0.01;   // time interval to write to disk

const double gam=1.4;

const int ieq=3;

const double eta = 0.1;

// for the schock-tube: coordinate of separation between initial states
const double X0 = 0.5;

// for advection: wave's velocity
const double A = 1.0;

const double DX = (XR-XL)/NX;      // Espaciamiento de la malla

// global variables
double U[ieq][NX+2];       // "current" conservative variables
double UP[ieq][NX+2];      // "advanced" conservative variables
double F[ieq][NX+2];       // physical fluxes
double P[ieq][NX+2];
double UT[ieq][NX+2];
double dt;            // time step
double t;          // current time
int it;               // current iteration
clock_t start;        // initial time
double tprint;        // time for the following output
int itprint;          // output number


/******************************************************************************/

// sets initial conditions
void initflow(double U[ieq][NX+2]) {

// initializes U in all the domain
// includes ghost cells
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

  // Initializes other variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;

}

/******************************************************************************/

// writes to disk the state of the simulation
void output(double U[ieq][NX+2]) {

  // generates output file's name
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // opens the file
  fstream fout(fname, ios::out);

  // writes U values to disk
  double x;
  for (int i=0; i <= NX; i++) {
    x = XL + i*DX;
    fout << x << " " << P[0][i] << " " << P[1][i] << " " << P[2][i] << " " << endl;
  }

  // closes the file
  fout.close();

  printf("Output %s\n", fname);

  itprint = itprint + 1;
  tprint = itprint * dtprint;

}

/******************************************************************************/

// applies boundary conditions to ghost cells
void boundary(double U[ieq][NX+2]) {

  for(int iieq=0; iieq<=2;iieq++){
    U[iieq][0]=U[iieq][1];
    U[iieq][NX+1]=U[iieq][NX];
  }
}

/******************************************************************************/

// computes primitives, including ghost cells
void primitivas(double U[ieq][NX+2], double P[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {
    //F[i]=A*U[i];
    P[0][i]=U[0][i];
    P[1][i]=U[1][i]/U[0][i];
    P[2][i]=(gam-1)*(U[2][i]-pow(U[1][i],2)/(2*U[0][i]));

  }

}


// computes physical fluxes, including ghost cells
void fluxes(double P[ieq][NX+2], double F[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {

    F[0][i]=P[0][i]*P[1][i];
    F[1][i]=P[0][i]*pow(P[1][i],2)+P[2][i];
    F[2][i]=P[1][i]*(0.5*P[0][i]*pow(P[1][i],2)+P[2][i]*gam/(gam-1));
  }

}

/******************************************************************************/

// computes new time step resulting from the CFL condition
double timestep(double P[ieq][NX+2]) {

  double dt;

  // for advection eq., max_u is simply A
  // double max_speed = abs(A);

  // for other cases, we shall compute the maximum value abs(vel)
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

// applies MacCormack's method to obtain the numerical intercell fluxes
void mac(double P[ieq][NX+2], double UT[ieq][NX+2], double U[ieq][NX+2], double F[ieq][NX+2], double UP[ieq][NX+2]) {

  for (int i=1; i<=NX; i++){
    for (int iieq=0; iieq<=2; iieq++){
      UT[iieq][i]=U[iieq][i]-dt/DX*(F[iieq][i+1]-F[iieq][i]);
      //UP[iieq][i]=(U[iieq][i+1]+U[iieq][i-1])/2.0-dt/(2*DX)*(F[iieq][i+1]-F[iieq][i-1]);
    }
  }

  boundary(UT);

  primitives(UT,P);

  fluxes(P,F);

  for (int i=1; i<=NX; i++){
    for (int iieq=0; iieq<=2; iieq++){
      UP[iieq][i]=(U[iieq][i]+UT[iieq][i])/2.0-dt/(2.0*DX)*(F[iieq][i]-F[iieq][i-1]);
      //UP[iieq][i]=(U[iieq][i+1]+U[iieq][i-1])/2.0-dt/(2*DX)*(F[iieq][i+1]-F[iieq][i-1]);
    }
  }

}

/******************************************************************************/

// this represents one complete time step
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

  // initial conditions and initializes variables
  initflow(U);

  primitives(U,P);

  // writes initial conditions to disk
  output(U);

  // simulation's initial time
  start = clock();
  while (t <= TFIN) {

    primitives(U,P); // updates primitives

    // updates time step
    dt = timestep(P);

    // updates phyisical fluxes
    fluxes(P, F);

    // applies MacCormack's method
    mac(P,UT, U, F, UP);

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
