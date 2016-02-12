#include <iostream>
#include "include/armadillo"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

double initial_fcc_position(double b, double N_c, mat (&r))
{
    //Unit-cell-vectors
    vec ux(3); ux = {1,0,0}; ux *= b;
    vec uy(3); uy = {0,1,0}; uy *= b;
    vec uz(3); uz = {0,0,1}; uz *= b;
    r.zeros();
    vec origo(3);
    int number = -1;
    for (int ix = 0; ix < N_c; ix++)
    {
        for (int iy = 0; iy < N_c; iy++)
        {
            for (int iz = 0; iz < N_c; iz++)
            {
                origo = ix*ux + iy*uy + iz*uz;
                    //cout << origo(0) << setw(10) << origo(1) << setw(10) << origo(2) << endl;
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i); //lign. (3)
                }
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i) + 0.5*ux(i) + 0.5*uy(i); //lign. (4)
                }
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i) + 0.5*uy(i) + 0.5*uz(i); //lign. (5)
                }
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i) + 0.5*ux(i) + 0.5*uz(i); //lign. (6)
                }
            }
        }
    }
}

// generating velocity with gaussian distribution around 0 with standard deviation sqrt(k_B T / m)
void gaussian_velocity_generator(mat (&velocity), int number_of_particles, double Temperature, double mass)
{
  //standard_deviation = pow(k_B*Temperature / mass, 0.5);
  velocity.zeros();
  srand(time(NULL));
  double k_B = 10;
  double standard_deviation = pow(k_B*Temperature*mass,0.5);
  for (int i = 0; i < number_of_particles; i++)
  {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  for (int j = 0; j < 3; j++)
  {
    do{
      v1 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      v2 = 2.*((double) rand() / (RAND_MAX)) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    velocity(i,j) = standard_deviation*v2*fac;
  }
  }
} // end function for gaussian deviates

int main()
{
    double b = 5.260; //lattice constant for Argon in Å
    double N_c = 8; //Number of cells in x, y and z direction (cubic unit cell)
    double T = 100.0;
    double mass = 39.948;
    mat r(N_c*N_c*N_c*4,3);
    mat v(N_c*N_c*N_c*4,3);
    gaussian_velocity_generator(v,N_c*N_c*N_c*4,T,mass);
    initial_fcc_position(b,N_c,r);

    ofstream myfile ("a_initial_state.xyz");
            if (myfile.is_open())
            {
                myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
                myfile << "Position of argon atoms in fcc cell" << endl;
                for (int i = 0; i < N_c*N_c*N_c*4; i++)
                {
                    myfile << "Ar" << setw(10) << r(i,0) << setw(10) << r(i,1) << setw(10) << r(i,2) << setw(10) << v(i,0) << setw(10) << v(i,1) << setw(10) << v(i,2) << endl;
                }
            }

    return 0;
}

