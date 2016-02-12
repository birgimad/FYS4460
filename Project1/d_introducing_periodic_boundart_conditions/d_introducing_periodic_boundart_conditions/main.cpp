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
  double k_B = 1.38 * pow(10,-23);
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
    double dt = 10;
    mat r(N_c*N_c*N_c*4,3);
    mat v(N_c*N_c*N_c*4,3);
    gaussian_velocity_generator(v,N_c*N_c*N_c*4,T,mass);
    initial_fcc_position(b,N_c,r);

    //saving initial state
    ofstream myfile ("DataFile_for_c3_initial_state.xyz");
            if (myfile.is_open())
            {
                myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
                myfile << "Position of argon atoms in fcc cell after one integration" << endl;
                for (int i = 0; i < N_c*N_c*N_c*4; i++)
                {
                    myfile << "Ar" << setw(20) << r(i,0) << setw(20) << r(i,1) << setw(20) << r(i,2) << setw(20) << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
                }
            }


    double dr, r_inverted, r_inverted_6, r_inverted_12;
    vec r_vec(3), F_temp(3);
    F_temp.zeros();
    mat F(N_c*N_c*N_c*4,3);  //
    F.zeros();

    int FileNumber = 0;
    int number_of_timesteps = 10;
    for (int filecounter = 0; filecounter < number_of_timesteps; filecounter++)
    {

    for (int i = 0; i < N_c*N_c*N_c*4; i++)     //Loop for one integration over all particles
    {
        for (int j = 0; j < N_c*N_c*N_c*4; j++)
        {
            dr = (r(i,0)-r(j,0))*(r(i,0)-r(j,0)) + (r(i,1)-r(j,1))*(r(i,1)-r(j,1)) + (r(i,2)-r(j,2))*(r(i,2)-r(j,2));
            if (dr != 0)
            {
            r_inverted = 1 / dr;
            r_inverted_6 = r_inverted * r_inverted * r_inverted;
            r_inverted_12 = r_inverted_6 * r_inverted_6;
            for (int k = 0; k < 3; k++)
            {
                r_vec(k) = (r(j,k) - r(i,k));
                F_temp(k) = 24 * (2 * r_inverted_12 - r_inverted_6) * r_inverted * r_vec(k); //do not mult by 24 in the loop
                F(i,k) += F_temp(k);
            }
            }
        }
        for (int k = 0; k < 3; k++)
        {
            v(i,k) += F(i,k) / (2*mass) *dt;
            r(i,k) += v(i,k) *dt;
            F(i,k) = 0; //setting U_i = 0, hnece do not include (9)
        }
    }

    //generating datafiles:
    ofstream myfile;
        FileNumber++;
            string fileName = "DataFile_for_c3" + to_string(FileNumber) + ".xyz";

            cout << "File Name is " << fileName << endl;
            myfile.open(fileName, ios::app);    //open file
            myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
            myfile << "Position of argon atoms in fcc cell after one integration" << endl;
            for (int i = 0; i < N_c*N_c*N_c*4; i++)
            {
                myfile << "Ar" << setw(20) << r(i,0) << setw(20) << r(i,1) << setw(20) << r(i,2) << setw(20) << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
            }
            myfile.close(); //close file

    }

    //cout << r << endl;

/*
    ofstream myfile ("c_state_after_one_integration_2.xyz");
            if (myfile.is_open())
            {
                myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
                myfile << "Position of argon atoms in fcc cell after one integration" << endl;
                for (int i = 0; i < N_c*N_c*N_c*4; i++)
                {
                    myfile << "Ar" << setw(20) << r(i,0) << setw(20) << r(i,1) << setw(20) << r(i,2) << setw(20) << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
                }
            }
*/

    return 0;
}


