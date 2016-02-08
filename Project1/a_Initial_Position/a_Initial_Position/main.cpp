#include <iostream>
#include "include/armadillo"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

int main()
{
    double b = 5.260; //lattice constant for Argon in Å
    double N_c = 8; //Number of cells in x, y and z direction (cubic unit cell)
    //Unit-cell-vectors
    vec ux(3); ux = {1,0,0}; ux *= b;
    vec uy(3); uy = {0,1,0}; uy *= b;
    vec uz(3); uz = {0,0,1}; uz *= b;
    mat r(N_c*N_c*N_c*4,3);
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
                    r(number,i) = origo(i) + (b/2)*ux(i) + (b/2)*uy(i); //lign. (4)
                }
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i) + (b/2)*uy(i) + (b/2)*uz(i); //lign. (5)
                }
                number += 1;
                for (int i = 0; i < 3; i++)
                {
                    r(number,i) = origo(i) + (b/2)*ux(i) + (b/2)*uz(i); //lign. (6)
                }
            }
        }
    }
    cout << r << endl;
    ofstream myfile ("a_initial_position.xyz");
            if (myfile.is_open())
            {
                myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
                myfile << "Position of argon atoms in fcc cell" << endl;
                for (int i = 0; i < N_c*N_c*N_c*4; i++)
                {
                    myfile << "Ar" << setw(10) << r(i,0) << setw(10) << r(i,1) << setw(10) << r(i,2) << endl;
                }
            }
    return 0;
}

