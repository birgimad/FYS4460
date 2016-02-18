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
  //double k_B = 1.38 * pow(10,-23);
  double epsilon = 119.8; //without kB
  double standard_deviation = pow((Temperature*mass)/(epsilon),0.5);  //unitless ??
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

//Generate uniformly distributed velocities
void uniform_vel_generator(mat (&velocity), int Number_of_atoms, double std_avvik)
{
vec x(Number_of_atoms);
vec y(Number_of_atoms);
vec z(Number_of_atoms);
int minus_plus_numberX;
int minus_plus_numberY;
int minus_plus_numberZ;
srand(time(NULL));

for (int i=0;i<Number_of_atoms;i++){

        x(i) = ((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)
        y(i) = ((double) rand() / (RAND_MAX));
        z(i) = ((double) rand() / (RAND_MAX));
    }
for (int i=0;i<Number_of_atoms;i++){
    minus_plus_numberX = (1.99999)*((double) rand() / (RAND_MAX));
    minus_plus_numberY = (1.99999)*((double) rand() / (RAND_MAX));
    minus_plus_numberZ = (1.99999)*((double) rand() / (RAND_MAX));
    velocity(i,0)=x(i)*std_avvik*pow(-1,minus_plus_numberX);
    velocity(i,1)=y(i)*std_avvik*pow(-1,minus_plus_numberY);
    velocity(i,2)= z(i)*std_avvik*pow(-1,minus_plus_numberZ);
    }
}

double force_calculation(mat (&F), mat r, double N_c, double b, double sigma, int i)
//calculates total force F(i,k) on the i'th particle due to the presence of all other particles
{
    double dr, r_inverted, r_inverted_6, r_inverted_12, drx, dry, drz, cell_length, cell_length_half;
    cell_length = N_c * b;
    cell_length_half = cell_length / 2;
    vec r_vec(3), F_temp(3);
    F_temp.zeros();
    for (int j = i; j < N_c*N_c*N_c*4; j++)  //
    {
        //Consider only shortest distance bewteen two particles
        drx = r(i,0)-r(j,0);
        if (abs(drx) > cell_length_half)
        {
            drx -= cell_length;
        }
        dry = r(i,1)-r(j,1);
        if (abs(dry) > cell_length_half)
        {
            dry -= cell_length;
        }
        drz = r(i,2)-r(j,2);
        if (abs(drz) > cell_length_half)
        {
            drz -= cell_length;
        }

        dr = drx*drx + dry*dry + drz*drz;

        if (dr != 0 && dr < 3) //include critical distance. However, this is not the day to do it! Make lists instead!!
        {
        r_inverted = 1 / dr;
        r_inverted_6 = r_inverted * r_inverted * r_inverted;
        r_inverted_12 = r_inverted_6 * r_inverted_6;
        for (int k = 0; k < 3; k++)
        {
            r_vec(k) = (r(j,k) - r(i,k));
            F_temp(k) = 24 * (2 * r_inverted_12 - r_inverted_6) * r_inverted * r_vec(k);    //do not mult by 24 in the loop
            F(i,k) += F_temp(k);
            F(j,k) += F_temp(k);
        }
        }
    }
}

int main()
{
    double k_B = 1.38 * pow(10,-23);
    double sigma = 3.405;   //*pow(10,-10)
    double epsilon = 119.8;     //without k_B
    double m_Ar = 39.948;
    double F0 = epsilon/sigma;
    double t0 = sigma*pow(m_Ar/epsilon,0.5);
    double b = 5.260/sigma; //lattice constant for Argon in Å
    double N_c = 8; //Number of cells in x, y and z direction (cubic unit cell), 8
    double T = 100.0/epsilon;
    double mass = 39.948/m_Ar;
    double dt = pow(10,-7)/t0;
    cout << "dt = " << dt << endl;
    int number = 1;
    mat r(N_c*N_c*N_c*4,3);     //mult volume by 4, since there are 4 atoms in each fcc cell
    mat v(N_c*N_c*N_c*4,3);
    double std_avvik = 10;
    //gaussian_velocity_generator(v,N_c*N_c*N_c*4,T,mass);
    uniform_vel_generator(v,N_c*N_c*N_c*4,std_avvik);
    initial_fcc_position(b,N_c,r);

    int FileNumber = 0;
    int number_of_timesteps = 5000;







/*
    //saving initial state
    ofstream myfile ("DataFile_for_d_initial_state_with_periodic_bound.xyz");
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

    ofstream myfile ("DataFile_Velocities_initial_state.txt");
            if (myfile.is_open())
            {
                for (int i = 0; i < N_c*N_c*N_c*4; i++)
                {
                    myfile << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
                }
            }

    mat F(N_c*N_c*N_c*4,3);
    F.zeros();
    double unit_cell_size = N_c*b; //For periodic boundary condition

    for (int filecounter = 0; filecounter < number_of_timesteps; filecounter++) //for number_of_timesteps_integrations
    {

    for (int i = 0; i < N_c*N_c*N_c*4; i++)     //Loop for one integration over all particles
    {
        force_calculation(F,r,N_c,b,sigma,i);
        for (int k = 0; k < 3; k++) //looping over the three spacial coordinates
        {
            v(i,k) += (F(i,k) / (2*mass)) *dt;    //eq (7). Temp v_i
            r(i,k) += v(i,k) *dt;   //eq (8). New r_i

                //Periodic boundary conditions:
                //Check whether particle has gone through side of box every time the position is updated.
                if (r(i,k) < 0)
                {
                    r(i,k) = fmod(r(i,k),unit_cell_size) + unit_cell_size;     //r_i % unit_cell_size;
                }
                if (r(i,k) >= N_c*b )
                {
                    r(i,k) = fmod(r(i,k),unit_cell_size);     //r_i % unit_cell_size;
                }

            force_calculation(F,r,N_c,b,sigma,i); //eq (9). Calculate F_i(t+dt)
            v(i,k) += (F(i,k) / (2*mass)) *dt;
        }

    }
/*
    ofstream myfile;
        FileNumber++;
            string fileName = "DataFile_for_d_with_periodic_bound" + to_string(FileNumber) + ".xyz";

            cout << "File Name is " << fileName << endl;
            myfile.open(fileName, ios::app);
            myfile << N_c*N_c*N_c*4 << endl;    //number of atoms
            myfile << "Position of argon atoms in fcc cell after one integration" << endl;
            for (int i = 0; i < N_c*N_c*N_c*4; i++)
            {
                myfile << "Ar" << setw(20) << r(i,0) << setw(20) << r(i,1) << setw(20) << r(i,2) << setw(20) << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
            }
            myfile.close();
*/
    if ((number % 500) == 0)
    {
    ofstream myfile;
        FileNumber++;
            string fileName = "Datafile_velocities" + to_string(number) + ".txt";
            cout << "File Name is " << fileName << endl;
            myfile.open(fileName, ios::app);
            for (int i = 0; i < N_c*N_c*N_c*4; i++)
            {
                myfile << v(i,0) << setw(20) << v(i,1) << setw(20) << v(i,2) << endl;
            }
            myfile.close();
    }
    number += 1;
    }

    cout << "end" << endl;
    return 0;
}

