#include <iostream>
#include "include/armadillo"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

//generating initial position
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

double force_calculation(mat (&F), mat r, double N_c, double b, double sigma, int i, int j)
//calculates force between particle i and j (fra j på i)
{
    double dr, r_inverted, r_inverted_6, r_inverted_12, drx, dry, drz, cell_length, cell_length_half;
    cell_length = N_c * b;
    cell_length_half = cell_length / 2;
    vec r_vec(3), F_temp(3);
    F_temp.zeros();

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
            //F(j,k) += F_temp(k);
        }
        }
}

int main()
{
    double k_B = 1.38 * pow(10,-23), sigma = 3.405, epsilon = 119.8, m_Ar = 39.948;
    double F0 = epsilon/sigma, t0 = sigma*pow(m_Ar/epsilon,0.5), b = 5.260/sigma, T = 100.0/epsilon, mass = 39.948/m_Ar;
    double N_c = 8; //Number of cells in x, y and z direction (cubic unit cell), 8
    double dt = pow(10,-7)/t0;  //time step length
    int number_of_atoms = 4*N_c*N_c*N_c;
    mat r(number_of_atoms,3);       //3-D
    r.zeros();
    mat F(number_of_atoms,3);
    F.zeros();
    mat v(number_of_atoms,3);
    v.zeros();
    double std_avvik = 10;
    //gaussian_velocity_generator(v,N_c*N_c*N_c*4,T,mass);
    uniform_vel_generator(v,N_c*N_c*N_c*4,std_avvik);
    initial_fcc_position(b,N_c,r);

    double unit_cell_size = N_c*b; //For periodic boundary condition
    double length_of_rc_cell = 2*b; //approx 3*sigma
    double length_of_unit_cell = unit_cell_size;
    int number_of_cells = length_of_unit_cell / length_of_rc_cell;
    int number_of_timesteps = 10;

    int i_temp, j_temp;
    int Nabo_x_value, Nabo_y_value, Nabo_z_value;
//---||---
    vector<int> firstList;
    vector< vector <int> > secondList;
    vector< vector< vector <int> > > thirdList;
    vector< vector< vector< vector <int> > > > fourthList;
    for (int t=0; t<number_of_cells; t++)
    {
        secondList.push_back(firstList);
    }
    for (int t=0; t<number_of_cells;t++)
    {
        thirdList.push_back(secondList);
    }
    for (int t=0; t< number_of_cells;t++)
    {
        fourthList.push_back(thirdList);
    }
//---||---

    int cell_number_x, cell_number_y, cell_number_z;

//START!
    for (int timestep = 0; timestep < number_of_timesteps; timestep++)
    {
        for (int x = 0; x < number_of_cells; x++)
        {
            for (int y = 0; y < number_of_cells; y++)
            {
                for (int z = 0; z < number_of_cells; z++)
                {
                    fourthList[x][y][z].clear();
                }
            }
        }
    //putting atoms in correct cells
    for (int i = 0; i < number_of_atoms; i++)
    {
    cell_number_x = (r(i,0) / length_of_unit_cell) *number_of_cells;
    cell_number_y = (r(i,1) / length_of_unit_cell) *number_of_cells;
    cell_number_z = (r(i,2) / length_of_unit_cell) *number_of_cells;
    fourthList[cell_number_x][cell_number_y][cell_number_z].push_back(i);
    }

    for (int box_x = 0; box_x < number_of_cells; box_x++)
    {
    for (int box_y = 0; box_y < number_of_cells; box_y++)
    {
    for (int box_z = 0; box_z < number_of_cells; box_z++)
    {
        for (int atom_in_box = 0; atom_in_box < fourthList[box_x][box_y][box_z].size(); atom_in_box++)
        {
            i_temp = fourthList[box_x][box_y][box_z][atom_in_box];

        for (int index_nabo_x = -1; index_nabo_x < 2; index_nabo_x++)
        {
            Nabo_x_value = (number_of_cells + box_x + index_nabo_x) % number_of_cells;
        for (int index_nabo_y = -1; index_nabo_y < 2; index_nabo_y++)
        {
            Nabo_y_value = (number_of_cells +box_y + index_nabo_y) % number_of_cells;
        for (int index_nabo_z = -1; index_nabo_z < 2; index_nabo_z++)
        {
            Nabo_z_value = (number_of_cells +box_z + index_nabo_z) % number_of_cells;

            for (int atom_in_nabo_box = 0; atom_in_nabo_box < fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value].size(); atom_in_nabo_box++)
            {
                j_temp = fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box];
                if (j_temp != i_temp)
                {
                force_calculation(F,r,N_c,b,sigma,i_temp,j_temp);
                }
            }
        }
        }
        }
            for (int k = 0; k < 3; k++) //looping over the three spacial coordinates
            {
                v(i_temp,k) += (F(i_temp,k) / (2*mass)) *dt;    //eq (7). Temp v_i
                r(i_temp,k) += v(i_temp,k) *dt;   //eq (8). New r_i

                    //Periodic boundary conditions:
                    //Check whether particle has gone through side of box every time the position is updated.
                    if (r(i_temp,k) < 0)
                    {
                        r(i_temp,k) = fmod(r(i_temp,k),unit_cell_size) + unit_cell_size;     //r_i % unit_cell_size;
                    }
                    if (r(i_temp,k) >= N_c*b )
                    {
                        r(i_temp,k) = fmod(r(i_temp,k),unit_cell_size);     //r_i % unit_cell_size;
                    }
            }
            for (int index_nabo_x = -1; index_nabo_x < 2; index_nabo_x++)
            {
                Nabo_x_value = (number_of_cells + box_x + index_nabo_x) % number_of_cells;
            for (int index_nabo_y = -1; index_nabo_y < 2; index_nabo_y++)
            {
                Nabo_y_value = (number_of_cells +box_y + index_nabo_y) % number_of_cells;
            for (int index_nabo_z = -1; index_nabo_z < 2; index_nabo_z++)
            {
                Nabo_z_value = (number_of_cells +box_z + index_nabo_z) % number_of_cells;
            for (int atom_in_nabo_box = 0; atom_in_nabo_box < fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value].size(); atom_in_nabo_box++)
            {
                j_temp = fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box];
                if (j_temp != i_temp)
                {
                force_calculation(F,r,N_c,b,sigma,i_temp,j_temp);
                }
            }
            }
            }
            }
            for (int k = 0; k < 3; k++)
            {
                v(i_temp,k) += (F(i_temp,k) / (2*mass)) *dt;
            }
        }
    }
    }
    }
    cout << timestep << endl;
    } //end loop

    return 0;
}

