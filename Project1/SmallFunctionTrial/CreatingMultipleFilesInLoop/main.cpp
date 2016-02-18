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
    int number_of_atoms = 5;
    mat r(number_of_atoms,3);       //3-D
    r.zeros();
    r(0,0) = 0.7; r(0,1) = 0.5; r(0,2) = 2;
    r(1,0) = 2.1; r(1,1) = 1.2; r(1,2) = 0;
    r(2,0) = 2.1; r(2,1) = 2.2; r(2,2) = 0;
    r(3,0) = 2.4; r(3,1) = 2.7; r(3,2) = 0;
    r(4,0) = 1.2; r(4,1) = 3.1; r(4,2) = 0;
    double length_of_unit_cell = 4;
    double length_of_rc_cell = 1;
    int number_of_cells = length_of_unit_cell / length_of_rc_cell;

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
    //putting atoms in correct cells
    int cell_number_x, cell_number_y, cell_number_z;
    for (int i = 0; i < number_of_atoms; i++)
    {
    cell_number_x = (r(i,0) / length_of_unit_cell) *number_of_cells;
    cell_number_y = (r(i,1) / length_of_unit_cell) *number_of_cells;
    cell_number_z = (r(i,2) / length_of_unit_cell) *number_of_cells;
    fourthList[cell_number_x][cell_number_y][cell_number_z].push_back(i);
    }

int Nabo_x_value, Nabo_y_value, Nabo_z_value, Nabo_number;
double dr_x, dr_y, dr_z, dr;

    for (int box_x = 0; box_x < number_of_cells; box_x++)
    {
    for (int box_y = 0; box_y < number_of_cells; box_y++)
    {
    for (int box_z = 0; box_z < 1; box_z++) //number_of_cells
    {
        for (int atom_in_box = 0; atom_in_box < fourthList[box_x][box_y][box_z].size(); atom_in_box++)
        {

            cout << "Naboer til: " << fourthList[box_x][box_y][box_z][atom_in_box] << endl;

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
                if (fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box] != fourthList[box_x][box_y][box_z][atom_in_box])
                {
                dr_x = (r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],0)-r(fourthList[box_x][box_y][box_z][atom_in_box],0));
                dr_y = (r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],1)-r(fourthList[box_x][box_y][box_z][atom_in_box],1));
                dr_z = (r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],2)-r(fourthList[box_x][box_y][box_z][atom_in_box],2));

                dr = dr_x*dr_x +dr_y*dr_y + dr_z*dr_z;
                cout << "(" << r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],0) << "," << r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],1) << "," << r(fourthList[Nabo_x_value][Nabo_y_value][Nabo_z_value][atom_in_nabo_box],2) << ")" << setw(10) << dr << endl;
                }
            }
        }
        }
        }
        }

    }
    }
    }


    return 0;
}


/*vec r(5);
    r(0) = 1; r(1) = 1.5; r(2) = 2; r(3) = 7; r(4) = 0;
    double r_c = 2;
    r /= r_c;
    double length_of_unit_cell = 10;
    double number_of_boxes = length_of_unit_cell / r_c;
    int t1 = 0;
    int number_of_particles = 5;
    for (int i = 0; i < number_of_particles; i++)
    {
    t1 = 0;
    for (int t = 0; t < number_of_boxes; t++)
    {
        if (t1 <= r(i) && r(i) < t1+1 )
        {
            cout << t1 << setw(10) << r(i) << setw(10) << i << endl;
        }
        t1 += 1;
    }
    }
*/

/*
    typedef struct node
    {
    int data;
    node *head = NULL;             //empty linked list
    node *temp;             //create a temporary node
    temp = (node*)malloc(sizeof(node)); //allocate space for node
    temp->data = info;             // store data(first field)
    temp->next=head;  // store the address of the pointer head(second field)
    head = temp;                  // transfer the address of 'temp' to 'head'
    }
*/

double periodic_bound_con()
{
    vec r(5);
    r.zeros();
    r(0) = 0;
    r(1) = 2.5;
    double dr = r(1) -  r(0);
    double cell_length = 5;
    double cell_length_half = cell_length/2;
    if (dr > cell_length_half)
    {
        dr -= cell_length;
    }
    cout << dr << endl;
}

double uniform_velocity()
{
    mat velocity(10,3);
    double std_avvik = 3;
    int Number_of_atoms = 10;
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

//Periodic boundary conditions:
//Check whether particle has gone through side of box every time the position is updated.
double periodic_boundary()
{
    double N_c = 3;
    double b = 1;
    vec r(3);
    r(0) = -7.1;
    r(1) = 1.7;
    r(2) = 3;

    int unit_cell_size = N_c*b;
    for (int i = 0; i < 3; i++) //looping over the three spacial coordinates
    {
        if (r(i) < 0)
        {
            r(i) = fmod(r(i),unit_cell_size) + unit_cell_size;     //r_i % unit_cell_size;
        }
        if (r(i) >= N_c*b )
        {
            r(i) = fmod(r(i),unit_cell_size);     //r_i % unit_cell_size;
        }
    }
}

//Generating multiple files in loop with name: afile1.txt, afile2.txt etc.
double Generating_multiple_files()
{
    int FileNumber = 0;
    int number_of_timesteps = 3;
    for (int filecounter = 0; filecounter < number_of_timesteps; filecounter++)
    {
    ofstream myfile;
        FileNumber++;
            string fileName = "afile" + to_string(FileNumber) + ".txt";

            cout << "File Name is " << fileName << endl;
            myfile.open(fileName, ios::app);
            myfile << "Blah blah blah! File sequence is " << FileNumber << "\n";
            myfile.close();
    }
}

/*
//generating Nabo-list
double Genarating_nabo_list()
{
    int number_of_atoms = 5;
    mat r(number_of_atoms,2);       //2-D
    r.zeros();
    r(0,0) = 0.7; r(0,1) = 0.5; //r(0,2) = 0.1;
    r(1,0) = 2.1; r(1,1) = 1.2; //r(1,2) = 0.3;
    r(2,0) = 2.1; r(2,1) = 2.2; //r(2,2) = 0.3;
    r(3,0) = 2.4; r(3,1) = 2.7; //r(3,2) = 0.3;
    r(4,0) = 1.2; r(4,1) = 3.1;

    r(0,0) = 0.7; r(0,1) = 0.5; r(0,2) = 0.1;
    r(1,0) = 0.1; r(1,1) = 0.2; r(1,2) = 0.3;
    r(2,0) = 0.1; r(2,1) = 0.2; r(2,2) = 0.3;
    r(3,0) = 0.1; r(3,1) = 0.2; r(3,2) = 0.3;

    //r(2,0) = 2.1; r(2,1) = 0.2; r(2,2) = 1.3;
    //r(3,0) = 2.1; r(3,1) = 0.2; r(3,2) = 3.3;
    double length_of_unit_cell = 4;
    double length_of_rc_cell_x = 1;
    double length_of_rc_cell_y = 1;
    double length_of_rc_cell_z = 1;
    int number_of_cells_x = length_of_unit_cell / length_of_rc_cell_x;
    int number_of_cells_y = length_of_unit_cell / length_of_rc_cell_y;
    int number_of_cells_z = length_of_unit_cell / length_of_rc_cell_z;

    int cell_number_x, cell_number_y, cell_number_z;

    vector<int> firstList;
    vector< vector <int> > secondList;
    vector< vector< vector <int> > > thirdList;
    vector< vector< vector< vector <int> > > > fourthList;
    for (int t=0; t<number_of_cells_x; t++)
    {
        secondList.push_back(firstList);
    }
    for (int t=0; t<number_of_cells_x;t++)
    {
        thirdList.push_back(secondList);
    }

    //fourthList.push_back(thirdList);
    //fourthList.push_back(thirdList);
    //fourthList.push_back(thirdList);
    //fourthList.push_back(thirdList);

    for (int i = 0; i < 5; i++)
    {
    cell_number_x = (r(i,0) / length_of_unit_cell) *number_of_cells_x;
    cell_number_y = (r(i,1) / length_of_unit_cell) *number_of_cells_y;
    //cell_number_z = (r(i,2) / length_of_unit_cell) *number_of_cells_z;

    //cout << r(i,0) << setw(10) << cell_number_x << endl;
    //cout << r(i,1) << setw(10) << cell_number_y << endl;
    //cout << r(i,2) << setw(10) << cell_number_z << endl;
    thirdList[cell_number_x][cell_number_y].push_back(i);  //[cell_number_z]
    }

    cout << fourthList[0][0][0][0] << endl;
    cout << fourthList[0][0][0][1] << endl;
    cout << fourthList[2][0][1][0] << endl;
    cout << fourthList[2][0][3][0] << endl;


    for (int k = 0; k <  fourthList[0][0][0].size(); k++)
    {
        cout << fourthList[0][0][0][k] << endl;
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "elements in (" << i << "," << j << ")" << endl;
            for (int k = 0; k <  thirdList[i][j].size(); k++)
            {
                cout << thirdList[i][j][k] << endl;
            }
        }
    }

imat Nabo_boxes(9,2);
int Nabo_x_value, Nabo_y_value, Nabo_number;
double dr_x, dr_y, dr;

    for (int box_x = 0; box_x < number_of_cells_x; box_x++)
    {
    for (int box_y = 0; box_y < number_of_cells_y; box_y++)
    {
        for (int atom_in_box = 0; atom_in_box <  thirdList[box_x][box_y].size(); atom_in_box++)
        {
            cout << thirdList[box_x][box_y][atom_in_box] << endl;
            Nabo_number = 0;
        for (int index_nabo_x = -1; index_nabo_x < 2; index_nabo_x++)
        {
            Nabo_x_value = (number_of_cells_x + box_x + index_nabo_x) % number_of_cells_x;
        for (int index_nabo_y = -1; index_nabo_y < 2; index_nabo_y++)
        {
            Nabo_y_value = (number_of_cells_y +box_y + index_nabo_y) % number_of_cells_y;
            Nabo_boxes(Nabo_number,0) = Nabo_x_value;
            Nabo_boxes(Nabo_number,1) = Nabo_y_value;
            Nabo_number += 1;
            for (int atom_in_nabo_box = 0; atom_in_nabo_box < thirdList[Nabo_x_value][Nabo_y_value].size(); atom_in_nabo_box++)
            {
                if (thirdList[Nabo_x_value][Nabo_y_value][atom_in_nabo_box] != thirdList[box_x][box_y][atom_in_box])
                {
                dr_x = (r(thirdList[Nabo_x_value][Nabo_y_value][atom_in_nabo_box],0)-r(thirdList[box_x][box_y][atom_in_box],0));
                dr_y = (r(thirdList[Nabo_x_value][Nabo_y_value][atom_in_nabo_box],1)-r(thirdList[box_x][box_y][atom_in_box],1));
                dr = dr_x*dr_x +dr_y*dr_y;
                cout << "(" << r(thirdList[Nabo_x_value][Nabo_y_value][atom_in_nabo_box],0) << "," << r(thirdList[Nabo_x_value][Nabo_y_value][atom_in_nabo_box],1) << ")" << setw(10) << dr << endl;
                }
            }
        }
        }
            //cout << Nabo_boxes << endl;
        }

    }
    }
}
*/
