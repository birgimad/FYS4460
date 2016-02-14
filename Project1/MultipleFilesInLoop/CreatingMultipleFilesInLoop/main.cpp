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
    //Periodic boundary conditions:
    //Check whether particle has gone through side of box every time the position is updated.
    double N_c = 3;
    double b = 1;
    vec r(3);
    r(0) = -7.1;
    r(1) = 1.7;
    r(2) = 3;

    cout << r << endl;

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

    cout << r << endl;
    //Generating multiple files in loop with name: afile1.txt, afile2.txt etc.
    /*
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
    */
    return 0;
}

