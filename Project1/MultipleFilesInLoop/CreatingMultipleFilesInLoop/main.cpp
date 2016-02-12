#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
using namespace std;


int main()
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
    return 0;
}

