#include <iostream>
#include "program_runs.hpp"
//#include "fun_strings.hpp"

/*
###########################################################################################################################
This is all I want to see in main(?)
###########################################################################################################################
*/

int determine_mode(int argc, char const* argv[])
{
    int i = 1;
    int mode = -1;
    while (i < argc)
    {
        std::string s(argv[i]);
        if (s.compare("-m")==0){
            if (i+1 < argc){
                mode = std::stoi(argv[i+1]);
            }
            break;
        }
        i+=1;
    }
    return mode;
}


int main(int argc, char const* argv[])
{
    int mode = determine_mode(argc, argv);

    if (mode == -1)
    {
        std::cout << "No legal mode detected\n";
    } 
    else if (mode == 2)
    {
        run_mode_2(argc, argv);
    } 
    else
    {
        std::cout << "No legal mode detected\n";
    } 
    std::cout << "Program run finished, thx 4 using and cya nextime\n";
}

/*
###########################################################################################################################
Old code below, copy from there to new locations
###########################################################################################################################
*/
