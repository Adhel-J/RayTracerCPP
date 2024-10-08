
#include <iostream>
using namespace std;





int main(int argc, char const *argv[])
{
    constexpr int Cw = 300;
    constexpr int Ch = 300;

    int recursionDepth = 3;

    constexpr int positionx = 0;
    constexpr int positiony = 0;
    constexpr int positionz = 0;
    constexpr int anglex = 0;
    constexpr int angley = 0;
    constexpr int anglez = 0;

    for (int x = 0; x < Cw; x++){
        for (int y = 0; y < Ch; y++){
            std::cout << x << ", " << y << endl;
        }
    }

    return 0;
}


//
// Created by adhel on 10/8/2024.
//
