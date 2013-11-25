#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>

//In the code given below:
//t stands for protected,
//v stands for private,
//m stands for member-variable
//p stands for pointer.

const double PI = 3.141592653589793238462;

typedef struct {
    int numPixels;
    int numProjections;
    int numSlices;
    int imgN;
    int imgLen;
    int detLen;
    } reconVars;

class art {
public:
    art(reconVars* params,
        float *pAngle,
        float *pCenter,
        float *pInput,
        float *pOutput);
    void reconstruct();

private:
    float *vpAngle;
    float *vpCenter;
    float *vpInput;
    float *vpOutput;

    reconVars *vpVars;
    int vNumPixels;
    int vNumSlices;
    int vNumProjections;
    int vImgN;
    int vImgLen;
    int vDetLen;
    };
