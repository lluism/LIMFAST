#include "SOURCES.H"
#include <stdio.h>

double test_fstar(double z, double M){return 0.1;}
double test_fstarIII(double z, double M){return 0.1;}
double test_fesc(double z, double M){return 0.1;}
double test_fescIII(double z, double M){return 0.1;}
double test_Nion(double z, double M){return 4000;}
double test_NionIII(double z, double M){return 4000;}
double test_minMass(double z){return 1.0e8;}
double test_minMassIII(double z){return -1;}

int main()
{
    init_ps();
    sources s;
    double res;
    //s = setSources(test_fesc, test_fescIII, test_Nion, test_NionIII, test_fstar, test_fstarIII, test_minMass, test_minMassIII);
    s = defaultSources();
    res = ionEff(7, s);
    printf("%le\n", res);
    return 0;
}
