#ifndef HEXAEDRON_CELL_H
#define HEXAEDRON_CELL_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>



class HexaedronCell
{
private:
  double p0[3];
  double p1[3];
  double p2[3];
  double p3[3];
  double p4[3];
  double p5[3];
  double p6[3];
  double p7[3];

public:
  HexaedronCell (double p0[], double p1[], double p2[], double p3[], double p4[], double p5[], double p6[], double p7[])
  {
    memcpy(this->p0,p0,sizeof(double)*3);
    memcpy(this->p1,p1,sizeof(double)*3);
    memcpy(this->p2,p2,sizeof(double)*3);
    memcpy(this->p3,p3,sizeof(double)*3);
    memcpy(this->p3,p3,sizeof(double)*3);
    memcpy(this->p4,p4,sizeof(double)*3);
    memcpy(this->p5,p5,sizeof(double)*3);
    memcpy(this->p6,p6,sizeof(double)*3);
    memcpy(this->p7,p7,sizeof(double)*3);
  }

};


#endif //MONOALG3D_UTILS_H_H
