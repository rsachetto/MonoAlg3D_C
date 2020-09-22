#ifndef CALC_PROPAGATION_VELOCITY_UTILS_H
#define CALC_PROPAGATION_VELOCITY_UTILS_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>

#include "file_utils.h"

#define PRINT_LINE "====================================================================================================================="

void usage (const char pname[]);

class Point_3D
{
private:
    double x, y, z;
public:
    Point_3D (const double x, const double y, const double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    bool operator <(const Point_3D& pt) const
    {
        return (this->x < pt.x) || \
               ((!(pt.x < this->x)) && (this->y < pt.y)) || \
               ((!(pt.x < this->x)) && (!(pt.y < this->y)) && (this->z < pt.z));
    }
    friend std::ostream& operator<<(std::ostream& os, const Point_3D& obj);
};

#endif //MONOALG3D_UTILS_H_H