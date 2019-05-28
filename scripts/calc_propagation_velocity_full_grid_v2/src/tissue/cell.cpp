#include "cell.h"

bool is_inside_bounds (const double center[],\
                    const double min_x, const double max_x,\
                    const double min_y, const double max_y,\
                    const double min_z, const double max_z)
{
    bool ignore = center[0] < min_x || center[0] > max_x ||\
                center[1] < min_y || center[1] > max_y ||\
                center[2] < min_z || center[2] > max_z;
    
    return !ignore;
}

Cell::Cell ()
{
    this->at = 0.0;
    this->cv = 0.0;
    
    this->in_boundary_1 = false;
    this->in_boundary_2 = false;
    this->in_boundary_3 = false;
    this->in_boundary_4 = false;
}

Cell::Cell (const uint32_t id, const double x, const double y, const double z)
{
    this->id = id;

    this->x = x;
    this->y = y;
    this->z = z;

    this->at = 0.0;
    this->cv = 0.0;

    this->in_boundary_1 = false;
    this->in_boundary_2 = false;
    this->in_boundary_3 = false;
    this->in_boundary_4 = false;
}

Cell::Cell (const uint32_t id, const double x, const double y, const double z,\
        const double p0[], const double p1[], const double p2[], const double p3[], const double p4[],\
        const double p5[], const double p6[], const double p7[])
{
    this->id = id;

    this->x = x;
    this->y = y;
    this->z = z;

    this->at = 0.0;
    this->cv = 0.0;

    memcpy(this->p0,p0,sizeof(double)*3);
    memcpy(this->p1,p1,sizeof(double)*3);
    memcpy(this->p2,p2,sizeof(double)*3);
    memcpy(this->p3,p3,sizeof(double)*3);
    memcpy(this->p4,p4,sizeof(double)*3);
    memcpy(this->p5,p5,sizeof(double)*3);
    memcpy(this->p6,p6,sizeof(double)*3);
    memcpy(this->p7,p7,sizeof(double)*3);

    this->in_boundary_1 = false;
    this->in_boundary_2 = false;
    this->in_boundary_3 = false;
    this->in_boundary_4 = false;
}

void Cell::setCenter (const double center[])
{
    this->x = center[0];
    this->y = center[1];
    this->z = center[2];
}

void Cell::setVertex (const double p0[], const double p1[], const double p2[], const double p3[],\
                    const double p4[], const double p5[], const double p6[], const double p7[])
{
    memcpy(this->p0,p0,sizeof(double)*3);
    memcpy(this->p1,p1,sizeof(double)*3);
    memcpy(this->p2,p2,sizeof(double)*3);
    memcpy(this->p2,p2,sizeof(double)*3);
    memcpy(this->p3,p3,sizeof(double)*3);
    memcpy(this->p4,p4,sizeof(double)*3);
    memcpy(this->p5,p5,sizeof(double)*3);
    memcpy(this->p6,p6,sizeof(double)*3);
    memcpy(this->p7,p7,sizeof(double)*3);
}

void Cell::setActivationTime (const double value)
{
    this->at = value;
}

void Cell::setIndex (const uint32_t id)
{
    this->id = id;
}

void Cell::setBoundaries (const uint32_t i, const uint32_t j, const uint32_t nx, const uint32_t ny)
{
    // Inside east boundary
    if (i == 0)
        this->in_boundary_1 = true;
    // Inside south boundary
    if (j == (ny-1))
        this->in_boundary_2 = true;
    // Inside west boundary
    if (i == (nx-1))
        this->in_boundary_3 = true;
    // Inside north boundary
    if (j == 0)
        this->in_boundary_4 = true;
}

void Cell::print ()
{
    printf("****** Cell %u ******\n",id);
    printf("Center = (%g,%g,%g)\n",x,y,z);
    printf("Activation time = %g\n",at);
    printf("Conduction velocity = %g\n",cv);
    printf("Boundary 1 = %u\n",in_boundary_1);
    printf("Boundary 2 = %u\n",in_boundary_2);
    printf("Boundary 3 = %u\n",in_boundary_3);
    printf("Boundary 4 = %u\n",in_boundary_4);
    printf("------ Vertex -------\n");
    printf("\tVertex 0 = (%g,%g,%g)\n",p0[0],p0[1],p0[2]);
    printf("\tVertex 1 = (%g,%g,%g)\n",p1[0],p1[1],p1[2]);
    printf("\tVertex 2 = (%g,%g,%g)\n",p2[0],p2[1],p2[2]);
    printf("\tVertex 3 = (%g,%g,%g)\n",p3[0],p3[1],p3[2]);
    printf("\tVertex 4 = (%g,%g,%g)\n",p4[0],p4[1],p4[2]);
    printf("\tVertex 5 = (%g,%g,%g)\n",p5[0],p5[1],p5[2]);
    printf("\tVertex 6 = (%g,%g,%g)\n",p6[0],p6[1],p6[2]);
    printf("\tVertex 7 = (%g,%g,%g)\n",p7[0],p7[1],p7[2]);
    printf("---------------------\n");
    printf("*********************\n");
}

void print_cells (const Cell *cells, const uint32_t num_cells)
{
    for (uint32_t i = 0; i < num_cells; i++)
    {
        printf("Cell %u -- (%g,%g,%g)\n",cells[i].id,cells[i].x,cells[i].y,cells[i].z);
    }
}