#include "purkinje.h"

struct purkinje_network* new_purkinje_network ()
{
    struct purkinje_network *pk = (struct purkinje_network*)malloc(sizeof(struct purkinje_network));

    pk->num_points = 0;
    pk->num_lines = 0;

    pk->points = NULL;
    pk->lines = NULL;

    return pk;
}

void read_purkinje_network_from_vtp (struct purkinje_network *the_purkinje_network, vtkPolyData *polydata)
{
    uint32_t num_points = polydata->GetNumberOfPoints();
    uint32_t num_cells = polydata->GetNumberOfCells();

    the_purkinje_network->points = (struct point*)malloc(sizeof(struct point)*num_points);
    the_purkinje_network->lines = (struct line*)malloc(sizeof(struct line)*num_cells);
    the_purkinje_network->point_data = (double*)malloc(sizeof(double)*num_points);

    the_purkinje_network->num_points = num_points;
    the_purkinje_network->num_lines = num_cells;

    struct point *the_points = the_purkinje_network->points;
    struct line *the_lines = the_purkinje_network->lines;
    double *the_data = the_purkinje_network->point_data;

    // Read the points
    for (uint32_t i = 0; i < num_points; i++)
    {
        double p[3];
        polydata->GetPoint(i,p);
        
        the_points[i].id = i;
        the_points[i].x = p[0];
        the_points[i].y = p[1];
        the_points[i].z = p[2];
    }

    // Read the cells
    for (uint32_t i = 0; i < num_cells; i++)
    {
        vtkCell *cell = polydata->GetCell(i);

        vtkLine *line = dynamic_cast<vtkLine*>(cell);

        uint32_t src = line->GetPointIds()->GetId(0);
        uint32_t dest = line->GetPointIds()->GetId(1);
        
        the_lines[i].id = i;
        the_lines[i].src = src;
        the_lines[i].dest = dest;
    }    

    // Read point data
    std::string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray(array_name.c_str()));

    if(array)
    {
        //std::cout << "Got array " << array_name << " with " << num_cells << " values" << std::endl;
        for (int i = 0; i < num_points; i++)
        {
            double value;
            value = array->GetValue(i);

            the_data[i] = value;
            //std::cout << i << ": " << value << std::endl;
        }
    }
    else
    {
        std::cout << " Not found a PointData array named " << array_name << std::endl;
    }
}

void print_purkinje_network (struct purkinje_network *pk)
{
    uint32_t np = pk->num_points;
    uint32_t nl = pk->num_lines;

    for (uint32_t i = 0; i < np; i++)
        std::cout << "Point " << i << " = (" << pk->points[i].x << "," << pk->points[i].y << "," << pk->points[i].z << ") || data = " << pk->point_data[i] << std::endl;

    for (uint32_t i = 0; i < nl; i++)
        std::cout << "Line " << i << " = (" << pk->lines[i].src << "," << pk->lines[i].dest << ")" << std::endl;     


    std::cout << "There are " << pk->num_points << std::endl;
    std::cout << "There are " << pk->num_lines << std::endl;

}

