// Author: Lucas Berg
// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Program that reads the LAT of the tissue and Purkinje network and generates candidate points to be new PMJs from the Purkinje terminals.
// ------------------------------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "purkinje/purkinje.h"
#include "graph/graph.h"

#define PRINT_LINE "=================================================================================================================="

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <vtp_file> <vtu_file> <root_nodes_file>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<vtp_file> = VTP file with the Purkinje network LAT" << endl;
        cerr << "<vtu_file> = VTU file with the tissue LAT" << endl;
        cerr << "<root_nodes_file> = Root nodes coordinate positions" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " inputs/purkinje_LAT.vtp inputs/tissue_LAT.vtu inputs/root_nodes_coordinates.txt" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string pk_filename = argv[1];
    string tiss_filename = argv[2];
    string root_nodes_filename = argv[3];

    // Read Purkinje data
    printf("[purkinje] Reading Purkinje data from '%s'\n",pk_filename.c_str());
    vtkSmartPointer<vtkXMLPolyDataReader> pk_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    pk_reader->SetFileName(pk_filename.c_str());
    pk_reader->Update();

    vtkPolyData* polydata = pk_reader->GetOutput();

    printf("[purkinje] Parsing Purkinje network to a graph data structure\n");
    struct purkinje_network *the_purkinje_network = new_purkinje_network();
    read_purkinje_network_from_vtp(the_purkinje_network,polydata);
    //print_purkinje_network(the_purkinje_network);

    struct graph *the_graph = new_graph();
    build_graph_from_purkinje_network(the_graph,the_purkinje_network);

    printf("[purkinje] Calculating Purkinje network terminals\n");
    uint32_t num_terminals;
    struct terminal *the_terminals;
    the_terminals = calculate_terminals(the_graph,num_terminals);
    //print_terminals(the_terminals,num_terminals);

    // Read Tissue data
    printf("[tissue] Reading Tissue data from '%s'\n",tiss_filename.c_str());
    vtkSmartPointer<vtkXMLUnstructuredGridReader> tiss_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    tiss_reader->SetFileName(tiss_filename.c_str());
    tiss_reader->Update();
    vtkUnstructuredGrid *ugrid = tiss_reader->GetOutput();
    string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> lat_array = vtkFloatArray::SafeDownCast(ugrid->GetCellData()->GetArray(array_name.c_str()));
    if(!lat_array)
    {
        fprintf(stderr,"[-] ERROR! Cannot find CellData '%s'\n",array_name.c_str());
        exit(EXIT_FAILURE);
    }

    // Read root nodes data
    printf("[root_nodes] Reading RootNodes data from '%s'\n",root_nodes_filename.c_str());
    uint32_t num_root_nodes;
    struct terminal *root_nodes = read_root_nodes(root_nodes_filename.c_str(),num_root_nodes);
    //print_terminals(root_nodes,num_root_nodes);

    printf("[pmj] Building CellLocator for the tissue mesh\n");
    vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
    cell_locator->SetDataSet(ugrid);
    cell_locator->BuildLocator();

    // Search for the closest tissue cell to each root node
    double max_lat_root_node = __DBL_MIN__;
    for (uint32_t i = 0; i < num_root_nodes; i++)
    {
        int subId;
        vtkIdType cellId;
        double closest_dist;
        double pos[3], closest_pos[3];
        pos[0] = root_nodes[i].x;
        pos[1] = root_nodes[i].y;
        pos[2] = root_nodes[i].z;

        cell_locator->FindClosestPoint(pos,closest_pos,cellId,subId,closest_dist);
        root_nodes[i].value = lat_array->GetValue(cellId);
        if (lat_array->GetValue(cellId) > max_lat_root_node)
            max_lat_root_node = lat_array->GetValue(cellId);
        
        //printf("[root_node] id = %u // pos = (%g %g %g) // closest_pos = (%g %g %g) // closest_dist = %g // closest_cell_id = %u\n",i,pos[0],pos[1],pos[2],closest_pos[0],closest_pos[1],closest_pos[2],closest_dist,cellId);
    }
    
    // Now, search for the closest tissue cell to each Purkinje terminal node  
    //const double lat_offset = 5.0;    // low
    //const double lat_offset = 10.0;   // medium
    const double lat_offset = 15.0;   // high
    vector<uint32_t> new_pmjs_ids;
    for (uint32_t i = 0; i < num_terminals; i++)
    {
        int subId;
        vtkIdType cellId;
        double closest_dist;
        double pos[3], closest_pos[3];
        pos[0] = the_terminals[i].x;
        pos[1] = the_terminals[i].y;
        pos[2] = the_terminals[i].z;

        cell_locator->FindClosestPoint(pos,closest_pos,cellId,subId,closest_dist);
        double pk_lat = the_terminals[i].value;
        double tiss_lat = lat_array->GetValue(cellId);
        if (pk_lat+lat_offset >= tiss_lat)
        {
            //printf("Terminal %u is activated by the Purkinje+Offset at %g ms and by the tissue at %g ms\n",i,pk_lat+lat_offset,tiss_lat);
            new_pmjs_ids.push_back(i);
        }
    }
    printf("There is a total of %u/%u candidates for PMJs\n",new_pmjs_ids.size(),num_terminals);

    // Fill the new pmjs array
    uint32_t total_num_pmjs = new_pmjs_ids.size();
    struct terminal *the_pmjs = (struct terminal*)malloc(sizeof(struct terminal)*(total_num_pmjs+num_root_nodes));
    for (uint32_t i = 0; i < total_num_pmjs; i++)
    {
        uint32_t id = new_pmjs_ids[i];
        the_pmjs[i].id = i;
        the_pmjs[i].x = the_terminals[id].x;
        the_pmjs[i].y = the_terminals[id].y;
        the_pmjs[i].z = the_terminals[id].z;
        the_pmjs[i].value = the_terminals[id].value;
    }
    // Add the original root nodes in the back of the array
    for (uint32_t i = total_num_pmjs, j = 0; i < total_num_pmjs+num_root_nodes; i++, j++)
    {
        the_pmjs[i].id = i;
        the_pmjs[i].x = root_nodes[j].x;
        the_pmjs[i].y = root_nodes[j].y;
        the_pmjs[i].z = root_nodes[j].z;
        the_pmjs[i].value = root_nodes[j].value;
    }

    // Define a percentage of the candidate points to be selected
    const double percentage = 1.0;
    write_terminals_to_vtk(the_pmjs,total_num_pmjs,num_root_nodes,percentage);

    free_purkinje_network(the_purkinje_network);
    free_graph(the_graph);
    free(the_terminals);
    free(root_nodes);
    free(the_pmjs);

    return 0;
}
