// Author: Lucas Berg

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkCellLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPointData.h>

using namespace std;

// Return the tissue conduction velocity: (CV_x, CV_y, CV_z)
double* calculate_conduction_velocity_from_simulation ()
{
    string filename = "outputs/purkinje_cuboid/tissue_activation_time_map_pulse_it_0.vtu";

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid* unstructured_grid = reader->GetOutput();
    uint32_t num_points = unstructured_grid->GetNumberOfPoints();
    uint32_t num_cells = unstructured_grid->GetNumberOfCells();

    // Read points scalar values
    string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(unstructured_grid);
    cellLocator->BuildLocator();

    double p1[3] = {5000,5000,5000};
    double p2[3] = {5000,9000,5000};
    double p3[3] = {5000,1000,5000};
    double p4[3] = {9000,5000,5000};
    double p5[3] = {1000,5000,5000};
    double p6[3] = {5000,5000,9000};
    double p7[3] = {5000,5000,1000};

    double closestPoint[3];
    double closestPointDist;
    vtkIdType cellId[7];
    int subId;

    cellLocator->FindClosestPoint(p1,closestPoint,cellId[0],subId,closestPointDist);
    cellLocator->FindClosestPoint(p2,closestPoint,cellId[1],subId,closestPointDist);
    cellLocator->FindClosestPoint(p3,closestPoint,cellId[2],subId,closestPointDist);
    cellLocator->FindClosestPoint(p4,closestPoint,cellId[3],subId,closestPointDist);
    cellLocator->FindClosestPoint(p5,closestPoint,cellId[4],subId,closestPointDist);
    cellLocator->FindClosestPoint(p6,closestPoint,cellId[5],subId,closestPointDist);
    cellLocator->FindClosestPoint(p7,closestPoint,cellId[6],subId,closestPointDist);

    double *cv = new double[3]();
    if (array)
    {
        const double delta_s = 4000.0;

        // X
        cv[0] = delta_s / (array->GetValue(cellId[3]) - array->GetValue(cellId[0]))*0.001; // m/s
        // Y
        cv[1] = delta_s / (array->GetValue(cellId[1]) - array->GetValue(cellId[0]))*0.001; // m/s
        // Z
        cv[2] = delta_s / (array->GetValue(cellId[5]) - array->GetValue(cellId[0]))*0.001; // m/s
    }
    else
    {
        delete [] cv;
        cerr << "[!] ERROR! No Scalar_value found for the points!" << endl;
        exit(EXIT_FAILURE);
    }

    return cv;
}

// Purkinje coupling cuboid simulation
void write_configuration_file (const double rpmj)
{
    FILE *file = fopen("/home/berg/Github/MonoAlg3D_C/scripts/tunePMJ/configs/purkinje_cuboid.ini","w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=8\n");
    fprintf(file,"dt_pde=0.02\n");
    fprintf(file,"simulation_time=150.0\n");
    fprintf(file,"abort_on_no_activity=false\n");
    fprintf(file,"use_adaptivity=false\n");
    fprintf(file,"quiet=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[update_monodomain]\n");
    fprintf(file,"main_function=update_monodomain_default\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n");
    fprintf(file,"\n");
    
    fprintf(file,"[save_result]\n");
    fprintf(file,"print_rate=10\n");
    fprintf(file,"output_dir=/home/berg/Github/MonoAlg3D_C/scripts/tunePMJ/outputs/purkinje_cuboid\n");
    fprintf(file,"save_pvd=true\n");
    fprintf(file,"file_prefix=V\n");
    fprintf(file,"save_activation_time=true\n");
    fprintf(file,"save_apd=false\n");
    fprintf(file,"save_purkinje_velocity=true\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_save_mesh_purkinje.so\n");
    fprintf(file,"main_function=save_purkinje_coupling_with_activation_times\n");
    fprintf(file,"init_function=init_save_purkinje_coupling_with_activation_times\n");
    fprintf(file,"end_function=end_save_purkinje_coupling_with_activation_times\n");
    fprintf(file,"remove_older_simulation=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_coupling_fvm\n");
    fprintf(file,"sigma_x=0.0001334\n");
    fprintf(file,"sigma_y=0.0001334\n");
    fprintf(file,"sigma_z=0.0001334\n");
    //fprintf(file,"sigma_x=0.00015\n");
    //fprintf(file,"sigma_y=0.000085\n");
    //fprintf(file,"sigma_z=0.000075\n");
    fprintf(file,"sigma_purkinje=0.002567\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libpurkinje_coupling_matrix_assembly.so\n");
    fprintf(file,"main_function=purkinje_coupling_assembly_matrix\n");
    fprintf(file,"\n");
    
    fprintf(file,"[linear_system_solver]\n");
    fprintf(file,"tolerance=1e-16\n");
    fprintf(file,"use_preconditioner=no\n");
    fprintf(file,"max_iterations=500\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n");
    fprintf(file,"main_function=conjugate_gradient\n");
    fprintf(file,"init_function=init_conjugate_gradient\n");
    fprintf(file,"end_function=end_conjugate_gradient\n");
    fprintf(file,"\n");

    fprintf(file,"[purkinje_linear_system_solver]\n");
    fprintf(file,"tolerance=1e-16\n");
    fprintf(file,"use_preconditioner=no\n");
    fprintf(file,"max_iterations=500\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n");
    fprintf(file,"main_function=conjugate_gradient\n");
    fprintf(file,"init_function=init_conjugate_gradient\n");
    fprintf(file,"end_function=end_conjugate_gradient\n");
    fprintf(file,"\n");
    
    fprintf(file,"[ode_solver]\n");
    //fprintf(file,"adaptive=true\n");
    fprintf(file,"dt=0.02\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"gpu_id=0\n");
    //fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libtrovato_2019.so\n");
    //fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libToRORd_fkatp_endo.so\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libten_tusscher_3_endo.so\n");
    fprintf(file,"\n");

    fprintf(file,"[purkinje_ode_solver]\n");
    fprintf(file,"adaptive=true\n");
    fprintf(file,"use_gpu=no\n");
    fprintf(file,"gpu_id=0\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libtrovato_2019.so\n");
    fprintf(file,"\n");

    fprintf(file,"[purkinje]\n");
    fprintf(file,"name=Purkinje cable\n");
    fprintf(file,"dx=100.0\n");
    fprintf(file,"rpmj=%g\n",rpmj);
    fprintf(file,"asymm_ratio=1.0\n");
    fprintf(file,"pmj_scale=800.0\n");
    fprintf(file,"nmin_pmj=60\n");
    fprintf(file,"nmax_pmj=65\n");
    fprintf(file,"retro_propagation=true\n");
    fprintf(file,"network_file=/home/berg/Github/MonoAlg3D_C/local_networks/shocker_paper/cable_10cm.vtk\n");
    fprintf(file,"pmj_location_file=/home/berg/Github/MonoAlg3D_C/local_networks/shocker_paper/cable_pmj.vtk\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n");
    fprintf(file,"main_function=initialize_purkinje_coupling_with_custom_mesh\n");
    fprintf(file,"\n");
    
    fprintf(file,"[domain]\n");
    fprintf(file,"name=Cuboid Mesh\n");
    fprintf(file,"start_dx=400.0\n");
    fprintf(file,"start_dy=400.0\n");
    fprintf(file,"start_dz=400.0\n");
    fprintf(file,"side_length_x=10000\n");
    fprintf(file,"side_length_y=10000\n");
    fprintf(file,"side_length_z=10000\n");
    fprintf(file,"main_function=initialize_grid_with_cuboid_mesh\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_domains.so\n");
    fprintf(file,"\n");

    fprintf(file,"[purkinje_stim_his]\n");
    fprintf(file,"start = 0.0\n");
    fprintf(file,"duration = 5.0\n");
    fprintf(file,"current = -40.0\n");
    fprintf(file,"id_limit = 25\n");
    fprintf(file,"main_function=stim_if_id_less_than\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_stimuli.so\n");
    fprintf(file,"\n");

    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        cerr << "=============================================================================" << endl;
        cerr << "Usage:> " << argv[0] << " <R_PMJ>" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "<R_PMJ> = PMJ resistance in M.ohm" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " 1.0" << endl;
        cerr << "=============================================================================" << endl;
        
        exit(EXIT_FAILURE);
    }

    double rpmj = atof(argv[1]);
    
    write_configuration_file(rpmj);
        
    // Run the simulation
    system("/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D -c /home/berg/Github/MonoAlg3D_C/scripts/tunePMJ/configs/purkinje_cuboid.ini");
    
    double *cv = calculate_conduction_velocity_from_simulation();

    printf("|| CV_x = %g m/s || CV_y = %g m/s || CV_z = %g m/s ||\n",cv[0],cv[1],cv[2]);

    delete [] cv;

    return 0;
}
