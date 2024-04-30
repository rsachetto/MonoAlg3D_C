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
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkCellLocator.h>

using namespace std;

double* calculate_conduction_velocity_from_benchmark_simulation ()
{
    string filename = "outputs/benchmark/tissue_activation_time_map_pulse_it_0.vtu";

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid* unstructuredGrid = reader->GetOutput();
    uint32_t num_points = unstructuredGrid->GetNumberOfPoints();
    uint32_t num_lines = unstructuredGrid->GetNumberOfCells();

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(unstructuredGrid);
    cellLocator->BuildLocator();

    // Points of interest
    double point_x_0[3] = {4000, 1000, 1000};
    double point_x_1[3] = {18000, 1000, 1000};
    double point_y_0[3] = {1000, 4000, 1000};
    double point_y_1[3] = {1000, 8000, 1000};
    double point_z_0[3] = {1000, 1000, 4000};
    double point_z_1[3] = {1000, 1000, 8000};

    // Find (closest points) Cell indexes for CV computation
    vtkIdType cellId_x_0; // the cell id of the cell containing the closest point
    vtkIdType cellId_x_1; // the cell id of the cell containing the closest point
    vtkIdType cellId_y_0; // the cell id of the cell containing the closest point
    vtkIdType cellId_y_1; // the cell id of the cell containing the closest point
    vtkIdType cellId_z_0; // the cell id of the cell containing the closest point
    vtkIdType cellId_z_1; // the cell id of the cell containing the closest point

    double closestPoint_x_0[3];
    double closestPoint_x_1[3];
    double closestPoint_y_0[3];
    double closestPoint_y_1[3];
    double closestPoint_z_0[3];
    double closestPoint_z_1[3];

    int subId; // not needed
    double closestPointDist2; // not needed
    
    cellLocator->FindClosestPoint(point_x_0, closestPoint_x_0, cellId_x_0, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_x_1, closestPoint_x_1, cellId_x_1, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_y_0, closestPoint_y_0, cellId_y_0, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_y_1, closestPoint_y_1, cellId_y_1, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_z_0, closestPoint_z_0, cellId_z_0, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_z_1, closestPoint_z_1, cellId_z_1, subId, closestPointDist2);

    double delta_s_x = sqrt(pow(closestPoint_x_0[0]-closestPoint_x_1[0], 2)); // only changes over x-axis
    double delta_s_y = sqrt(pow(closestPoint_y_0[1]-closestPoint_y_1[1], 2)); // only changes over y-axis
    double delta_s_z = sqrt(pow(closestPoint_z_0[2]-closestPoint_z_1[2], 2)); // only changes over z-axis

    cout << delta_s_x << endl;
    cout << delta_s_y << endl;
    cout << delta_s_z << endl;

    // Read points scalar values
    string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructuredGrid->GetCellData()->GetArray(array_name.c_str()));

    double cv_x = -1.0;
    double cv_y = -1.0;
    double cv_z = -1.0;

    if(array)
    {
                
        double delta_lat_x = (array->GetValue(cellId_x_1) - array->GetValue(cellId_x_0)); // ms
        double delta_lat_y = (array->GetValue(cellId_y_1) - array->GetValue(cellId_y_0)); // ms
        double delta_lat_z = (array->GetValue(cellId_z_1) - array->GetValue(cellId_z_0)); // ms

	cout << delta_lat_x << endl;
	cout << array->GetValue(cellId_x_0) << endl;
	cout << array->GetValue(cellId_x_1) << endl;
	

        cv_x = (delta_s_x / delta_lat_x)*0.001;     // {m/s}
        cv_y = (delta_s_y / delta_lat_y)*0.001;     // {m/s}
        cv_z = (delta_s_z / delta_lat_z)*0.001;     // {m/s}
    }
    else
    {
        cerr << "[!] ERROR! No Scalar_value found for the points!" << endl;
        exit(EXIT_FAILURE);
    }

    double *cv = new double[3];
    cv[0] = cv_x;
    cv[1] = cv_y;
    cv[2] = cv_z;

    return cv;
}

// TODO: Maybe pass a pre-configured config file as an input parameter with the cellular model setup that the user will use
void write_configuration_file (const double sigma_x, const double sigma_y, const double sigma_z)
{
    FILE *file = fopen("/home/Julia/MonoAlg3D_C/scripts/evaluateBenchmarkCV/configs/benchmark.ini","w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=6\n");
    fprintf(file,"dt_pde=0.01\n");
    fprintf(file,"simulation_time=50.0\n"); // CAREFUL DON'T USE A VALUE THAT'S TOO SMALL!
    fprintf(file,"abort_on_no_activity=false\n");
    fprintf(file,"use_adaptivity=false\n");
    fprintf(file,"quiet=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[update_monodomain]\n");
    fprintf(file,"main_function=update_monodomain_default\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n");
    fprintf(file,"\n");
   
    // For saving the LATs in a format that can be read for calculating the CVs
    fprintf(file,"[save_result]\n");
    fprintf(file,"print_rate=1\n");
    fprintf(file,"output_dir=/home/Julia/MonoAlg3D_C/scripts/evaluateBenchmarkCV/outputs/benchmark\n");
    fprintf(file,"save_pvd=true\n");
    fprintf(file,"file_prefix=V\n");
    fprintf(file,"save_activation_time=true\n");
    fprintf(file,"save_apd=false\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_save_mesh_purkinje.so\n");
    fprintf(file,"main_function=save_tissue_with_activation_times\n");
    fprintf(file,"init_function=init_save_tissue_with_activation_times\n");
    fprintf(file,"end_function=end_save_tissue_with_activation_times\n");
    fprintf(file,"remove_older_simulation=true\n");
    fprintf(file,"\n");
  
/*
    // For saving the VMs for debugging
    fprintf(file,"[save_result]\n");
    fprintf(file,"print_rate=100\n");
    fprintf(file,"output_dir=/home/Julia/MonoAlg3D_C/scripts/evaluateBenchmarkCV/outputs/benchmark\n");
    fprintf(file,"add_timestamp=false\n");
    fprintf(file,"binary=true\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n");
    fprintf(file,"main_function=save_as_ensight\n");
    fprintf(file,"remove_older_simulation=true\n");
    fprintf(file,"\n");
*/

    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_fvm\n");
    fprintf(file,"sigma_x=%g\n",sigma_x);
    fprintf(file,"sigma_y=%g\n",sigma_y);
    fprintf(file,"sigma_z=%g\n",sigma_z);
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_matrix_assembly.so\n");
    fprintf(file,"main_function=homogeneous_sigma_assembly_matrix\n");
    fprintf(file,"\n");
    
    fprintf(file,"[linear_system_solver]\n");
    fprintf(file,"tolerance=1e-15\n");
    fprintf(file,"use_preconditioner=no\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"max_iterations=200\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n");
    fprintf(file,"init_function=init_conjugate_gradient\n");
    fprintf(file,"end_function=end_conjugate_gradient\n");
    fprintf(file,"main_function=conjugate_gradient\n");
    fprintf(file,"\n");
    
    fprintf(file,"[domain]\n");  
    fprintf(file,"name=N-Version Benchmark\n");
    fprintf(file,"start_discretization=500.0\n");
    fprintf(file,"maximum_discretization=500.0\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_domains.so\n");
    fprintf(file,"main_function=initialize_grid_with_benchmark_mesh\n");
    fprintf(file,"side_length_x=20000\n");
    fprintf(file,"side_length_y=10000\n");
    fprintf(file,"side_length_z=10000\n");
    fprintf(file,"\n");
    
    fprintf(file,"[ode_solver]\n");
    fprintf(file,"dt=0.01\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"gpu_id=0\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libToRORd_fkatp_mixed_endo_mid_epi.so\n");
    fprintf(file,"\n");
    
    fprintf(file,"[stim_benchmark]\n");
    fprintf(file,"start = 0.0\n");
    fprintf(file,"duration = 2.0\n");
    fprintf(file,"current = -50.0\n");
    fprintf(file, "min_x = 0.0\n");
    fprintf(file, "max_x = 3000.0\n");
    fprintf(file, "min_y = 0.0\n");
    fprintf(file, "max_y = 3000.0\n");
    fprintf(file, "min_z = 0.0\n");
    fprintf(file, "max_z = 3000.0\n");
    fprintf(file,"main_function=stim_x_y_z_limits\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_stimuli.so\n");
    fprintf(file,"\n");

    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 6)
    {
        cerr << "=============================================================================" << endl;
        cerr << "Usage:> " << argv[0] << " <target_cv_x>" << " <target_cv_y>" << " <target_cv_z>" << \
 " <sigma_x>" <<  " <sigma_y>" <<   " <sigma_z>" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "<target_CV> = Target conduction velocity in m/s" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "Example:" << endl;
        cerr << argv[1] << " 0.67 (Longitudinal normal direction ventricle)" << endl;
        cerr << argv[2] << " 0.33 (Transversal normal direction ventricle)" << endl;
        cerr << argv[3] << " 0.17 (Sheet normal direction ventricle)" << endl;
        cerr << "=============================================================================" << endl;
        
        exit(EXIT_FAILURE);
    }

    double target_cv_x = atof(argv[1]);
    double target_cv_y = atof(argv[2]);
    double target_cv_z = atof(argv[3]);

    double sigma_x = atof(argv[4]);
    double sigma_y = atof(argv[5]);
    double sigma_z = atof(argv[6]);

    write_configuration_file(sigma_x, sigma_y, sigma_z);
        
    // Run the simulation
    system("/home/Julia/MonoAlg3D_C/bin/MonoAlg3D -c /home/Julia/MonoAlg3D_C/scripts/evaluateBenchmarkCV/configs/benchmark.ini");
    
    double* cv;
    cv = calculate_conduction_velocity_from_benchmark_simulation();

    double cv_x, cv_y, cv_z;
    cv_x = cv[0];
    cv_y = cv[1];
    cv_z = cv[2];
    
    printf("\n|| Target CV_x = %g m/s || Computed CV_x = %g m/s || Error CV_x = %g m/s \n\n", target_cv_x,cv_x, cv_x-target_cv_x);
    printf("\n|| Target CV_y = %g m/s || Computed CV_y = %g m/s || Error CV_y = %g m/s \n\n", target_cv_y,cv_y, cv_y-target_cv_y);
    printf("\n|| Target CV_z = %g m/s || Computed CV_z = %g m/s || Error CV_z = %g m/s \n\n", target_cv_z,cv_z, cv_z-target_cv_z);


    return 0;
}
