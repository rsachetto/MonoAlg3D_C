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

const double TOLERANCE = 1.0e-02; // 1 cm/s

double calculate_conduction_velocity_from_cable_simulation ()
{
    string filename = "outputs/cable/tissue_activation_time_map_pulse_it_0.vtu";

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
    double point_x_0[3] = {4000, 250, 250};
    double point_x_1[3] = {18000, 250, 250};

    // Find (closest points) Cell indexes for CV computation
    vtkIdType cellId_x_0; // the cell id of the cell containing the closest point
    vtkIdType cellId_x_1; // the cell id of the cell containing the closest point
    
    double closestPoint_x_0[3];
    double closestPoint_x_1[3];
    
    int subId; // not needed
    double closestPointDist2; // not needed
    
    cellLocator->FindClosestPoint(point_x_0, closestPoint_x_0, cellId_x_0, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_x_1, closestPoint_x_1, cellId_x_1, subId, closestPointDist2);
    
    double delta_s_x = sqrt(pow(closestPoint_x_0[0]-closestPoint_x_1[0], 2)); // only changes over x-axis
    
    cout << delta_s_x << endl;
    
    // Read points scalar values
    string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructuredGrid->GetCellData()->GetArray(array_name.c_str()));

    double cv_x = -1.0;

    if(array)
    {
                
        double delta_lat_x = (array->GetValue(cellId_x_1) - array->GetValue(cellId_x_0)); // ms
        
	cout << delta_lat_x << endl;
	cout << array->GetValue(cellId_x_0) << endl;
	cout << array->GetValue(cellId_x_1) << endl;
	
        cv_x = (delta_s_x / delta_lat_x)*0.001;     // {m/s}
    }
    else
    {
        cerr << "[!] ERROR! No Scalar_value found for the points!" << endl;
        exit(EXIT_FAILURE);
    }

    return cv_x;
}

// TODO: Maybe pass a pre-configured config file as an input parameter with the cellular model setup that the user will use
void write_configuration_file (const double sigma)
{
    FILE *file = fopen("/home/Julia/MonoAlg3D_C/scripts/tuneCVbenchmark/configs/cable.ini","w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=6\n");
    fprintf(file,"dt_pde=0.01\n");
    fprintf(file,"simulation_time=100.0\n");
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
    fprintf(file,"output_dir=/home/Julia/MonoAlg3D_C/scripts/tuneCVbenchmark/outputs/cable\n");
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
    
    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_fvm\n");
    fprintf(file,"sigma_x=%g\n",sigma);
    fprintf(file,"sigma_y=%g\n",sigma);
    fprintf(file,"sigma_z=%g\n",sigma);
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
    fprintf(file,"name=Simple Cable\n");
    fprintf(file,"start_dx=500.0\n");
    fprintf(file,"start_dy=500.0\n");
    fprintf(file,"start_dz=500.0\n");
    fprintf(file,"cable_length=20000.0\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_domains.so\n");
    fprintf(file,"main_function=initialize_grid_with_cable_mesh\n");
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
    fprintf(file, "max_x = 500.0\n");
    fprintf(file, "min_y = 0.0\n");
    fprintf(file, "max_y = 500.0\n");
    fprintf(file, "min_z = 0.0\n");
    fprintf(file, "max_z = 3000.0\n");
    fprintf(file,"main_function=stim_x_y_z_limits\n");
    fprintf(file,"library_file=/home/Julia/MonoAlg3D_C/shared_libs/libdefault_stimuli.so\n");
    fprintf(file,"\n");
    
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        cerr << "=============================================================================" << endl;
        cerr << "Usage:> " << argv[0] << " <target_CV>" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "<target_CV> = Target conduction velocity in m/s" << endl;
        cerr << "=============================================================================" << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " 0.67 (Longitudinal normal direction ventricle)" << endl;
        cerr << argv[0] << " 0.33 (Transversal normal direction ventricle)" << endl;
        cerr << argv[0] << " 0.17 (Sheet normal direction ventricle)" << endl;
        cerr << argv[0] << " 1.90 (Purkinje fiber)" << endl;
        cerr << "=============================================================================" << endl;
        
        exit(EXIT_FAILURE);
    }

    double cv, factor;
    double target_cv = atof(argv[1]);
    double sigma = 0.0002;

    do
    {
        write_configuration_file(sigma);
        
        // Run the simulation
        system("/home/Julia/MonoAlg3D_C/bin/MonoAlg3D -c /home/Julia/MonoAlg3D_C/scripts/tuneCVbenchmark/configs/cable.ini");
        
        cv = calculate_conduction_velocity_from_cable_simulation();
        factor = pow(target_cv/cv,2);
        sigma = sigma*factor;

        printf("\n|| Target CV = %g m/s || Computed CV = %g m/s || Factor = %g || Adjusted sigma = %g mS/um ||\n\n",target_cv,cv,factor,sigma);

    }while ( fabs(cv-target_cv) > TOLERANCE );
    
    printf("\n[+] Target conductivity = %g mS/um\n",sigma);

    return 0;
}
