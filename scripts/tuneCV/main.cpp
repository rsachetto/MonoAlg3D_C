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

using namespace std;

const double TOLERANCE = 1.0e-03;

double calculate_conduction_velocity_from_simulation ()
{
    string filename = "outputs/cable/purkinje_activation_time_map_pulse_it_0.vtp";

    // Read all the data from the file
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkPolyData* polydata = reader->GetOutput();
    uint32_t num_points = polydata->GetNumberOfPoints();
    uint32_t num_lines = polydata->GetNumberOfLines();

    // Read points scalar values
    string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray(array_name.c_str()));

    double cv = 0.0;
    if(array)
    {
        // Cell indexes for CV computation
        uint32_t middle_cell = num_points/2;
        uint32_t prev_cell = middle_cell-25;
        uint32_t next_cell = middle_cell+25;
        
        double delta_s = 5000.0;
        double delta_t = (array->GetValue(next_cell) - array->GetValue(prev_cell));
        cv = (delta_s / delta_t)*0.001;     // {m/s}
    }
    else
    {
        cerr << "[!] ERROR! No Scalar_value found for the points!" << endl;
        exit(EXIT_FAILURE);
    }

    return cv;
}

// TODO: Maybe pass a pre-configured config file as an input parameter with the cellular model setup that the user will use
void write_configuration_file (const double sigma)
{
    FILE *file = fopen("/home/berg/Github/MonoAlg3D_C/scripts/tuneCV/configs/cable.ini","w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=4\n");
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
    fprintf(file,"print_rate=1\n");
    fprintf(file,"output_dir=/home/berg/Github/MonoAlg3D_C/scripts/tuneCV/outputs/cable\n");
    fprintf(file,"save_pvd=true\n");
    fprintf(file,"file_prefix=V\n");
    fprintf(file,"save_activation_time=true\n");
    fprintf(file,"save_apd=true\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_save_mesh_purkinje.so\n");
    fprintf(file,"main_function=save_purkinje_with_activation_times\n");
    fprintf(file,"init_function=init_save_purkinje_with_activation_times\n");
    fprintf(file,"end_function=end_save_purkinje_with_activation_times\n");
    fprintf(file,"remove_older_simulation=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_fvm\n");
    fprintf(file,"sigma_purkinje=%g\n",sigma);
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libpurkinje_matrix_assembly.so\n");
    fprintf(file,"main_function=purkinje_fibers_assembly_matrix\n");
    fprintf(file,"\n");
    
    fprintf(file,"[linear_system_solver]\n");
    fprintf(file,"tolerance=1e-16\n");
    fprintf(file,"use_preconditioner=yes\n");
    fprintf(file,"max_iterations=500\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n");
    fprintf(file,"main_function=conjugate_gradient\n");
    fprintf(file,"\n");
    
    fprintf(file,"[purkinje]\n");
    fprintf(file,"name=Simple Purkinje\n");
    //fprintf(file,"dx=100.0\n");
    fprintf(file,"dx=400.0\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n");
    fprintf(file,"main_function=initialize_purkinje_with_custom_mesh\n");
    //fprintf(file,"network_file=/home/berg/Github/MonoAlg3D_C/networks/simple_cable_10cm.vtk\n");
    fprintf(file,"network_file=/home/berg/Github/MonoAlg3D_C/networks/simple_cable_40cm.vtk\n");
    fprintf(file,"\n");
    
    fprintf(file,"[ode_solver]\n");
    //fprintf(file,"adaptive=true\n");
    fprintf(file,"dt=0.02\n");
    fprintf(file,"use_gpu=no\n");
    fprintf(file,"gpu_id=0\n");
    //fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libtrovato_2019.so\n");
    //fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libToRORd_fkatp_endo.so\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libten_tusscher_3_endo.so\n");
    fprintf(file,"\n");
    
    fprintf(file,"[purkinje_stim_his]\n");
    fprintf(file,"start = 0.0\n");
    fprintf(file,"duration = 5.0\n");
    //fprintf(file,"current = -40.0\n");
    fprintf(file,"current = -53.0\n");
    fprintf(file,"id_limit = 25\n");
    fprintf(file,"main_function=stim_if_id_less_than\n");
    fprintf(file,"library_file=/home/berg/Github/MonoAlg3D_C/shared_libs/libdefault_stimuli.so\n");
    fprintf(file,"\n");

    fclose(file);
}

double calculate_error (const double cv, const double target_cv)
{
    return sqrt(pow(cv,2)-pow(target_cv,2))/target_cv;
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
        system("/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D -c /home/berg/Github/MonoAlg3D_C/scripts/tuneCV/configs/cable.ini");
        
        cv = calculate_conduction_velocity_from_simulation();
        factor = pow(target_cv/cv,2);
        sigma = sigma*factor;

        printf("\n|| Target CV = %g m/s || Computed CV = %g m/s || Factor = %g || Adjusted sigma = %g mS/um ||\n\n",target_cv,cv,factor,sigma);

    }while ( fabs(cv-target_cv) > TOLERANCE );
    
    printf("\n[+] Target conductivity = %g mS/um\n",sigma);

    return 0;
}
