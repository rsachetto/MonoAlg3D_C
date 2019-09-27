#include "file_utils.h"

void read_vtu_files_from_folder(const std::string folder_path,\
                            std::vector<struct vtu_file> &vtu_files)
{
    DIR *dir;
    std::string extension_name = "vtu";

    dir = opendir(folder_path.c_str());
    if ( dir )
    {
        printf("[+] Folder '%s' has been open !\n",folder_path.c_str());

        // Filter files with a specific extension and sort them by their number
        filter_files_inside_folder_by_extension(dir,folder_path,extension_name,vtu_files);
        
        // Sort the files inside the vector by their timestep number
        std::sort(vtu_files.begin(),vtu_files.end(),sort_by_timestep);

    }
    else
    {
        printf("[-] ERROR! Cannot open '%s' folder !\n",folder_path.c_str());
    }
    closedir(dir);
}

void filter_files_inside_folder_by_extension (DIR *dir,\
                            const std::string folder_path,const std::string extension_name,\
                            std::vector<struct vtu_file> &vtu_files)
{
    struct dirent *ent;

    while ( (ent = readdir(dir)) )
    {
        std::string filename = ent->d_name;
        if (filename.size() > 3)    // Avoid the '.' and '..'
        {
            std::string extension = filename.substr(filename.size()-3,filename.size()-1);

            // Insert this file into the vector
            if (extension_name == extension && filename != "activation-map.vtu")
            {
                struct vtu_file file;
                
                file.name = filename;
                file.timestep = get_timestep_from_filename(filename);

                vtu_files.push_back(file);
            }
                

        }
    }
}

uint32_t get_timestep_from_filename (const std::string filename)
{
    uint32_t start, end, offset;

    start = filename.find_last_of("_");
    end = filename.find_first_of(".");
    offset = end - start - 1;

    std::string number = filename.substr(start+1,offset);
    
    size_t aux;
    return stoi(number,&aux);
}

bool sort_by_timestep (struct vtu_file a, struct vtu_file b) 
{ 
    return a.timestep < b.timestep; 
}