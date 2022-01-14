#!/bin/bash

tar -zcf custom_files.tar.gz ../src/domains_library/custom_* ../src/extra_data_library/custom_*
tar -zcf local_files.tar.gz ../local_meshes ../local_networks ../private_configs ../private_models
