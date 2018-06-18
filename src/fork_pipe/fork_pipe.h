
#ifdef _WIN32
void run_child_process_and_process_output (char * program_with_path,  void (*function_to_apply)(void*));
#else
void run_child_process_and_process_output (char * program_with_path,  void (*function_to_apply)(void*), pid_t *pid);
#endif
