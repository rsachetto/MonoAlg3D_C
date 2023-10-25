import sys

def read_acm_file (filename, target_center_x, target_center_y, target_center_z):
    file = open(filename)
    counter_lines = 0
    for line in file:
        if (counter_lines > 1):
            line = line.strip()
            new_line = line.replace(",", " ")
            
            # Cell geometry section
            first_open_bracket_id = new_line.find('[')          # find() -> lowest index
            first_close_bracket_id = new_line.find(']')         
            second_open_bracket_id = new_line.rfind('[')        # rfind() -> hightest index
            second_close_bracket_id = new_line.rfind(']')
            cell_geometry_data = new_line[0:first_open_bracket_id]
            tokens = cell_geometry_data.split()
            center_x, center_y, center_z, dx, dy, dz = float(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5])
            #print("%g %g %g %g %g %g" % (center_x, center_y, center_y, dx, dy, dz))

            # Activation time section
            activation_time_values = []
            activation_time_data = new_line[first_open_bracket_id+1:first_close_bracket_id].split()
            for value in activation_time_data:
                activation_time_values.append(float(value))
            #print(activation_time_values)

            # APD section
            apd_values = []
            apd_data = new_line[second_open_bracket_id+1:second_close_bracket_id].split()
            for value in apd_data:
                apd_values.append(float(value))
            #print(apd_values)

            # SUCESS: Found the target cell!
            if (center_x == target_center_x and center_y == target_center_y and center_z == target_center_z):
                return activation_time_values, apd_values
        
        counter_lines = counter_lines + 1
    file.close()
    return activation_time_values, apd_values

def main():
	
    if len(sys.argv) != 5:
        print("-------------------------------------------------------------------------")
        print("Usage:> python %s <input_file> <center_x> <center_y> <center_z>" % sys.argv[0])
        print("-------------------------------------------------------------------------")
        print("<input_file> = Input \".acm\" file with the activation time and APD data")
        print("<center_x> = Cell center x position")
        print("<center_y> = Cell center y position")
        print("<center_z> = Cell center z position")
        print("-------------------------------------------------------------------------")
        return 1

    input_file = sys.argv[1]
    target_center_x = float(sys.argv[2])
    target_center_y = float(sys.argv[3])
    target_center_z = float(sys.argv[4])
    
    activation_time_values, apd_values = read_acm_file(input_file, target_center_x, target_center_y, target_center_z)
    print("Activation times:")
    print(activation_time_values)
    print()
    print("APDs:")
    print(apd_values)

if __name__ == "__main__":
	main()