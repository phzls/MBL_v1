__author__ = 'liangshengzhang'

import subprocess
import os
import fnmatch

class bcolors:
    # Used for colorful output on consoles
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def para_check(name, var, var_temp, tasks_models):
    """
    Check validity of var_temp according to the type of var.
    :param name: The name of this parameter
    :param var: A base parameter used to give type
    :param var_temp: A parameter whose type and value range are checked
    :param tasks_models: Gives all possible task and model names
    :return: Whether the temp_var is valid and if not, an error message about
    its validity
    """
    valid_var = True
    var_option = None
    if name == "model":
        valid_var = False
        var_option = "Model name not found\n" + "Possible model names:\n"
        for n in tasks_models.models:
            if n == var_temp:
                valid_var = True
            var_option += n + '\n'
    elif name == "task":
        valid_var = False
        var_option = "Task name not found\n" + "Possible task names:\n"
        for n in tasks_models.tasks:
            if n == var_temp:
                valid_var = True
            var_option += n + '\n'
    elif name == "init_func_name":
        valid_var = False
        var_option = "Initial state function name not found\n" + "Possible initial state names:\n"
        for n in tasks_models.inits:
            if n == var_temp:
                valid_var = True
            var_option += n + '\n'
    elif (var == "true") or (var == "false"):
        if (var_temp != "true") and (var_temp != "false"):
            valid_var = False
            var_option = "Argument must be true or false"
    elif var.isdigit():
        try:
            n = int(var_temp)
            if n < 0:
                valid_var = False
                var_option = "Argument must be non-negative integer"
        except ValueError:
            valid_var = False
            var_option = "Argument must be non-negative integer"
    return valid_var, var_option


class DATA(object):
    """
    Record data that is used to generate submitting script
    """
    def __init__(self):
        self.task = None
        self.model = None
        self.num_threads = None
        self.size = None
        self.task_para = "" # String for task parameters
        self.model_para = "" # String for model parameters

class PARA_TEMP_DATA(object):
    """
    Record data when reading parameters which will be used for various decisions when
    generating parameters in the same file
    """
    def __init__(self):
        self.log_time = None # Whether log scale of time is taken in evolution
        self.ent_cal = None # Whether entropy calculation is performed
        self.ent_half_chain = False # Whether need to ask about the half chain size

class RunTime(object):
    """
    Record parameters for calculation, which are estimated running time
    and estimated memory for Feynman cluster
    """
    def __init__(self):
        self.hour = "00"
        self.min = "00"
        self.mem = "500mb"

value_table = {} # Convert some strings to corresponding values
value_table["true"] = True
value_table["false"] = False

def para_gen(filename, tasks_models, data, count, modify_words = None):
    """
    Generate parameters from command line. Possible parameters are first read from file with filename.
    Then these parameters are asked, with default values being the values given in the file. If user
    directly presses enter, then this default value is used. These parameters are then updated with
    same filename.
    :param filename: Filename for the file
    :param tasks_models: Data structure which records possible tasks and models
    :param data: Data struture which records parameters used for generating submitting script
    :param count: Current count of running. Used as suffix for data files
    :param modify_words: Items in the file that possibly needs to be modified if appeared. For other
    item, previous values are used. If it is None, then every item will be asked.
    :return: None
    """
    # Generate parameters from command line
    if count > 1:
        f_old = open("./parameters/" + filename + "_" + str(count-1) + ".dat",'r')
    else:
        f_old = open("./parameters/" + filename + ".dat",'r')

    para_temp_data = PARA_TEMP_DATA()

    task_file = False # Whether this is a task file
    model_file = False # Whether this is a model file

    if (data.task is not None) and (filename == data.task.lower()):
        task_file = True
    elif (data.model is not None) and (filename == data.model.lower()):
        model_file = True

    f_new = open("./parameters/" + filename + "_" + str(count) + ".dat",'w')
    print bcolors.BOLD + '\n' + filename.upper() + bcolors.ENDC
    print bcolors.FAIL + "By directly pressing enter, default/previous values will be used.\n" + bcolors.ENDC

    for line in f_old:
        if (line.startswith("//")) is True:
            f_new.write(line)
        else:
            comment_start = line.find("//")

            # Get lines after //
            comment = " ".join(line[comment_start+2:].split())

            name = line.split()[0]
            var = line.split()[1]

            valid_var = False
            var_option = None

            if modify_words is not None: # Check whether it needs modification
                if name not in modify_words:
                    if name != "left_size":
                        valid_var = True
                        var_temp = var

            if valid_var is False:
                choice = exception(name,data, para_temp_data)
                if choice is not None: # This variable may be skipped
                    valid_var = True
                    if choice == "Previous":
                        var_temp = var
                    else:
                        var_temp = choice
                        valid_var, var_option = para_check(name, var, var_temp, tasks_models)

            while valid_var is not True:
                var_temp = raw_input(name+": (previous: " + var +" )  " + bcolors.OKBLUE +
                                     comment + bcolors.ENDC+ '\n')
                if var_temp == '':
                    var_temp = var
                valid_var, var_option = para_check(name, var, var_temp, tasks_models)
                if valid_var is not True:
                    print bcolors.FAIL + "Error: " + var_option + bcolors.ENDC

            if name == "task":
                data.task = var_temp
            elif name == "model":
                data.model = var_temp
            elif name == "threads_N":
                data.num_threads = var_temp
            elif name == "size":
                data.size = var_temp
            elif name == "log_time":
                para_temp_data.log_time = value_table[var_temp]
            elif name == "Entropy_Per_Model":
                para_temp_data.ent_cal = value_table[var_temp]
                if (modify_words is not None) and ("left_size" not in modify_words):
                    para_temp_data.ent_half_chain = True

            # Obtain strings for parameters
            if (model_file is True):
                if (data.model.find("Flo") > -1) and (name == "J"):
                    data.model_para += "_J_" + str(var_temp)
            elif (task_file is True):
                if data.task.find("Time_Evolution") > -1:
                    if name == "time_step":
                        data.task_para += "_time_step_" + str(var_temp)
                    elif name == "init_func_name":
                        data.task_para += "_init_" + var_temp
                    if data.task.find("Multi") > -1:
                        if name == "model_num":
                            data.task_para += "_model_num_" + str(var_temp)


            f_new.write(name + " " + var_temp + " // " + comment + '\n')

    f_old.close()
    f_new.close()


def slurm_file_gen(data, run_time, count, exe_folder = None):
    """
    Generate a slurm submitting file.
    :param data: Data structure which gives parameters in the submitting script
    :param run_time: Data structure which gives estimated running time
    :param count: The current count of running, which needs to be passed into main program
    :param exe_folder: Executable folder in which the executable lies
    :return: program name
    """
    f_old = open("submit_template.run",'r')

    filename = data.model.lower() + "_" + data.size + data.model_para + "_" \
               + data.task.lower() + data.task_para
    f_new = open(filename + ".run",'w')

    threads = "--ntasks-per-node="
    time = "#SBATCH -t "
    prog = "./"
    dire = "cd "
    server = "Server"

    proc = subprocess.Popen("pwd", stdout=subprocess.PIPE, shell=True)
    address = proc.communicate()[0]

    for line in f_old:
        threads_start = line.find(threads)
        time_start = line.find(time)
        prog_start = line.find(prog)
        dir_start = line.find(dire)
        server_start = line.find(server)

        new_line = line

        new_line_write = True
        if threads_start > -1:
            new_line = line[:threads_start + len(threads)] + data.num_threads + '\n'
        elif time_start > -1:
            new_line = line[:time_start + len(time)] + run_time.hour + ":" + run_time.min + ":00\n"
        elif prog_start > -1:
            if exe_folder is None:
                new_line = line[:prog_start + len(prog)] + filename + " " + str(count) \
                           + " > " + filename + ".out\n"
            else:
                new_line = line[:prog_start + len(prog)] + exe_folder + '/' + filename + " " + str(count) \
                           + " > " + filename + ".out\n"
        elif dir_start > -1:
            new_line = line[:dir_start + len(dire)] + address + '\n'
        elif server_start > -1:
            new_line_write = False

        if new_line_write:
            f_new.write(new_line)

    f_old.close()
    f_new.close()

    return filename


def qsub_file_gen(data, run_time, count, exe_folder = None):
    """
    Generate a qsub submitting file.
    :param data: Data structure which gives parameters in the submitting script
    :param run_time: Data structure which gives estimated running time and memory
    :param count: The current count of running, which needs to be passed into main program
    :param exe_folder: Executable folder in which the executable lies
    :return: program name
    """
    f_old = open("submit_template.run",'r')

    filename = data.model.lower() + "_" + data.size + data.model_para + "_" \
               + data.task.lower() + data.task_para
    f_new = open(filename + ".run",'w')


    threads = "ppn="
    memory = "mem="
    prog = "./"
    dire = "cd "
    server = "Server"

    proc = subprocess.Popen("pwd", stdout=subprocess.PIPE, shell=True)
    address = proc.communicate()[0]

    for line in f_old:
        threads_start = line.find(threads)
        prog_start = line.find(prog)
        dir_start = line.find(dire)
        server_start = line.find(server)
        mem_start = line.find("mem=")

        new_line = line

        new_line_write = True
        if threads_start > -1:
            # Threads and template are in the same line
            new_line = line[:threads_start + len(threads)] + data.num_threads + \
                       ",cput=" + run_time.hour + ":" + run_time.min + ":00\n"
        elif prog_start > -1:
            if exe_folder is None:
                new_line = line[:prog_start + len(prog)] + filename + " " + str(count) \
                           + " > " + filename + ".out\n"
            else:
                new_line = line[:prog_start + len(prog)] + exe_folder + '/' + filename + " " + str(count) \
                           + " > " + filename + ".out\n"
        elif dir_start > -1:
            new_line = line[:dir_start + len(dire)] + address + '\n'
        elif server_start > -1:
            new_line_write = False
        elif mem_start > -1:
            new_line = line[:mem_start + len(memory)] + run_time.mem + "\n"

        if new_line_write:
            f_new.write(new_line)

    f_old.close()
    f_new.close()

    return filename

def exception(name, data, para_temp_data):
    """
    Exception in asking for input according to tasks and models. If it returns not None, then this
    variable will be skipped when process parameter files, and the output value is used.
    :param name: The variable name
    :param data: Data structure which holds tasks and models
    :param para_temp_data: Recording of some parameters used for decision making
    :return: The possible choice for this value. If it is None, then it cannot be skipped. If it
    is "Previous", then the previous value should be used. Otherwise the returned value should be used.
    """
    exc = None
    if data.task is not None:
        if (data.task == "Disorder_Transition") and (name == "J"):
            exc = "Previous"
        if data.model is not None:
           if (data.task == "Disorder_Transition") and (data.model.find("Flo") != -1) and (name == "mid_half_spectrum"):
               exc = "Previous"

    if data.size is not None:
        if name == "left_size":
            if para_temp_data.ent_cal is False:
                exc = "Previous"
            elif para_temp_data.ent_cal is True:
                half_size_choice = None
                if para_temp_data.ent_half_chain:
                    half_size_choice = True
                while half_size_choice is None:
                    ans = raw_input("Compute entanglement entropy for half chain?\n")
                    if ans == "Yes" or ans == "yes" or ans == "Y" or ans == "y" or ans == '':
                        half_size_choice = True
                    elif ans == "No" or ans == "no" or ans == "N" or ans == "n":
                        half_size_choice = False
                    else:
                        print bcolors.FAIL + "Answer must be Yes(Y) or No(N)." + bcolors.ENDC
                if half_size_choice:
                    exc = str(int(data.size)/2)
                    if para_temp_data.ent_half_chain is not True:
                        print "Half chain size: " + exc

    if data.model is not None:
        if (name == "step_size") and (data.model.find("Flo") != -1):
            exc = str(1)
    if (name == "log_time_jump") and (para_temp_data.log_time is False):
        exc = "Previous"
    elif (name == "jump") and (para_temp_data.log_time is True):
        exc = "Previous"

    return exc

def count_obtain():
    """
    Get the count for the current running. This is used to generate a unique number for
    all data files.
    :return: current new count
    """
    file = "./parameters/count.txt"
    file_new = "./parameters/count_temp.txt"
    if os.path.exists(file):
        f_old = open(file,'r')
        f_new = open(file_new,'w')

        for line in f_old:
            count = int(line.split()[0])
        f_old.close()

        count += 1
        f_new.write(str(count)+'\n')
        f_new.close()

    else:
        f = open(file_new,'w')
        count = 1
        f.write(str(count) + '\n')
        f.close()

    return count

def file_clean(count = None, file_modify = None):
    """
    Clean up generated data files and count files when program is
    not run successfully.
    :param count: Given count of this run
    :param file_modify: A map which records whether files have been modified
    :return: None
    """

    filename = "./parameters/count_temp.txt"
    if os.path.exists(filename):
        subprocess.call("rm " + filename, shell=True)

    if count is not None:
        match_name = "_" + str(count) + ".dat"
        for file in os.listdir('./parameters'):
            if fnmatch.fnmatch(file, "*" + match_name):
                start = file.find(match_name)
                file_name = file[:start]
                if file_modify[file_name] is False:
                    new_file = file_name + "_" + str(count-1) + ".dat"
                    subprocess.call("mv ./parameters/" + file +
                                    " ./parameters/" + new_file, shell=True)
                else:
                    subprocess.call("rm ./parameters/" + file, shell=True)

def display_file(filename, data, count):
    """
    Display parameters in the file
    :param filename: Filename for the file
    :param data: Data structure which gives parameters in the submitting script
    :param count: Current count of running. Used as suffix for data files
    :return: None
    """
    # Generate parameters from command line
    if count > 1:
        f = open("./parameters/" + filename + "_" + str(count) + ".dat",'r')
    else:
        f = open("./parameters/" + filename + ".dat",'r')

    para_temp_data = PARA_TEMP_DATA()
    for line in f:
        if (line.startswith("//")) is not True:
            comment_start = line.find("//")

            # Get lines after //
            comment = " ".join(line[comment_start+2:].split())

            name = line.split()[0]
            var = line.split()[1]

            if exception(name, data, para_temp_data) is not "Previous":
                print name + ": " + var + " (" + bcolors.OKBLUE + comment + bcolors.ENDC+ ')\n'

            if name == "log_time":
                para_temp_data.log_time = value_table[var]
            elif name == "Entropy_Per_Model":
                para_temp_data.ent_cal = value_table[var]
                if para_temp_data.ent_cal is True:
                    para_temp_data.ent_cal = None # So that left_size will be directly displayed

    f.close()

def Model_Name_Gen(model):
    """
    Generates filename from the model name
    :param model: Name of the model
    :return: Filename for the corresponding model
    """
    if (model.startswith("XXZ_") and model.endswith("Z_Random_Shift_Real_Flo")):
        return "xxz_general_z_random_shift_real_flo"
    elif (model.startswith("XXZ_") and model.endswith("Random_Field_Shift_Real_Flo")):
        return "xxz_general_random_field_shift_real_flo"
    else:
        return data.model.lower()


