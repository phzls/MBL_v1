__author__ = 'liangshengzhang'

import subprocess

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


class RunTime(object):
    """
    Record parameters for calculation, which are estimated running time
    and estimated memory for Feynman cluster
    """
    def __init__(self):
        self.hour = "00"
        self.min = "00"
        self.mem = "500mb"


def para_gen(filename, tasks_models, data = None):
    """
    Generate parameters from command line. Possible parameters are first read from file with filename.
    Then these parameters are asked, with default values being the values given in the file. If use
    directly presses enter, then this default value is used. These parameters are then updated with
    same filename.
    :param filename: Filename for the file
    :param tasks_models: Data structure which records possible tasks and models
    :param data: Data struture which records parameters used for generating submitting script
    :return: None
    """
    # Generate parameters from command line
    f_old = open("./parameters/" + filename + ".dat",'r')
    f_new = open("./parameters/" + filename + "_temp.dat",'w')
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

            if exception(name, data): # This variable can be skipped
                valid_var = True
                var_temp = var

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

            f_new.write(name + " " + var_temp + " // " + comment + '\n')

    f_old.close()
    f_new.close()

    subprocess.call("rm " + "./parameters/" + filename + ".dat", shell=True)
    subprocess.call("mv " + "./parameters/" + filename + "_temp.dat " +
                    " ./parameters/" + filename + ".dat", shell=True)


def slurm_file_gen(data, run_time):
    """
    Generate a slurm submitting file.
    :param data: Data structure which gives parameters in the submitting script
    :param run_time: Data structure which gives estimated running time
    :return: program name
    """
    f_old = open("submit_template.run",'r')

    filename = data.model.lower() + "_" + data.size + "_" + data.task.lower()
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
            new_line = line[:prog_start + len(prog)] + filename + " > " + filename + ".out\n"
        elif dir_start > -1:
            new_line = line[:dir_start + len(dire)] + address + '\n'
        elif server_start > -1:
            new_line_write = False

        if new_line_write:
            f_new.write(new_line)

    f_old.close()
    f_new.close()

    return filename


def qsub_file_gen(data, run_time):
    """
    Generate a qsub submitting file.
    :param data: Data structure which gives parameters in the submitting script
    :param run_time: Data structure which gives estimated running time and memory
    :return: program name
    """
    f_old = open("submit_template.run",'r')

    filename = data.model.lower() + "_" + data.size + "_" + data.task.lower()
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

        new_line = line

        new_line_write = True
        if threads_start > -1:
            # Threads and template are in the same line
            new_line = line[:threads_start + len(threads)] + data.num_threads + \
                       ",cput=" + run_time.hour + ":" + run_time.min + ":00\n"
        elif prog_start > -1:
            new_line = line[:prog_start + len(prog)] + filename + " > " + filename + ".out\n"
        elif dir_start > -1:
            new_line = line[:dir_start + len(dire)] + address + '\n'
        elif server_start > -1:
            new_line_write = False

        if new_line_write:
            f_new.write(new_line)

    f_old.close()
    f_new.close()

    return filename

def exception(name, data):
    """
    Exception in asking for input according to tasks and models. If it returns true, then this
    variable will be skipped when process parameter files, and previous values are used.
    :param name: The variable name
    :param data: Data structure which holds tasks and models
    :return: Whether this variable is an exception
    """
    exc = False
    if data.task is not None:
        if (data.task == "Disorder_Transition") and (name == "J"):
            exc = True
    return exc
