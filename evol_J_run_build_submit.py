__author__ = 'liangshengzhang'

import file_process as fp
import structure as structure

def evol_J_run_file_read(tasks_models,  count, file_modify):
    """
    This function gives the function to handle reading/modifying parameters for time evolution
    project.
    :param tasks_models: The possible tasks and models
    :param count: The count number used in file readings
    :param file_modify: A map recording which parameter file is read
    :return: A data structure which contains information regarding this execution
    """
    data = fp.DATA()
    modify_words = ["size", "num_realizations", "J", "init_func_name", "time_step", "log_time"]
    # Generate parameters
    fp.para_gen("generic_para", tasks_models, data, count, modify_words)
    file_modify["generic_para"] = True

    # Generate Model data
    fp.para_gen(data.model.lower(), tasks_models, data, count, modify_words)
    file_modify[data.model.lower()] = True

    # Generate Task data
    fp.para_gen(data.task.lower(), tasks_models, data, count, modify_words)
    file_modify[data.task.lower()] = True

    print fp.bcolors.BOLD, "System size is " + str(data.size) + "\n", fp.bcolors.ENDC

    # Ask whether user want running details
    run_detail = None
    while run_detail is None:
        detail = raw_input("Show running details?")
        if (detail == '') or detail.lower().startswith('n'):
            run_detail = False
        elif detail.lower().startswith('y'):
            run_detail = True
            fp.display_file("generic_para", data, count)
            fp.display_file(data.model.lower(), data, count)
            fp.display_file(data.task.lower(), data, count)

    return data

structure.structure(evol_J_run_file_read)







