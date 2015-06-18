__author__ = 'liangshengzhang'


import file_process as fp
import structure as structure

def size_run_file_read(tasks_models,  count, file_modify):
    """
    This function gives the function to handle reading/modifying parameters for phase transition project
    :param tasks_models: The possible tasks and models
    :param count: The count number used in file readings
    :param file_modify: A map recording which parameter file is read
    :return: A data structure which contains information r
    """

    modify_words = ["size", "num_realizations"]
    data = fp.DATA()

    # Generate parameters
    fp.para_gen("generic_para", tasks_models, data, count, modify_words)

    file_modify["generic_para"] = True

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
            fp.display_file(data.model.lower(), data, count-1) # This file's name has not been updated yet
            fp.display_file(data.task.lower(), data, count-1) # This file's name has not been updated yet

    return data

structure.structure(size_run_file_read)



