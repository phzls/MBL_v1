__author__ = 'liangshengzhang'

import task_model as tm
import file_process as fp
import structure

def build_submit_file_read(tasks_models,  count, file_modify):
    """
    This function gives the general function to handle reading/modifying parameters
    :param tasks_models: The possible tasks and models
    :param count: The count number used in file readings
    :param file_modify: A map recording which parameter file is read
    :return: A data structure which contains information regarding this execution
    """

    data = fp.DATA()

    # Generate generic parameters
    fp.para_gen("generic_para", tasks_models, data, count)

    file_modify["generic_para"] = True

    """
    print "Tasks:"
    for n in tasks_models.tasks:
        print n
    print "Models:"
    for n in tasks_models.models:
        print n
        m = n.lower()
        model_name = ' '.join(n.split())
        f = open("./parameters/" + m+".dat",'w')
        f.write("//" + model_name + "Parameters\n")
        f.write("//\n")
        f.write("J 0.6 // Disorder Parameter\n")
        f.close()
    """

    print fp.bcolors.BOLD, "System size is " + str(data.size) + "\n", fp.bcolors.ENDC

    # Ask whether need to change output data
    output_a = False
    while output_a is False:
        answer = raw_input("Modify output parameters?\n")
        if (answer == "Yes") or (answer == "yes") or (answer == "Y") or (answer == "y"):
            fp.para_gen("output_para", tasks_models, data, count)
            output_a = True
            file_modify["output_para"] = True
        elif (answer == 'No') or (answer == 'no') or (answer == 'N') or (answer == 'n'):
            output_a = True
        else:
            print fp.bcolors.FAIL + "Answer must be Yes(Y) or No(N)." + fp.bcolors.ENDC

    # Generate Model data
    model_filename = fp.Model_Name_Gen(data.model)
    fp.para_gen(model_filename, tasks_models, data, count)
    file_modify[model_filename] = True

    # Generate Task data
    fp.para_gen(data.task.lower(), tasks_models, data, count)
    file_modify[data.task.lower()] = True

    return data

structure.structure(build_submit_file_read)


