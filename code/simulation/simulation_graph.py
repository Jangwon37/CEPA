import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def count_T(data):
    count = 0
    for i in range(data.shape[0]):
        if (data[i][0] < data[i][2]) & (data[i][2] < data[i][4]):
            count = count + 1
    return (count / data.shape[0])

def precedence_P(data):
    p = []
    for i in range(data.shape[0]):
        p.append(data[i][0]*data[i][1]*data[i][2])
    return p

def graph_result(sim_name, sim_path):
    count_data = []
    for i in range(1, 101):
        count_result_data = np.array(pd.read_csv("%s/%s/output.random_data_sim_X.csv" %(sim_path, i), index_col=0))
        count_data.append(count_T(count_result_data))

    copt_result_data = np.array(pd.read_csv("%s/prAB_result.csv" %sim_path, index_col=0))
    dabrowska_result_data = np.array(pd.read_csv("%s/prAB_result_dabrowska.csv" % sim_path, index_col=0))
    linying_result_data = np.array(pd.read_csv("%s/prAB_result_linying.csv" % sim_path, index_col=0))

    copt_data = precedence_P(copt_result_data)
    dabrowska_data = precedence_P(dabrowska_result_data)
    linying_data = precedence_P(linying_result_data)

    print("%s \ncopt: %s \ndabroska: %s \nlinying: %s \ncount: %s"
        % (sim_name, round(np.mean(copt_data), 4), round(np.mean(dabrowska_data), 4), round(np.mean(linying_data), 4), round(np.mean(count_data), 4)))

    fig, ax = plt.subplots()
    ax.boxplot([copt_data, dabrowska_data, linying_data, count_data], vert=True, sym='')
    ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.8)
    plt.xticks([1, 2, 3, 4], ['CEPA', 'Dabroska', 'Linying', 'Naive'], fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel('Likelihood', fontsize=16)
    plt.savefig('%s/%s_box.tiff' % (sim_path, sim_name), dpi=300)
    plt.show()
    plt.close()


path = "D:/project/CEPA/result/simulation_data/"

data_name = ['exp', 'logn', 'uniform', 'clayton']

path_exp = "%s/X_%s_simulation/" %(path, data_name[0])
path_logn = "%s/X_%s_simulation/" %(path, data_name[1])
path_uniform = "%s/X_%s_simulation/" %(path, data_name[2])
path_clayton = "%s/X_%s_simulation/" %(path, data_name[3])

graph_result(data_name[0], path_exp)
graph_result(data_name[1], path_logn)
graph_result(data_name[2], path_uniform)
graph_result(data_name[3], path_clayton)
