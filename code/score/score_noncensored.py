from itertools import *
import pandas as pd
import numpy as np

def calculate_probability(events, preceding_probabilities):
    total_permutations = list(permutations(events))
    total_cases = len(total_permutations)

    event_probabilities = []
    for case in total_permutations:
        event_prob = 1.0

        combinations_list = list(combinations(case, 2))
        for combination in combinations_list:
            event_prob *= preceding_probabilities[combination]
        event_probabilities.append(event_prob)
        # for event_pair in zip(case[:-1], case[1:]):
        #     event_prob *= preceding_probabilities[event_pair]
        # event_probabilities.append(event_prob)

    probability_dict = {}
    for case, prob in zip(total_permutations, event_probabilities):
        probability_dict[case] = prob / total_cases

    return probability_dict


events = ['INFE', 'EYE', 'RESP', 'DIGE', 'SKIN', 'INJU', 'EAR']
preceding_probabilities = {('INFE', 'EYE'): 0.7, ('INFE', 'RESP'): 0.12, ('INFE', 'DIGE'): 0.62, ('INFE', 'SKIN'): 0.48, ('INFE', 'INJU'): 0.81, ('INFE', 'EAR'): 0.67,
                           ('EYE', 'INFE'): 0.3, ('EYE', 'RESP'): 0.08, ('EYE', 'DIGE'): 0.42, ('EYE', 'SKIN'): 0.29, ('EYE', 'INJU'): 0.59, ('EYE', 'EAR'): 0.43,
                           ('RESP', 'INFE'): 0.88, ('RESP', 'EYE'): 0.92, ('RESP', 'DIGE'): 0.88, ('RESP', 'SKIN'): 0.79, ('RESP', 'INJU'): 0.99, ('RESP', 'EAR'): 0.95,
                           ('DIGE', 'INFE'): 0.38, ('DIGE', 'EYE'): 0.58, ('DIGE', 'RESP'): 0.12, ('DIGE', 'SKIN'): 0.37, ('DIGE', 'INJU'): 0.66, ('DIGE', 'EAR'): 0.52,
                           ('SKIN', 'INFE'): 0.52, ('SKIN', 'EYE'): 0.71, ('SKIN', 'RESP'): 0.21, ('SKIN', 'DIGE'): 0.63, ('SKIN', 'INJU'): 0.80, ('SKIN', 'EAR'): 0.67,
                           ('INJU', 'INFE'): 0.19, ('INJU', 'EYE'): 0.41, ('INJU', 'RESP'): 0.01, ('INJU', 'DIGE'): 0.34, ('INJU', 'SKIN'): 0.20, ('INJU', 'EAR'): 0.32,
                           ('EAR', 'INFE'): 0.33, ('EAR', 'EYE'): 0.57, ('EAR', 'RESP'): 0.05, ('EAR', 'DIGE'): 0.48, ('EAR', 'SKIN'): 0.33, ('EAR', 'INJU'): 0.68}

probability_dict = calculate_probability(events, preceding_probabilities)

pro_list = []
pro = []
for case, probability in probability_dict.items():
    pro_list.append(case)
    pro.append(probability)
    print(f'event case: {case}, prob: {probability}')

pro_list_df = pd.DataFrame(pro_list)
pro_df = pd.DataFrame(-np.log10(np.array(pro)))
result = pd.concat([pro_list_df, pro_df], axis=1)
result.to_csv('../../data/noncensored/nonensored_permutation.csv', header=None)