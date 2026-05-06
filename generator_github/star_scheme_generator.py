import numpy as np
import tqdm
import yaml
import pandas as pd
import datetime as dt
from xeger import Xeger

class StarGenerator:
    '''
    Generator is initialised with
    '''

    def __init__(self, path = 'test_star.yaml'):
        self.read_yaml(path)
        self.dim_data = []

    def read_yaml(self, path):
        with open(path) as stream:
            raw_vars = yaml.safe_load(stream)

        self.fact_table     = raw_vars['fact table']
        self.dim_tables     = raw_vars['dimension tables']

        self.dm_activities  = raw_vars['activities']
        self.dm_followed_by = raw_vars['followed by']
        self.dm_start_time  = raw_vars['meta']['start time']

    def generate_dim_data(self):
        for table in self.dim_tables:
            for _ in range(table['nr rows']):
                record = {}
                for col in table['columns']:
                    if col['gen type'] == 'regex':
                        value = Xeger().xeger(col['gen value'])
                    elif col['gen type'] == 'choice':
                        norm = sum(col['gen probability'])
                        prob = [x / norm for x in col['gen probability']]
                        value = np.random.choice(
                            a = col['gen value'],
                            p = prob
                        )

                    for corr_c in table['columns']:
                        if      'gen type'                   in corr_c \
                                and corr_c['gen type']       == 'correlate'\
                                and corr_c['correlate with'] == col['name']:

                            for a, b in zip (col['gen value'], corr_c ['gen value']):
                                if a == value:
                                    print ('col: ', a, ' , corr: ', b)
                                    record [corr_c['name']] = b 
                            # record[corr_c['name']] = [b for a, b in
                            #                                    zip(col['gen value'], corr_c['gen value'])
                            #                                    if a == value][0]
                    record[col['name']] = value

                self.dim_data.append(
                    {'table': table['name'],
                     'data':  record
                    }
                )
