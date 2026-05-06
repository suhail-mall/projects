"""
Created on 04/03/2026
Author: Suhail Mall

"""
import numpy as np
import tqdm
import yaml
import pandas as pd
import datetime as dt
from xeger import Xeger


class Generator:
    '''
    Generator is initialised with
    '''

    def __init__(self, path = 'test2.yaml'):
        self.read_yaml(path)


    def read_yaml(self, path):
        with open(path) as stream:
            raw_vars = yaml.safe_load(stream)

        self.dm_tables      = raw_vars['tables']
        self.dm_fks         = raw_vars['foreign keys']
        self.dm_activities  = raw_vars['activities']
        self.dm_followed_by = raw_vars['followed by']
        self.dm_start_time  = raw_vars['meta']['start time']

        self.all_data = [{'table': t, 'data': []} for t in self.table_names()]
        self.current_fks = []

    # Easily access high-level names
    def table_names (self):
        return [t['name'] for t in self.dm_tables]

    def activity_names (self):
        return [a['name'] for a in self.dm_activities]

    def table_columns(self, table):
        for t in self.dm_tables:
            if t['name'] == table:              
                return [c['name'] for c in t['columns']]

    def get_next_activity(self, curr_activity):
        curr_follow = [x for x in self.dm_followed_by if x['name'] == curr_activity][0]
        activities, probs = zip(*curr_follow['followed activities'])
        norm = sum(probs)
        next_activity = np.random.choice(
            a = activities,
            p = [x / norm for x in probs]
            )
        return next_activity

    def decide_new_record(self, proposed_table):

        return

    def propagate_record(self, curr_table, child_record):
        # curr table and column are the ones that we're propagating the record from
        # if child key is in parent key column in ,
        # and if parent has all parent columns populated in self.all_data,
        # then just populareh the fk in the child table.
        # fks = {parent table, parent column, value}

        for fk_reln in self.dm_fks:
            found = False
            if fk_reln['child'] == curr_table: # identifying which FK reln to look at
                # fks =
                for prev_rec in self.current_fks: # Go through the previous FKs thath ave been added for this case
                    if prev_rec['table'] == fk_reln['parent'] and prev_rec['column'] == fk_reln['parent key']: # is already in there for this case
                    # Take the value from current_fks and insert into hte curren trecord
                        found = True
                        child_record[fk_reln['child key']] = prev_rec['value'] # Populate the FK in the child table

                        # for t in self.all_data:
                        #     if t['table'] == fk_reln['parent']:
                        #         for parent_data in t['data']:
                        #             if parent_data[fk_reln['parent_key']] == prev_rec['']: # Find the record that was previously populated

                # 
                # if not found: # If this FK  hasn't been populated for this case before, generate a shared value and then propagate to the next table
                # 
                # 
                # 
                # 
                # 
                # else:


        # return

    def generate_case(self):
        # When generating a record, populate all parent tables according to specified FKs
        # For each record associated with the same case, if a previous record has joined two tables together,
        # use the same FK
        # e.g. if an order was made and then delivery

        self.current_fks = [] # {table, column, value}
        curr_time = dt.datetime.today() - dt.timedelta(
            np.random.normal(self.dm_start_time['mean'], self.dm_start_time['spread'])
            )

        finished = False
        curr_activity = self.get_next_activity('START')

        while not finished:
            curr_time = self.generate_event(curr_activity, curr_time) # function returns time of event



    def generate_event(self, curr_event, from_time):
        # Get the info for how the event will be generated as per the YAML file
        curr_event_info = [a for a in self.dm_activities if a['name'] == curr_event][0]
        new_records = []
        new_time = from_time + dt.timedelta(
                 hours=np.random.normal(curr_event_info['time mean'], curr_event_info['time spread']))
        new_records.append(
            {'table': curr_event_info['timestamp']['table'],
             'data': {curr_event_info['timestamp']['column']: new_time}
             }
        )

        for origin in curr_event_info['origins']:
            found = False
            for record in new_records: # modify in place
                if origin['table'] == record['table']:
                    for column in origin['columns']:
                        if 'value' in column:
                            record['data'][column['name']] = column['value']
                        else:
                            record['data'] = self.generate_record(
                                origin['table'],
                                column['name'],
                                record['data']
                                )
                        found = True
            if not found:
                record = {}
                for column in origin['columns']:
                    record = self.generate_record(
                        origin['table'],
                        column['name'],
                        record
                        )
                new_records.append({
                    'table': origin['table'],
                    'data':  record
                    })


        # [t['data'].append(record['data'])
        #     for record in new_records
        #     for t in self.all_data
        #     if record['table'] == t['table']]

        # record is {table, data}
        # For record in new_records, self.propagate_record -- PUT IN THE FOR LOOP ABOVE ??
        for record in new_records:
            for t in self.all_data:
                if record['table'] == t['table']:
                    t['data'].append(record['data'])


        return new_time



    def generate_record(self, table, column, existing_record = {}):
        for t in self.dm_tables:
            if t['name'] == table:
                for c in t['columns']:
                    if c['name'] == column:
                        found_t = t
                        found_c = c
        # c has name, gen_type, gen_val, unique
        if found_c['gen type'] == 'regex':
            value = Xeger().xeger(found_c['gen value'])
        elif found_c['gen type'] == 'choice':
            norm = sum(found_c['gen probability']) # Do this once for efficiency
            prob = [x / norm for x in found_c['gen probability']] # normalise probabilities
            value = np.random.choice(
                a = found_c['gen value'],
                p = prob
            )
            for corr_c in found_t['columns']: # If there are any correlating values generated in the same table
                if 'gen type' in corr_c and corr_c['gen type'] == 'correlate': # in case there is no gen type for that record
                    existing_record[corr_c['name']] = [b for a,b in
                                                       zip(found_c['gen value'], corr_c['gen value'])
                                                       if a == value][0]

        existing_record[column] = value
        # Something like this - THEN MOVE TO FUNTION ABOVE to increase efficiency
        # print(found_t['name'])
        for col in found_t['columns']:
            # print(col['name'])
            if col['name'] not in existing_record:
                existing_record[col['name']] = None
        return existing_record


