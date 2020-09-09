#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  Copyright 2019 The FATE Authors. All Rights Reserved.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import csv
import json
import os
import sys
import time
import traceback

from arch.api import eggroll
from arch.api.storage import save_data

CSV = 'csv'
LOAD_DATA_COUNT = 10000
MAX_PARTITION_NUM = 32


def list_to_str(input_list):
    str1 = ''
    size = len(input_list)
    for i in range(size):
        if i == size - 1:
            str1 += str(input_list[i])
        else:
            str1 += str(input_list[i]) + ','

    return str1


def read_data(input_file='', head=True):
    split_file_name = input_file.split('.')
    if CSV in split_file_name:
        print("file type is csv")
        with open(input_file) as csv_file:
            csv_reader = csv.reader(csv_file)
            if head is True:
                csv_head = next(csv_reader)

            for row in csv_reader:
                yield (row[0], list_to_str(row[1:]))
    else:
        print("file type is not known, read it as txt")
        with open(input_file, 'r') as fin:
            if head is True:
                head = fin.readline()

            lines = fin.readlines()
            for line in lines:
                values = line.replace("\n", "").replace("\t", ",").split(",")
                yield (values[0], list_to_str(values[1:]))


def generate_table_name(input_file_path):
    local_time = time.localtime(time.time())
    str_time = time.strftime("%Y%m%d%H%M%S", time.localtime())
    file_name = input_file_path.split(".")[0]
    file_name = file_name.split("/")[-1]
    return file_name, str_time


def data_to_eggroll_table(data, namespace, table_name, partition=1, work_mode=0):
    eggroll.init(mode=work_mode)
    data_table = eggroll.table(table_name, namespace, partition=partition, create_if_missing=True, error_if_exist=False)
    data_table.put_all(data)
    data_table_count = data_table.count()
    print("------------load data finish!-----------------")
    print("total data_count:" + str(data_table.count()))
    print("namespace:%s, table_name:%s" % (namespace, table_name))
    # for kv in data_table.collect():
    #    print(kv)


def load_file(load_file_path):
    try:
        # args.config = os.path.abspath(args.config)
        input_file_path = None
        head = True
        table_name = None
        namespace = None
        with open(load_file_path, 'r') as f:
            data = json.load(f)
            try:
                input_file_path = data['file']
            except:
                traceback.print_exc()

            try:
                read_head = data['head']
                if read_head == 0:
                    head = False
                elif read_head == 1:
                    head = True
            except:
                print("'head' in .json should be 0 or 1, set head to 1")

            try:
                partition = data['partition']
                if partition <= 0 or partition > MAX_PARTITION_NUM:
                    print("Error number of partition, it should between %d and %d" % (0, MAX_PARTITION_NUM))
                    sys.exit()
            except:
                print("set partition to 1")
                partition = 1

            try:
                table_name = data['table_name']
            except:
                print("not setting table_name or setting error, set table_name according to current time")

            try:
                namespace = data['namespace']
            except:
                print("not setting namespace or setting error, set namespace according to input file name")

            work_mode = data.get('work_mode')
            if work_mode is None:
                work_mode = 0
            else:
                work_mode = int(work_mode)

        if not os.path.exists(input_file_path):
            print("%s is not exist, please check the configure" % (input_file_path))
            sys.exit()

        input_data = read_data(input_file_path, head)
        if True:
            eggroll.init(mode=work_mode)
            _namespace, _table_name = generate_table_name(input_file_path)
            if namespace is None:
                namespace = _namespace
            if table_name is None:
                table_name = _table_name
            save_data(input_data, table_name, namespace, partition, work_mode)

    except ValueError:
        print('json parse error')
        exit(-102)
    except IOError:
        print('read file error')
        exit(-103)
