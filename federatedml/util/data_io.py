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
#
################################################################################
#
#
################################################################################

import functools
import numpy as np
from arch.api.utils import log_utils
from federatedml.feature.instance import Instance
from federatedml.feature.sparse_vector import SparseVector
from federatedml.util import consts
from federatedml.util.param_checker import DataIOParamChecker
from federatedml.util import abnormal_detection
from federatedml.statistic import data_overview
from arch.api.model_manager import manager
from arch.api.proto.data_io_meta_pb2 import DataIOMeta
from arch.api.proto.data_io_param_pb2 import DataIOParam
from arch.api.proto.feature_imputer_meta_pb2 import ImputerMeta
from arch.api.proto.feature_imputer_param_pb2 import ImputerParam
from arch.api.proto.feature_outlier_meta_pb2 import OutlierMeta
from arch.api.proto.feature_outlier_param_pb2 import OutlierParam
from arch.api import storage

LOGGER = log_utils.getLogger()


# =============================================================================
# DenseFeatureReader
# =============================================================================
class DenseFeatureReader(object):
    def __init__(self, data_io_param):
        DataIOParamChecker.check_param(data_io_param)
        self.delimitor = data_io_param.delimitor
        self.data_type = data_io_param.data_type
        self.missing_fill = data_io_param.missing_fill
        self.default_value = data_io_param.default_value
        self.missing_fill_method = data_io_param.missing_fill_method
        self.missing_impute = data_io_param.missing_impute
        self.outlier_replace = data_io_param.outlier_replace
        self.outlier_replace_method = data_io_param.outlier_replace_method
        self.outlier_impute = data_io_param.outlier_impute
        self.outlier_replace_value = data_io_param.outlier_replace_value
        self.with_label = data_io_param.with_label
        self.label_idx = data_io_param.label_idx
        self.label_type = data_io_param.label_type
        self.output_format = data_io_param.output_format
        self.header = None

    def generate_header(self, input_data_feature, table_name, namespace):
        self.header = storage.get_data_table_meta("header", table_name, namespace)

        if not self.header:
            feature_shape = data_overview.get_data_shape(input_data_feature)
            self.header = ["fid" + str(i) for i in range(feature_shape)]
        else:
            if self.with_label:
                self.header = self.header.split(self.delimitor, -1)[: self.label_idx] + \
                              self.header.split(self.delimitor, -1)[self.label_idx + 1:]
            else:
                self.header = self.header.split(self.delimitor, -1)

    def read_data(self, table_name, namespace, mode="fit"):
        input_data = storage.get_data_table(table_name, namespace)
        LOGGER.info("start to read dense data and change data to instance")
        
        abnormal_detection.empty_table_detection(input_data)
        
        input_data_features = None
        input_data_labels = None

        if self.with_label:
            if type(self.label_idx).__name__ != "int":
                raise ValueError("label index should be integer")

            data_shape = data_overview.get_data_shape(input_data)
            if not data_shape or self.label_idx >= data_shape:
                raise ValueError("input data's value is empty, it does not contain a label")

            input_data_features = input_data.mapValues(
                lambda value: [] if data_shape == 1 else value.split(self.delimitor, -1)[:self.label_idx] + value.split(self.delimitor, -1)[
                                                                                 self.label_idx + 1:])
            input_data_labels = input_data.mapValues(lambda value: value.split(self.delimitor, -1)[self.label_idx])

        else:
            input_data_features = input_data.mapValues(lambda value: [] if not value else value.split(self.delimitor, -1))

        if mode == "fit":
            data_instance = self.fit(input_data_features, input_data_labels, table_name, namespace)
        else:
            data_instance = self.transform(input_data_features, input_data_labels)

        set_schema(data_instance, self.header)

        return data_instance

    def fit(self, input_data_features, input_data_labels, table_name, namespace):
        input_data_features = self.fill_missing_value(input_data_features, "fit")
        input_data_features = self.replace_outlier_value(input_data_features, "fit")

        self.generate_header(input_data_features, table_name, namespace)

        data_instance = self.gen_data_instance(input_data_features, input_data_labels)

        return data_instance

    def transform(self, input_data_features, input_data_labels):
        input_data_features = self.fill_missing_value(input_data_features, "transform")
        input_data_features = self.replace_outlier_value(input_data_features, "transform")

        data_instance = self.gen_data_instance(input_data_features, input_data_labels)
        return data_instance

    def fill_missing_value(self, input_data_features, mode="fit"):
        if self.missing_fill:
            from federatedml.feature.imputer import Imputer
            imputer_processor = Imputer(self.missing_impute)
            if mode == "fit":
                input_data_features, self.default_value = imputer_processor.fit(input_data_features,
                                                                                replace_method=self.missing_fill_method,
                                                                                replace_value=self.default_value)
                if self.missing_impute is None:
                    self.missing_impute = imputer_processor.get_imputer_value_list()
            else:
                input_data_features = imputer_processor.transform(input_data_features,
                                                                  replace_method=self.missing_fill_method,
                                                                  transform_value=self.default_value)

            if self.missing_impute is None:
                self.missing_impute = imputer_processor.get_imputer_value_list()

        return input_data_features

    def replace_outlier_value(self, input_data_features, mode="fit"):
        if self.outlier_replace:
            from federatedml.feature.imputer import Imputer
            imputer_processor = Imputer(self.outlier_impute)
            if mode == "fit":
                input_data_features, self.outlier_replace_value = \
                    imputer_processor.fit(input_data_features,
                                          replace_method=self.outlier_replace_method,
                                          replace_value=self.outlier_replace_value)

                if self.outlier_impute is None:
                    self.outlier_impute = imputer_processor.get_imputer_value_list()
            else:
                input_data_features = imputer_processor.transform(input_data_features,
                                                                  replace_method=self.outlier_replace_method,
                                                                  transform_value=self.outlier_replace_value)

        return input_data_features

    def gen_data_instance(self, input_data_features, input_data_labels):
        if self.with_label:
            data_instance = input_data_features.join(input_data_labels,
                                                     lambda features, label:
                                                     self.to_instance(features, label))
        else:
            data_instance = input_data_features.mapValues(lambda features: self.to_instance(features))

        return data_instance

    def to_instance(self, features, label=None):
        if self.with_label:
            if self.label_type == 'int':
                label = int(label)
            elif self.label_type in ["float", "float64"]:
                label = float(label)

            features = DenseFeatureReader.gen_output_format(features, self.data_type, self.output_format,
                                                            missing_impute=self.missing_impute)

        else:
            features = DenseFeatureReader.gen_output_format(features, self.data_type, self.output_format,
                                                            missing_impute=self.missing_impute)

        return Instance(inst_id=None,
                        features=features,
                        label=label)

    @staticmethod
    def gen_output_format(features, data_type='float', output_format='dense', missing_impute=None):

        if output_format not in ["dense", "sparse"]:
            raise ValueError("output format {} is not define".format(output_format))

        if output_format == "dense":
            return np.asarray(features, dtype=data_type)

        indices = []
        data = []
        column_shape = len(features)
        non_zero = 0

        for i in range(column_shape):
            if (missing_impute is not None and features[i] in missing_impute) or \
                    (missing_impute is None and features[i] in ['', 'NULL', 'null', "NA"]):
                continue

            if data_type in ['float', 'float64']:
                if np.fabs(float(features[i])) < consts.FLOAT_ZERO:
                    continue

                indices.append(i)
                data.append(float(features[i]))
                non_zero += 1

            elif data_type in ['int']:
                if int(features[i]) == 0:
                    continue
                indices.append(i)
                data.append(int(features[i]))

            else:
                indices.append(i)
                data.append(features[i])

        return SparseVector(indices, data, column_shape)

    def save_model(self, model_table, model_namespace):
        model_types = []

        if model_table is None or model_namespace is None:
            LOGGER.info("data io model-meta can't save, it will not reuse in transform stage")
            return []

        save_data_io_model(input_format="dense",
                           delimitor=self.delimitor,
                           data_type=self.data_type,
                           with_label=self.with_label,
                           label_idx=self.label_idx,
                           label_type=self.label_type,
                           output_format=self.output_format,
                           header=self.header,
                           model_name="DenseFeatureReader",
                           model_table=model_table,
                           model_namespace=model_namespace)
        model_types.append(("DenseFeatureReader.meta", "DenseFeatureReader.param"))

        save_missing_imputer_model(self.missing_fill,
                                   self.missing_fill_method,
                                   self.missing_impute,
                                   self.default_value,
                                   self.header,
                                   "Imputer",
                                   model_table,
                                   model_namespace)
        model_types.append(("Imputer.meta", "Imputer.param"))

        save_outlier_model(self.outlier_replace,
                           self.outlier_replace_method,
                           self.outlier_impute,
                           self.outlier_replace_value,
                           self.header,
                           "Outlier",
                           model_table,
                           model_namespace)
        model_types.append(("Outlier.meta", "Outlier.param"))

        return model_types

    def load_model(self, model_table, model_namespace):
        if model_table is None or model_namespace is None:
            LOGGER.info("data io model-meta can't reuse, model table name or namespace is null")
            return
        
        self.delimitor, self.data_type, _1, _2, self.with_label, \
        self.label_idx, self.label_type, self.output_format, self.header = load_data_io_model("DenseFeatureReader",
                                                                                              model_table,
                                                                                              model_namespace)

        self.missing_fill, self.missing_fill_method, \
        self.missing_impute, self.default_value = load_missing_imputer_model(self.header,
                                                                             "Imputer",
                                                                             model_table,
                                                                             model_namespace)

        self.outlier_replace, self.outlier_replace_method, \
        self.outlier_impute, self.outlier_replace_value = load_outlier_model(self.header,
                                                                             "Outlier",
                                                                             model_table,
                                                                             model_namespace)


# =============================================================================
# SparseFeatureReader: mainly for libsvm input format
# =============================================================================
class SparseFeatureReader(object):
    def __init__(self, data_io_param):
        DataIOParamChecker.check_param(data_io_param)
        self.delimitor = data_io_param.delimitor
        self.data_type = data_io_param.data_type
        self.label_type = data_io_param.label_type
        self.output_format = data_io_param.output_format
        self.header = None

    def get_max_feature_index(self, line, delimitor=' '):
        if line.strip() == '':
            raise ValueError("find an empty line, please check!!!")

        cols = line.split(delimitor, -1)
        if len(cols) <= 1:
            return -1

        return max([int(fid_value.split(":", -1)[0]) for fid_value in cols[1:]])

    def generate_header(self, max_feature):
        self.header = [str(i) for i in range(max_feature + 1)]

    def read_data(self, table_name, namespace, mode="fit"):
        input_data = storage.get_data_table(table_name, namespace)
        LOGGER.info("start to read sparse data and change data to instance")

        abnormal_detection.empty_table_detection(input_data)
        
        if not data_overview.get_data_shape(input_data):
            raise ValueError("input data's value is empty, it does not contain a label")
        
        if mode == "fit":
            data_instance = self.fit(input_data)
        else:
            data_instance = self.transform(input_data)

        set_schema(data_instance, self.header)
        return data_instance

    def fit(self, input_data):
        get_max_fid = functools.partial(self.get_max_feature_index, delimitor=self.delimitor)
        max_feature = input_data.mapValues(get_max_fid).reduce(lambda max_fid1, max_fid2: max(max_fid1, max_fid2))

        if max_feature == -1:
            raise ValueError("no feature value in input data, please check!")

        self.generate_header(max_feature)

        data_instance = self.gen_data_instance(input_data, max_feature)
        return data_instance

    def transform(self, input_data):
        max_feature = len(self.header)

        data_instance = self.gen_data_instance(input_data, max_feature)
        return data_instance

    def gen_data_instance(self, input_data, max_feature):
        params = [self.delimitor, self.data_type,
                  self.label_type,
                  self.output_format, max_feature]

        to_instance_with_param = functools.partial(self.to_instance, params)
        data_instance = input_data.mapValues(to_instance_with_param)

        return data_instance

    @staticmethod
    def to_instance(param_list, value):
        delimitor = param_list[0]
        data_type = param_list[1]
        label_type = param_list[2]
        output_format = param_list[3]
        max_fid = param_list[4]

        if output_format not in ["dense", "sparse"]:
            raise ValueError("output format {} is not define".format(output_format))

        cols = value.split(delimitor, -1)

        label = cols[0]
        if label_type == 'int':
            label = int(label)
        elif label_type in ["float", "float64"]:
            label = float(label)

        fid_value = []
        for i in range(1, len(cols)):
            fid, val = cols[i].split(":", -1)

            fid = int(fid)
            if data_type in ["float", "float64"]:
                val = float(val)
            elif data_type in ["int", "int64"]:
                val = int(val)

            fid_value.append((fid, val))

        if output_format == "dense":
            features = [0 for i in range(max_fid + 1)]
            for fid, val in fid_value:
                features[fid] = val

            features = np.asarray(features, dtype=data_type)

        else:
            indices = []
            data = []
            for fid, val in fid_value:
                indices.append(fid)
                data.append(val)

            features = SparseVector(indices, data, max_fid + 1)

        return Instance(inst_id=None,
                        features=features,
                        label=label)

    def save_model(self, model_table, model_namespace):
        model_types = []

        save_data_io_model(input_format="sparse",
                           delimitor=self.delimitor,
                           data_type=self.data_type,
                           label_type=self.label_type,
                           output_format=self.output_format,
                           header=self.header,
                           model_name="SparseFeatureReader",
                           model_table=model_table,
                           model_namespace=model_namespace)
        model_types.append(("SparseFeatureReader.meta", "SparseFeatureReader.param"))

        save_missing_imputer_model(missing_fill=False,
                                   model_name="Imputer",
                                   model_table=model_table,
                                   model_namespace=model_namespace)
        model_types.append(("Imputer.meta", "Imputer.param"))

        save_outlier_model(outlier_replace=False,
                           model_name="Outlier",
                           model_table=model_table,
                           model_namespace=model_namespace)
        model_types.append(("Outlier.meta", "Outlier.param"))

        return model_types

    def load_model(self, model_table, model_namespace):
        self.delimitor, self.data_type, _1, _2, _3, _4, \
        self.label_type, self.output_format, self.header = load_data_io_model("SparseFeatureReader",
                                                                                  model_table,
                                                                                  model_namespace)


# =============================================================================
# SparseTagReader: mainly for tag data
# =============================================================================
class SparseTagReader(object):
    def __init__(self, data_io_param):
        DataIOParamChecker.check_param(data_io_param)
        self.delimitor = data_io_param.delimitor
        self.data_type = data_io_param.data_type
        self.tag_with_value = data_io_param.tag_with_value
        self.tag_value_delimitor = data_io_param.tag_value_delimitor
        self.with_label = data_io_param.with_label
        self.label_type = data_io_param.label_type
        self.output_format = data_io_param.output_format
        self.header = None

    @staticmethod
    def agg_tag(kvs, delimitor=' ', with_label=True, tag_with_value=False, tag_value_delimitor=":"):
        tags_set = set()
        for key, value in kvs:
            if with_label:
                cols = value.split(delimitor, -1)[1 : ]
            else:
                cols = value.split(delimitor, -1)[0:]

            if tag_with_value is False:
                tags = cols
            else:
                tags = [fea_value.split(tag_value_delimitor, -1)[0] for fea_value in cols]

            tags_set |= set(tags)

        return tags_set

    def generate_header(self, tags):
        self.header = tags

    def read_data(self, table_name, namespace, mode="fit"):
        input_data = storage.get_data_table(table_name, namespace)
        LOGGER.info("start to read sparse data and change data to instance")

        abnormal_detection.empty_table_detection(input_data)
        
        if mode == "fit":
            data_instance = self.fit(input_data)
        else:
            data_instance = self.transform(input_data)

        set_schema(data_instance, self.header)
        return data_instance

    def fit(self, input_data):
        tag_aggregator = functools.partial(SparseTagReader.agg_tag, 
                                           delimitor=self.delimitor,
                                           with_label=self.with_label,
                                           tag_with_value=self.tag_with_value, 
                                           tag_value_delimitor=self.tag_value_delimitor)
        tags_set_list = list(input_data.mapPartitions(tag_aggregator).collect())
        tags_set = set()
        for _, _tags_set in tags_set_list:
            tags_set |= _tags_set
        tags = list(tags_set)

        tags = sorted(tags)
        tags_dict = dict(zip(tags, range(len(tags))))

        self.generate_header(tags)

        data_instance = self.gen_data_instance(input_data, tags_dict)
        return data_instance

    def transform(self, input_data):
        tags_dict = dict(zip(self.header, range(len(self.header))))

        data_instance = self.gen_data_instance(input_data, tags_dict)
        return data_instance

    def gen_data_instance(self, input_data, tags_dict):
        params = [self.delimitor,
                  self.data_type,
                  self.tag_with_value,
                  self.tag_value_delimitor,
                  self.with_label,
                  self.label_type,
                  self.output_format,
                  tags_dict]

        to_instance_with_param = functools.partial(self.to_instance, params)
        data_instance = input_data.mapValues(to_instance_with_param)

        return data_instance

    @staticmethod
    def to_instance(param_list, value):
        delimitor = param_list[0]
        data_type = param_list[1]
        tag_with_value = param_list[2]
        tag_value_delimitor = param_list[3]
        with_label = param_list[4]
        label_type = param_list[5]
        output_format = param_list[6]
        tags_dict = param_list[7]

        if output_format not in ["dense", "sparse"]:
            raise ValueError("output format {} is not define".format(output_format))

        cols = value.split(delimitor, -1)
        start_pos = 0
        label = None

        if with_label:
            start_pos = 1
            label = cols[0]
            if label_type == 'int':
                label = int(label)
            elif label_type in ["float", "float64"]:
                label = float(label)

        if output_format == "dense":
            features = [0 for i in range(len(tags_dict))]
            for fea in cols[start_pos:]:
                if tag_with_value:
                    _tag, _val = fea.split(tag_value_delimitor, -1)
                    features[tags_dict.get(_tag)] = _val
                else:
                    features[tags_dict.get(fea)] = 1

            features = np.asarray(features, dtype=data_type)
        else:
            indices = []
            data = []
            for fea in cols[start_pos:]:
                if tag_with_value:
                    _tag, _val = fea.split(tag_value_delimitor, -1)
                else:
                    _tag = fea
                    _val = 1
                indices.append(tags_dict.get(_tag))
                if data_type in ["float", "float64"]:
                    _val = float(_val)
                elif data_type in ["int", "int64", "long"]:
                    _val = int(_val)
                elif data_type == "str":
                    _val = str(_val)

                data.append(_val)

            features = SparseVector(indices, data, len(tags_dict))

        return Instance(inst_id=None,
                        features=features,
                        label=label)

    def save_model(self, model_table, model_namespace):
        model_types = []

        save_data_io_model(input_format="tag",
                           delimitor=self.delimitor,
                           data_type=self.data_type,
                           tag_with_value=self.tag_with_value,
                           tag_value_delimitor=self.tag_value_delimitor,
                           with_label=self.with_label,
                           label_type=self.label_type,
                           output_format=self.output_format,
                           header=self.header,
                           model_name="SparseTagReader",
                           model_table=model_table,
                           model_namespace=model_namespace)
        model_types.append(("SparseTagReader.meta", "SparseTagReader.param"))

        save_missing_imputer_model(missing_fill=False,
                                   model_name="Imputer",
                                   model_table=model_table,
                                   model_namespace=model_namespace)
        model_types.append(("Imputer.meta", "Imputer.param"))

        save_outlier_model(outlier_replace=False,
                           model_name="Outlier",
                           model_table=model_table,
                           model_namespace=model_namespace)
        model_types.append(("Outlier.meta", "Outlier.param"))

        return model_types

    def load_model(self, model_table, model_namespace):
        self.delimitor, self.data_type, self.tag_with_value, self.tag_value_delimitor, self.with_label, \
        _1, self.label_type, self.output_format, self.header = load_data_io_model("SparseTagReader",
                                                                                  model_table,
                                                                                  model_namespace)


def set_schema(data_instance, header):
    data_instance.schema = {"header": header}


def save_data_io_model(input_format="dense",
                       delimitor=",",
                       data_type="str",
                       tag_with_value=False,
                       tag_value_delimitor=":",
                       with_label=False,
                       label_idx=0,
                       label_type="int",
                       output_format="dense",
                       header=None,
                       model_name="DataIO",
                       model_table=None,
                       model_namespace=None):
    model_meta = DataIOMeta()
    model_param = DataIOParam()

    model_meta.input_format = input_format
    model_meta.delimitor = delimitor
    model_meta.data_type = data_type
    model_meta.tag_with_value = tag_with_value
    model_meta.tag_value_delimitor = tag_value_delimitor
    model_meta.with_label = with_label
    model_meta.label_idx = label_idx
    model_meta.label_type = label_type
    model_meta.output_format = output_format

    if header is not None:
        model_param.header.extend(header)

    manager.save_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.save_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)


def load_data_io_model(model_name="DataIO",
                       model_table=None,
                       model_namespace=None):
    model_meta = DataIOMeta()
    model_param = DataIOParam()

    manager.read_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.read_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)

    delimitor = model_meta.delimitor
    data_type = model_meta.data_type
    tag_with_value = model_meta.tag_with_value
    tag_value_delimitor = model_meta.tag_value_delimitor
    with_label = model_meta.with_label
    label_idx = model_meta.label_idx
    label_type = model_meta.label_type
    output_format = model_meta.output_format

    header = list(model_param.header)

    return delimitor, data_type, tag_with_value, tag_value_delimitor, with_label, label_idx, label_type, output_format, header


def save_missing_imputer_model(missing_fill=False,
                               missing_replace_method=None,
                               missing_impute=None,
                               missing_fill_value=None,
                               header=None,
                               model_name="Imputer",
                               model_table=None,
                               model_namespace=None):
    model_meta = ImputerMeta()
    model_param = ImputerParam()

    model_meta.is_imputer = missing_fill
    if missing_fill:
        if missing_replace_method:
            model_meta.strategy = str(missing_replace_method)

        if missing_impute is not None:
            if missing_impute is not None:
                model_meta.missing_value.extend(map(str, missing_impute))

        if missing_fill_value is not None:
            feature_value_dict = dict(zip(header, map(str, missing_fill_value)))

            model_param.missing_replace_value.update(feature_value_dict)

    manager.save_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.save_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)


def load_missing_imputer_model(header=None,
                               model_name="Imputer",
                               model_table=None,
                               model_namespace=None):
    model_meta = ImputerMeta()
    model_param = ImputerParam()
    manager.read_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.read_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)

    missing_fill = model_meta.is_imputer
    missing_replace_method = model_meta.strategy
    missing_value = model_meta.missing_value
    missing_fill_value = model_param.missing_replace_value

    if missing_fill:
        if not missing_replace_method:
            missing_replace_method = None

        if not missing_value:
            missing_value = None
        else:
            missing_value = list(missing_value)

        if missing_fill_value:
            missing_fill_value = [missing_fill_value.get(head) for head in header]
        else:
            missing_fill_value = None
    else:
        missing_replace_method = None
        missing_value = None
        missing_fill_value = None

    return missing_fill, missing_replace_method, missing_value, missing_fill_value


def save_outlier_model(outlier_replace=False,
                       outlier_replace_method=None,
                       outlier_impute=None,
                       outlier_replace_value=None,
                       header=None,
                       model_name="Outlier",
                       model_table=None,
                       model_namespace=None):
    model_meta = OutlierMeta()
    model_param = OutlierParam()

    model_meta.is_outlier = outlier_replace
    if outlier_replace:
        if outlier_replace_method:
            model_meta.strategy = str(outlier_replace_method)

        if outlier_impute:
            model_meta.outlier_value.extend(map(str, outlier_impute))

        if outlier_replace_value:
            outlier_value_dict = dict(zip(header, map(str, outlier_replace_value)))
            model_param.outlier_replace_value.update(outlier_value_dict)

    manager.save_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.save_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)


def load_outlier_model(header=None,
                       model_name="Outlier",
                       model_table=None,
                       model_namespace=None):
    model_meta = OutlierMeta()
    model_param = OutlierParam()

    manager.read_model(buffer_type=model_name + ".meta",
                       proto_buffer=model_meta,
                       name=model_table,
                       namespace=model_namespace)

    manager.read_model(buffer_type=model_name + ".param",
                       proto_buffer=model_param,
                       name=model_table,
                       namespace=model_namespace)

    outlier_replace = model_meta.is_outlier
    outlier_replace_method = model_meta.strategy
    outlier_value = model_meta.outlier_value
    outlier_replace_value = model_param.outlier_replace_value

    if outlier_replace:
        if not outlier_replace_method:
            outlier_replace_method = None

        if not outlier_value:
            outlier_value = None
        else:
            outlier_value = list(outlier_value)
        

        if outlier_replace_value:
            outlier_replace_value = [outlier_replace_value.get(head) for head in header]
        else:
            outlier_replace_value = None
    else:
        outlier_replace_method = None
        outlier_value = None
        outlier_replace_value = None
    
    return outlier_replace, outlier_replace_method, outlier_value, outlier_replace_value

