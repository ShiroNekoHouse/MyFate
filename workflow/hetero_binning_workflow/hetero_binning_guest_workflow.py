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

from arch.api.utils import log_utils
from federatedml.feature.hetero_feature_binning.hetero_binning_guest import HeteroFeatureBinningGuest
from federatedml.param import FeatureBinningParam
from federatedml.util import FeatureBinningParamChecker
from federatedml.util import ParamExtract
from federatedml.util import consts
from workflow import status_tracer_decorator
from workflow.workflow import WorkFlow

LOGGER = log_utils.getLogger()


class HeteroBinningGuestWorkflow(WorkFlow):
    def _initialize(self, config_path):
        self._initialize_role_and_mode()
        self._initialize_model(config_path)
        self._initialize_workflow_param(config_path)

    def _initialize_role_and_mode(self):
        self.role = consts.GUEST
        self.mode = consts.HETERO

    def _initialize_intersect(self, config):
        pass

    def _initialize_model(self, runtime_conf_path):
        binning_param = FeatureBinningParam()
        self.binning_param = ParamExtract.parse_param_from_config(binning_param, runtime_conf_path)
        FeatureBinningParamChecker.check_param(self.binning_param)
        self.model = HeteroFeatureBinningGuest(self.binning_param)
        LOGGER.debug("Guest model started")

    def save_binning_result(self):

        meta_table = self.model.save_model(self.workflow_param.model_table, self.workflow_param.model_namespace)
        return meta_table

    @status_tracer_decorator.status_trace
    def run(self):
        self._init_argument()

        if self.workflow_param.method == "binning":

            if self.binning_param.process_method == 'fit':
                train_data_instance = self.gen_data_instance(self.workflow_param.train_input_table,
                                                             self.workflow_param.train_input_namespace,
                                                             mode='fit')
                LOGGER.debug("After dataio, header is : {}".format(train_data_instance.schema))
                if self.binning_param.local_only:
                    self.model.fit_local(train_data_instance)
                else:
                    LOGGER.debug("Start model fit")
                    self.model.fit(train_data_instance)
                self.model.save_model(self.workflow_param.model_table, self.workflow_param.model_namespace)
                train_data_instance = self.one_hot_encoder_fit_transform(train_data_instance)

            else:
                train_data_instance = self.gen_data_instance(self.workflow_param.train_input_table,
                                                             self.workflow_param.train_input_namespace,
                                                             mode='transform')
                LOGGER.debug("After dataio, header is : {}".format(train_data_instance.schema))
                self.model.load_model(self.workflow_param.model_table, self.workflow_param.model_namespace)
                if self.binning_param.local_only:
                    self.model.transform_local(train_data_instance)
                else:
                    self.model.transform(train_data_instance)
                self.save_binning_result()
                train_data_instance = self.one_hot_encoder_transform(train_data_instance)
            self._show_data(train_data_instance)
        else:
            raise TypeError("method %s is not support yet" % (self.workflow_param.method))

        LOGGER.info("Task end")

    def _show_data(self, data_instances):
        local_data = data_instances.collect()
        LOGGER.debug("data header: {}".format(data_instances.schema))
        n = 0
        for k, v in local_data:
            LOGGER.debug("new data is :{}".format(v.features))
            n += 1
            if n >= 20:
                break


if __name__ == "__main__":
    workflow = HeteroBinningGuestWorkflow()
    workflow.run()
