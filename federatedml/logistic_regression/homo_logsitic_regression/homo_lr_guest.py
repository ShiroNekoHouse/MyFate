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

import functools

import numpy as np

from arch.api import federation
from arch.api.utils import log_utils
from federatedml.evaluation import Evaluation
from federatedml.logistic_regression.base_logistic_regression import BaseLogisticRegression
from federatedml.model_selection import MiniBatch
from federatedml.optim import Initializer
from federatedml.optim import Optimizer
from federatedml.optim import activation
from federatedml.optim.federated_aggregator.homo_federated_aggregator import HomoFederatedAggregator
from federatedml.optim.gradient import LogisticGradient
from federatedml.param import LogisticParam
from federatedml.util import consts
from federatedml.util.transfer_variable import HomoLRTransferVariable
from federatedml.statistic import data_overview

LOGGER = log_utils.getLogger()


class HomoLRGuest(BaseLogisticRegression):
    def __init__(self, params: LogisticParam):
        super(HomoLRGuest, self).__init__(params)
        self.learning_rate = params.learning_rate
        self.aggregator = HomoFederatedAggregator
        self.gradient_operator = LogisticGradient()
        self.party_weight = params.party_weight

        self.optimizer = Optimizer(learning_rate=self.learning_rate, opt_method_name=params.optimizer)
        self.transfer_variable = HomoLRTransferVariable()
        self.initializer = Initializer()
        self.classes_ = [0, 1]

        self.evaluator = Evaluation()
        self.header = []
        self.penalty = params.penalty
        self.loss_history = []
        self.is_converged = False

    def fit(self, data_instances):
        self._abnormal_detection(data_instances)

        self.header = data_instances.schema.get('header')  # ['x1', 'x2', 'x3' ... ]

        self.__init_parameters()

        self.__init_model(data_instances)

        mini_batch_obj = MiniBatch(data_inst=data_instances, batch_size=self.batch_size)

        for iter_num in range(self.max_iter):
            # mini-batch
            batch_data_generator = mini_batch_obj.mini_batch_data_generator()
            total_loss = 0
            batch_num = 0

            for batch_data in batch_data_generator:
                n = batch_data.count()

                f = functools.partial(self.gradient_operator.compute,
                                      coef=self.coef_,
                                      intercept=self.intercept_,
                                      fit_intercept=self.fit_intercept)
                grad_loss = batch_data.mapPartitions(f)

                grad, loss = grad_loss.reduce(self.aggregator.aggregate_grad_loss)

                grad /= n
                loss /= n

                if self.updater is not None:
                    loss_norm = self.updater.loss_norm(self.coef_)
                    total_loss += (loss + loss_norm)
                delta_grad = self.optimizer.apply_gradients(grad)

                self.update_model(delta_grad)
                batch_num += 1

            total_loss /= batch_num
            w = self.merge_model()
            self.loss_history.append(total_loss)
            LOGGER.info("iter: {}, loss: {}".format(iter_num, total_loss))
            # send model
            model_transfer_id = self.transfer_variable.generate_transferid(self.transfer_variable.guest_model,
                                                                           iter_num)
            federation.remote(w,
                              name=self.transfer_variable.guest_model.name,
                              tag=model_transfer_id,
                              role=consts.ARBITER,
                              idx=0)

            # send loss

            loss_transfer_id = self.transfer_variable.generate_transferid(self.transfer_variable.guest_loss, iter_num)
            federation.remote(total_loss,
                              name=self.transfer_variable.guest_loss.name,
                              tag=loss_transfer_id,
                              role=consts.ARBITER,
                              idx=0)

            # recv model
            model_transfer_id = self.transfer_variable.generate_transferid(
                self.transfer_variable.final_model, iter_num)
            w = federation.get(name=self.transfer_variable.final_model.name,
                               tag=model_transfer_id,
                               idx=0)

            w = np.array(w)
            self.set_coef_(w)

            # recv converge flag
            converge_flag_id = self.transfer_variable.generate_transferid(self.transfer_variable.converge_flag,
                                                                          iter_num)
            converge_flag = federation.get(name=self.transfer_variable.converge_flag.name,
                                           tag=converge_flag_id,
                                           idx=0)

            self.n_iter_ = iter_num
            LOGGER.debug("converge flag is :{}".format(converge_flag))

            if converge_flag:
                self.is_converged = True
                break

        self.show_meta()
        self.show_model()
        LOGGER.debug("in fit self coef: {}".format(self.coef_))
        return data_instances

    def __init_parameters(self):
        party_weight_id = self.transfer_variable.generate_transferid(
            self.transfer_variable.guest_party_weight
        )
        federation.remote(self.party_weight,
                          name=self.transfer_variable.guest_party_weight.name,
                          tag=party_weight_id,
                          role=consts.ARBITER,
                          idx=0)

        # LOGGER.debug("party weight sent")
        LOGGER.info("Finish initialize parameters")

    def __init_model(self, data_instances):
        model_shape = data_overview.get_features_shape(data_instances)

        LOGGER.info("Initialized model shape is {}".format(model_shape))

        w = self.initializer.init_model(model_shape, init_params=self.init_param_obj)
        if self.fit_intercept:
            self.coef_ = w[:-1]
            self.intercept_ = w[-1]
        else:
            self.coef_ = w
            self.intercept_ = 0

        # LOGGER.debug("Initialed model")
        return w

    def predict(self, data_instances, predict_param):
        LOGGER.debug("coef: {}, intercept: {}".format(self.coef_, self.intercept_))
        wx = self.compute_wx(data_instances, self.coef_, self.intercept_)
        pred_prob = wx.mapValues(lambda x: activation.sigmoid(x))
        pred_label = self.classified(pred_prob, predict_param.threshold)

        if predict_param.with_proba:
            predict_result = data_instances.mapValues(lambda x: x.label)
            predict_result = predict_result.join(pred_prob, lambda x, y: (x, y))
        else:
            predict_result = data_instances.mapValues(lambda x: (x.label, None))

        predict_result = predict_result.join(pred_label, lambda x, y: (x[0], x[1], y))
        return predict_result

    def set_flowid(self, flowid=0):
        self.transfer_variable.set_flowid(flowid)
