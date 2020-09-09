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

from arch.api import federation
from arch.api.utils import log_utils
from federatedml.logistic_regression.base_logistic_regression import BaseLogisticRegression
from federatedml.optim import DiffConverge
from federatedml.optim import activation
from federatedml.optim.federated_aggregator import HomoFederatedAggregator
from federatedml.param.param import LogisticParam
from federatedml.secureprotol import PaillierEncrypt, FakeEncrypt
from federatedml.util import consts
from federatedml.util.transfer_variable import HomoLRTransferVariable

LOGGER = log_utils.getLogger()


class HomoLRArbiter(BaseLogisticRegression):
    def __init__(self, params: LogisticParam):
        super(HomoLRArbiter, self).__init__(params)
        self.re_encrypt_batches = params.re_encrypt_batches
        self.aggregator = HomoFederatedAggregator()
        if params.converge_func == 'diff':
            self.convege_func = DiffConverge(eps=self.eps)
        else:
            raise RuntimeWarning("Cannot recognize converge_func, must be 'eps'.")
        self.transfer_variable = HomoLRTransferVariable()

        self.predict_threshold = params
        self.encrypt_param = params.encrypt_param
        self.classes_ = [0, 1]

        # To be initialized
        self.host_use_encryption = []
        self.re_encrypt_times = []  # Record the times needed for each host
        self.curt_re_encrypt_times = []
        self.host_encrypter = []
        self.party_weights = []  # The first one is guest weight, host weights for otherwise
        self.has_sychronized_encryption = False
        self.loss_history = []
        self.is_converged = False
        self.header = []

    def fit(self, data=None):
        LOGGER.debug("self.has_sychronized_encryption: {}".format(self.has_sychronized_encryption))
        self.__init_parameters()
        LOGGER.debug("self.has_sychronized_encryption: {}".format(self.has_sychronized_encryption))

        LOGGER.info("Finish init parameters")

        for iter_num in range(self.max_iter):
            # re_encrypt host models
            self.__re_encrypt(iter_num)

            # Part3: Aggregate models receive from each party
            final_model = self.aggregator.aggregate_model(transfer_variable=self.transfer_variable,
                                                          iter_num=iter_num,
                                                          party_weights=self.party_weights,
                                                          host_encrypter=self.host_encrypter)
            total_loss = self.aggregator.aggregate_loss(transfer_variable=self.transfer_variable,
                                                        iter_num=iter_num,
                                                        party_weights=self.party_weights,
                                                        host_use_encryption=self.host_use_encryption)
            self.loss_history.append(total_loss)
            LOGGER.info("Iter: {}, loss: {}".format(iter_num, total_loss))
            # send model
            final_model_id = self.transfer_variable.generate_transferid(self.transfer_variable.final_model, iter_num)
            federation.remote(final_model,
                              name=self.transfer_variable.final_model.name,
                              tag=final_model_id,
                              role=consts.GUEST,
                              idx=0)
            for idx, encrypter in enumerate(self.host_encrypter):
                encrypted_model = encrypter.encrypt_list(final_model)

                federation.remote(encrypted_model,
                                  name=self.transfer_variable.final_model.name,
                                  tag=final_model_id,
                                  role=consts.HOST,
                                  idx=idx)

            # send converge flag
            converge_flag = self.convege_func.is_converge(total_loss)
            converge_flag_id = self.transfer_variable.generate_transferid(
                self.transfer_variable.converge_flag,
                iter_num)

            federation.remote(converge_flag,
                              name=self.transfer_variable.converge_flag.name,
                              tag=converge_flag_id,
                              role=consts.GUEST,
                              idx=0)
            federation.remote(converge_flag,
                              name=self.transfer_variable.converge_flag.name,
                              tag=converge_flag_id,
                              role=consts.HOST,
                              idx=-1)
            self.set_coef_(final_model)
            self.n_iter_ = iter_num
            if converge_flag:
                self.is_converged = True
                break
        self._set_header()

    def predict(self, data=None, predict_param=None):
        # synchronize encryption information
        if not self.has_sychronized_encryption:
            self.__synchronize_encryption()
            self.__send_host_mode()

        for idx, use_encrypt in enumerate(self.host_use_encryption):
            if use_encrypt:
                encrypter = self.host_encrypter[idx]
                predict_wx_id = self.transfer_variable.generate_transferid(self.transfer_variable.predict_wx)
                predict_wx = federation.get(name=self.transfer_variable.predict_wx.name,
                                            tag=predict_wx_id,
                                            idx=idx
                                            )
                decrypted_wx = encrypter.distribute_decrypt(predict_wx)
                pred_prob = decrypted_wx.mapValues(lambda x: activation.sigmoid(x))
                pred_label = self.classified(pred_prob, predict_param.threshold)
                predict_result_id = self.transfer_variable.generate_transferid(self.transfer_variable.predict_result)
                federation.remote(pred_label,
                                  name=self.transfer_variable.predict_result.name,
                                  tag=predict_result_id,
                                  role=consts.HOST,
                                  idx=idx)
        LOGGER.info("Finish predicting, result has been sent back")
        return

    def __init_parameters(self):
        """
        This function is used to synchronized the parameters from each guest and host.
        :return:
        """
        # 1. Receive the party weight of each party
        # LOGGER.debug("To receive guest party weight")
        party_weight_id = self.transfer_variable.generate_transferid(
            self.transfer_variable.guest_party_weight
        )
        guest_weight = federation.get(name=self.transfer_variable.guest_party_weight.name,
                                      tag=party_weight_id,
                                      idx=0)

        # LOGGER.debug("Received guest_weight: {}".format(guest_weight))
        host_weight_id = self.transfer_variable.generate_transferid(
            self.transfer_variable.host_party_weight
        )
        host_weights = federation.get(name=self.transfer_variable.host_party_weight.name,
                                      tag=host_weight_id,
                                      idx=-1)
        weights = [guest_weight]
        weights.extend(host_weights)

        self.party_weights = [x / sum(weights) for x in weights]

        # 2. Synchronize encryption information
        self.__synchronize_encryption()

        # 3. Receive re-encrypt-times
        self.re_encrypt_times = [0] * len(self.host_use_encryption)
        for idx, use_encryption in enumerate(self.host_use_encryption):
            if not use_encryption:
                self.re_encrypt_times[idx] = 0
                continue
            re_encrypt_times_id = self.transfer_variable.generate_transferid(
                self.transfer_variable.re_encrypt_times
            )
            re_encrypt_times = federation.get(name=self.transfer_variable.re_encrypt_times.name,
                                              tag=re_encrypt_times_id,
                                              idx=idx)
            self.re_encrypt_times[idx] = re_encrypt_times
        LOGGER.info("re encrypt times for all parties: {}".format(self.re_encrypt_times))

    def __synchronize_encryption(self):
        """
        Communicate with hosts. Specify whether use encryption or not and transfer the public keys.
        """
        # 1. Use Encrypt: Specify which host use encryption
        host_use_encryption_id = self.transfer_variable.generate_transferid(
            self.transfer_variable.use_encrypt
        )
        host_use_encryption = federation.get(name=self.transfer_variable.use_encrypt.name,
                                             tag=host_use_encryption_id,
                                             idx=-1)
        self.host_use_encryption = host_use_encryption

        LOGGER.info("host use encryption: {}".format(self.host_use_encryption))
        # 2. Send pubkey to those use-encryption hosts
        for idx, use_encryption in enumerate(self.host_use_encryption):
            if not use_encryption:
                encrypter = FakeEncrypt()
            else:
                encrypter = PaillierEncrypt()
                encrypter.generate_key(self.encrypt_param.key_length)
                pub_key = encrypter.get_public_key()
                pubkey_id = self.transfer_variable.generate_transferid(self.transfer_variable.paillier_pubkey)
                federation.remote(pub_key, name=self.transfer_variable.paillier_pubkey.name,
                                  tag=pubkey_id, role=consts.HOST, idx=idx)
                # LOGGER.debug("send pubkey to host: {}".format(idx))

            self.host_encrypter.append(encrypter)
        self.has_sychronized_encryption = True

    def __send_host_mode(self):
        model = self.merge_model()
        final_model_id = self.transfer_variable.generate_transferid(self.transfer_variable.final_model, "predict")
        for idx, use_encrypt in enumerate(self.host_use_encryption):
            if use_encrypt:
                encrypter = self.host_encrypter[idx]
                final_model = encrypter.encrypt_list(model)
            else:
                final_model = model
            federation.remote(final_model,
                              name=self.transfer_variable.final_model.name,
                              tag=final_model_id,
                              role=consts.HOST,
                              idx=idx)

    def __re_encrypt(self, iter_num):
        # If use encrypt, model weight need to be re-encrypt every several batches.
        self.curt_re_encrypt_times = self.re_encrypt_times.copy()

        # Part2: re-encrypt model weight from each host
        batch_num = 0
        while True:
            batch_num += self.re_encrypt_batches

            to_encrypt_model_id = self.transfer_variable.generate_transferid(
                self.transfer_variable.to_encrypt_model, iter_num, batch_num
            )
            re_encrypted_model_id = self.transfer_variable.generate_transferid(
                self.transfer_variable.re_encrypted_model, iter_num, batch_num
            )
            for idx, left_times in enumerate(self.curt_re_encrypt_times):
                if left_times <= 0:
                    continue
                re_encrypt_model = federation.get(
                    name=self.transfer_variable.to_encrypt_model.name,
                    tag=to_encrypt_model_id,
                    idx=idx
                )
                encrypter = self.host_encrypter[idx]
                decrypt_model = encrypter.decrypt_list(re_encrypt_model)
                re_encrypt_model = encrypter.encrypt_list(decrypt_model)
                federation.remote(re_encrypt_model, name=self.transfer_variable.re_encrypted_model.name,
                                  tag=re_encrypted_model_id, role=consts.HOST, idx=idx)

                left_times -= 1
                self.curt_re_encrypt_times[idx] = left_times

            if sum(self.curt_re_encrypt_times) == 0:
                break

    def _set_header(self):
        self.header = ['head_' + str(x) for x in range(len(self.coef_))]
