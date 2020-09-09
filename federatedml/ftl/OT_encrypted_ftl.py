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

#port set as 7891

import numpy as np
from federatedml.ftl import mulmatOT
from arch.api.utils import log_utils
import time
import multiprocessing

LOGGER = log_utils.getLogger()

# from federatedml.optim.activation import sigmoid

class MulmatProcess(multiprocessing.Process):

    def __init__(self, role,port,b,expandsA,ROW_B,COL_B, queue):
        super().__init__()
        self.role=role
        self.port=port
        self.b=b
        self.expandsA=expandsA
        self.ROW_B=ROW_B
        self.COL_B=COL_B
        self.queue = queue

    def mulMat(self):
        res=mulmatOT.mulmatOT(self.role,self.port,self.b,self.expandsA,self.ROW_B,self.COL_B)
        # 将结果放入队列
        self.queue.put(res)

    def run(self):
        self.mulMat()


def sigmoid(x):
    return 1. / (1. + np.exp(-x))

#New!!把矩阵拉伸成向量
def stretching(matA):
    rowA=matA.shape[0]
    colA=matA.shape[1]
    matB=np.zeros((1,rowA*colA))
    for i in range(0,rowA):
        for j in range(0,colA):
            matB[0][i*colA+j]=matA[i][j]
    return matB

class PartyModelInterface(object):

    def send_components(self):
        pass

    def receive_components(self, components):
        pass

    def send_gradients(self):
        pass

    def receive_gradients(self, gradients):
        pass

    def predict(self, X):
        pass


class OTEncryptedFTLGuestModel(PartyModelInterface): #uB*phi, phi as a choice

    def __init__(self, local_model, model_param, is_trace=False):
        super(OTEncryptedFTLGuestModel, self).__init__()
        self.localModel = local_model
        self.feature_dim = local_model.get_encode_dim()
        self.alpha = model_param.alpha
        self.is_trace = is_trace
        self.logger = LOGGER
        self.expands=[] #new 记录所有参与矩阵运算的矩阵的expands和大小ROW,COL

    def set_batch(self, X, y, non_overlap_indexes=None, overlap_indexes=None):
        self.X = X
        self.y = y
        self.non_overlap_indexes = non_overlap_indexes
        self.overlap_indexes = overlap_indexes
        self.phi = None

    def __compute_phi(self, uA, y):
        length_y = len(y)
        return np.expand_dims(np.sum(y * uA, axis=0) / length_y, axis=0)

    #======= for loss:here we only need y_overlap, and phi; We need to calculate matmul(uB_overlap,phi.transpose) by OT matrix multiplication .and uA_overlap*uB_overlap by OT float multiplication
    #======= for gradient to host: we need y_overlap_2_phi_2, y_overlap_phi, mapping_comp_A
    #======= for gradient to guest: we need y_overlap_2_phi, and y_overlap. We need to calculate matmul(y_overlap_2_phi, uB_overlap_2) by OT matrix multiplication; 
    # and y_overlap * uB_overlap by OT float multiplication
    def _compute_components(self):
        self.uA = self.localModel.transform(self.X)
        # phi has shape (1, feature_dim)
        # phi_2 has shape (feature_dim, feature_dim)
        self.phi = self.__compute_phi(self.uA, self.y)
        self.phi_2 = np.matmul(self.phi.transpose(), self.phi)

        # y_overlap and y_overlap_2 have shape (len(overlap_indexes), 1)
        self.y_overlap = self.y[self.overlap_indexes]
        self.y_overlap_2 = self.y_overlap * self.y_overlap

        if self.is_trace:
            self.logger.debug("phi shape" + str(self.phi.shape))
            self.logger.debug("phi_2 shape" + str(self.phi_2.shape))
            self.logger.debug("y_overlap shape" + str(self.y_overlap.shape))
            self.logger.debug("y_overlap_2 shape" + str(self.y_overlap_2.shape))

        # following two parameters will be sent to host
        # y_overlap_2_phi_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # y_overlap_phi has shape (len(overlap_indexes), feature_dim)
        self.y_overlap_2_phi_2 = 0.25 * np.expand_dims(self.y_overlap_2, axis=2) * self.phi_2
        self.y_overlap_phi = -0.5 * self.y_overlap * self.phi

        self.uA_overlap = self.uA[self.overlap_indexes]
        # mapping_comp_A has shape (len(overlap_indexes), feature_dim)
        self.mapping_comp_A = - self.uA_overlap / self.feature_dim
        # self.special = - self.uA_overlap / self.feature_dim #先传过去看看New!!!

        if self.is_trace:
            self.logger.debug("y_overlap_2_phi_2 shape" + str(self.y_overlap_2_phi_2.shape))
            self.logger.debug("y_overlap_phi shape" + str(self.y_overlap_phi.shape))
            self.logger.debug("mapping_comp_A shape" + str(self.mapping_comp_A.shape))

        #compute the expands of the matrix and send it to the host
        uA_overlap = - self.uA_overlap / self.feature_dim
        self.uA_overlap_expands=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,uA_overlap)
        phit=self.phi.transpose()
        self.phit_expands=send_expand_matrix(self.feature_dim,1,phit)
        
        self.y_overlap_2_phi_expands=[]
        self.y_overlap_2_phi_2_expands=[]

        matrix_y_overlap=np.zeros((len(self.overlap_indexes),self.feature_dim))
        for i in range(0,len(self.overlap_indexes)): #倒过来
            for j in range(0,self.feature_dim):
                matrix_y_overlap[i][j]=self.y_overlap[i]
        self.y_overlap_expands=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,matrix_y_overlap)

        y_overlap_2_phi = np.expand_dims(self.y_overlap_2 * self.phi, axis=1)
        for i in range(0,len(self.overlap_indexes)):
            expands=send_expand_matrix(1,self.feature_dim,y_overlap_2_phi[i])
            self.y_overlap_2_phi_expands.append(expands)
        for i in range(0,len(self.overlap_indexes)):
            expands=send_expand_matrix(self.feature_dim,self.feature_dim,self.y_overlap_2_phi_2[i])
            self.y_overlap_2_phi_2_expands.append(expands)

    #new
    # #获取精度为4的放大缩小矩阵,在fit里循环调用这个函数来将所有需要的矩阵的expends和大小集合发送过去。
    # def send_expand_matrix(self,ROW,COL,matrix):
    #     matrixInt,matrixIntTimes=mulmatOT.ChaneToInt(ROW,COL,matrix)
    #     return [ROW,COL,matrixIntTimes]

    # #收取对方发送过来的expands集合
    # def receive_expand_matrix(self,expandsReceived):
    #     return expandsReceived

    # def compute_OT(self):
    #     #Loss_share1: matmul(uB_overlap,phi.transpose) ,guest as choice of OT,expands index=0
    #     uB_phi_share=mulmatOT.mulmatOT(1,phi.transpose,self.expands[0],self.feature_dim,1)
    #     loss_share1=-0.5 * np.sum(y_overlap*uB_phi_share)

    def send_components(self):
        #计算需要的变量
        self._compute_components()
        #计算参与OT的矩阵的放大缩小矩阵,ROW,COL
        #发送放大缩小矩阵,ROWB
        return [self.y_overlap_2_phi_2, self.y_overlap_phi, self.mapping_comp_A, self.uA_overlap_expands,self.phit_expands,self.y_overlap_2_phi_expands,self.y_overlap_2_phi_2_expands,self.y_overlap_expands]#New!!!这里应该把component全部换掉，矩阵的信息等
        #计算OT矩阵乘法与OT float，获得share，并加起来
        #发送share

    def receive_components(self, components):
        self.uB_overlap = components[0]
        self.uB_overlap_2 = components[1]
        self.uB_overlap_expands= components[2]
        self.uB_overlap_2_expands = components[3]
        self.uB_overlap_ex_expands = components[4]
        self.mapping_comp_B = components[5]
        # self._update_gradients()#New！！！这里要独立出去
        # self._update_loss()

    def _update_gradients(self):

        # y_overlap_2 have shape (len(overlap_indexes), 1),
        # phi has shape (1, feature_dim),
        # y_overlap_2_phi has shape (len(overlap_indexes), 1, feature_dim)
        y_overlap_2_phi = np.expand_dims(self.y_overlap_2 * self.phi, axis=1)

        # uB_overlap_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # loss_grads_const_part1 has shape (len(overlap_indexes), feature_dim)
        loss_grads_const_part1 = 0.25 * np.squeeze(np.matmul(y_overlap_2_phi, self.uB_overlap_2), axis=1)

        # loss_grads_const_part2 has shape (len(overlap_indexes), feature_dim)
        loss_grads_const_part2 = self.y_overlap * self.uB_overlap

        # print("reallossconst=-===")
        # print(loss_grads_const_part1)
        # print(loss_grads_const_part2)

        if self.is_trace:
            self.logger.debug("loss_grads_const_part1 shape" + str(loss_grads_const_part1.shape))
            self.logger.debug("loss_grads_const_part2 shape" + str(loss_grads_const_part2.shape))
            self.logger.debug("y_overlap shape" + str(self.y_overlap.shape))
            self.logger.debug("uB_overlap shape" + str(self.uB_overlap.shape))

        const = np.sum(loss_grads_const_part1, axis=0) - 0.5 * np.sum(loss_grads_const_part2, axis=0)
        # grad_A_nonoverlap has shape (len(non_overlap_indexes), feature_dim)
        # grad_A_overlap has shape (len(overlap_indexes), feature_dim)
        grad_A_nonoverlap = self.alpha * const * self.y[self.non_overlap_indexes] / len(self.y)
        grad_A_overlap = self.alpha * const * self.y_overlap / len(self.y) + self.mapping_comp_B

        loss_grad_A = np.zeros((len(self.y), self.uB_overlap.shape[1]))
        loss_grad_A[self.non_overlap_indexes, :] = grad_A_nonoverlap
        loss_grad_A[self.overlap_indexes, :] = grad_A_overlap
        loss_grads = loss_grad_A
        print("reallossgrads====")
        print(loss_grads)
        #self.localModel.backpropogate(self.X, self.y, loss_grad_A)

    #New!!!
    def prepare_gradient(self):
        # y_overlap_2 have shape (len(overlap_indexes), 1),
        # phi has shape (1, feature_dim),
        # y_overlap_2_phi has shape (len(overlap_indexes), 1, feature_dim)
        y_overlap_2_phi = np.expand_dims(self.y_overlap_2 * self.phi, axis=1)

        # uB_overlap_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # loss_grads_const_part1 has shape (len(overlap_indexes), feature_dim)
        res1=np.zeros((len(self.overlap_indexes), 1, self.feature_dim))
        processes=[]
        queues=[]
        for i in range(0,y_overlap_2_phi.shape[0]):
            matrixB=y_overlap_2_phi[i]
            # matrixA=self.uB_overlap_2[i]
            # msgfromA=send_expand_matrix(self.feature_dim,self.feature_dim,matrixA)
            # expandsA=receive_expand_matrix(msgfromA)
            expandsA=self.uB_overlap_2_expands[i]

            # temp=mulmatOT.mulmatOT(0,7891,matrixB,expandsA,1,self.feature_dim)
            # res1[i]=temp

            #speed up with process
            queue=multiprocessing.Queue()
            queues.append(queue)
            p=MulmatProcess(0,7891+i,matrixB,expandsA,1,self.feature_dim,queue)
            p.start()
            processes.append(p)

        for i,process in enumerate(processes):
            process.join()
            res1[i]=queues[i].get()
            
            # p=Process(target=process_P,args=(0,7891+i,matrixB,expandsA,1,self.feature_dim,res1,i))
            # p.start()
            # processes.append(p)

        # for process in processes:
        #     process.join()
        
        loss_grads_const_part1_share=0.25 * np.squeeze(res1, axis=1)

        # matrixA=np.zeros((1,1))#可优化
        # matrixB=np.zeros((1,1))
        # res2=np.zeros((len(self.overlap_indexes),self.feature_dim))
        # for i in range(0,self.uB_overlap.shape[0]): #倒过来
        #     for j in range(0,self.uB_overlap.shape[1]):
        #         matrixA[0][0]=self.y_overlap[i]
        #         matrixB[0][0]=self.uB_overlap[i][j]
        #         msgfromB=send_expand_matrix(1,1,matrixB)
        #         expandsB=receive_expand_matrix(msgfromB)
        #         temp=mulmatOT.mulmatOT(0,matrixA,expandsB,1,1)
        #         res2[i][j]=temp[0][0]
        matrixA=np.zeros((len(self.overlap_indexes),self.feature_dim))
        for i in range(0,len(self.overlap_indexes)): #倒过来
            for j in range(0,self.feature_dim):
                matrixA[i][j]=self.y_overlap[i]
        # matrixB=self.uB_overlap
        # msgfromB=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,matrixB)
        expandsB=self.uB_overlap_expands
        res2=mulmatOT.mulmatOT_wise(0,7891,matrixA,expandsB,len(self.overlap_indexes),self.feature_dim)

        loss_grads_const_part2_share=res2
        # self.const_share = [loss_grads_const_part1_share,loss_grads_const_part2_share]
        self.const_share=np.sum(loss_grads_const_part1_share, axis=0) - 0.5 * np.sum(loss_grads_const_part2_share, axis=0)

    def compute_gradients(self,gradient_share_from_host):
        # const_part1=self.const_share[0]+gradient_share_from_host[0]
        # const_part2=self.const_share[1]+gradient_share_from_host[1]
        # print("lossconst")
        # print(const_part1)
        # print(const_part2)
        # const = np.sum(const_part1, axis=0) - 0.5 * np.sum(const_part2, axis=0)

        const=self.const_share+gradient_share_from_host
        # grad_A_nonoverlap has shape (len(non_overlap_indexes), feature_dim)
        # grad_A_overlap has shape (len(overlap_indexes), feature_dim)
        grad_A_nonoverlap = self.alpha * const * self.y[self.non_overlap_indexes] / len(self.y)
        grad_A_overlap = self.alpha * const * self.y_overlap / len(self.y) + self.mapping_comp_B

        loss_grad_A = np.zeros((len(self.y), self.feature_dim))
        loss_grad_A[self.non_overlap_indexes, :] = grad_A_nonoverlap
        loss_grad_A[self.overlap_indexes, :] = grad_A_overlap
        
        self.loss_grads = loss_grad_A
        print("loss_grads")
        print(self.loss_grads)
        self.localModel.backpropogate(self.X, self.y, loss_grad_A)

    def send_gradient_shares(self):
        self.assist_gradient()
        return self.l1_grad_B_share

    #assist gradient for host
    def assist_gradient(self):
        # uB_overlap_ex has shape (len(overlap_indexes), 1, feature_dim)
        # uB_overlap_ex = np.expand_dims(self.uB_overlap, axis=1)
        res1=np.zeros((len(self.overlap_indexes), 1, self.feature_dim))
        processes=[]
        queues=[]
        for i in range(0,len(self.overlap_indexes)):
            # matrixB=uB_overlap_ex[i]
            matrixA=self.y_overlap_2_phi_2[i]
            # msgfromB=send_expand_matrix(1,self.feature_dim,matrixB)
            expandsB=self.uB_overlap_ex_expands[i]
            # temp=mulmatOT.mulmatOT(1,7891,matrixA,expandsB,self.feature_dim,self.feature_dim)
            # res1[i]=temp

            #speed up with process
            queue=multiprocessing.Queue()
            queues.append(queue)
            p=MulmatProcess(1,7891+i,matrixA,expandsB,self.feature_dim,self.feature_dim,queue)
            p.start()
            processes.append(p)
        for i,process in enumerate(processes):
            process.join()
            res1[i]=queues[i].get()
            
        #     p=Process(target=process_P,args=(1,7891+i,matrixA,expandsB,self.feature_dim,self.feature_dim,res1,i))
        #     p.start()
        #     processes.append(p)

        # for process in processes:
        #     process.join()
            
        self.l1_grad_B_share=np.squeeze(res1, axis=1)+self.y_overlap_phi

    def send_loss(self):
        return self.loss

    def receive_loss(self, loss):
        self.loss = loss

    def _update_loss(self):
        uA_overlap = - self.uA_overlap / self.feature_dim
        loss_overlap = np.sum(uA_overlap * self.uB_overlap)
        # print("=======realLoss")
        # print(loss_overlap)
        loss_y = self.__compute_loss_y(self.uB_overlap, self.y_overlap, self.phi)
        #self.loss = self.alpha * loss_y + loss_overlap
        loss = self.alpha * loss_y + loss_overlap
        print("realloss:===")
        print(loss)

    #New!!
    def prepare_loss(self):
        self.prepare_loss_part1()
        self.prepare_loss_part2()

    def prepare_loss_part1(self):
        phit=self.phi.transpose()
        # msgfromB=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,self.uB_overlap)
        expandsB=self.uB_overlap_expands

        self.uB_phi_share=mulmatOT.mulmatOT(1,7891,phit,expandsB,self.feature_dim,1)
        # print("------------------result------------------")
        # print(self.uB_phi_share)

        # print("------------------phi-------------------")
        # print(phit)
        # print("------------------uB-------------------")
        # print(self.uB_overlap)

        # outres=np.matmul(self.uB_overlap,phit)
        # print("-=======================================-")
        # print(outres)

    def prepare_loss_part2(self):
        uA_overlap = - self.uA_overlap / self.feature_dim  #这是什么东西？
        # uB_overlap=self.uB_overlap

        # matrixA=np.zeros((1,1))#可优化
        # matrixB=np.zeros((1,1))
        # self.loss_overlap_share=0
        # for i in range(0,uA_overlap.shape[0]):
        #     for j in range(0,uA_overlap.shape[1]):
        #         matrixA[0][0]=uA_overlap[i][j]
        #         matrixB[0][0]=uB_overlap[i][j]
        #         msgfromB=send_expand_matrix(1,1,matrixB)
        #         expandsB=receive_expand_matrix(msgfromB)
        #         res=mulmatOT.mulmatOT(1,matrixA,expandsB,1,1)
        #         self.loss_overlap_share=self.loss_overlap_share+res[0][0]

        # msgfromB=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,uB_overlap)
        expandsB=self.uB_overlap_expands
        res=mulmatOT.mulmatOT_wise(1,7891,uA_overlap,expandsB,len(self.overlap_indexes),self.feature_dim)
        self.loss_overlap_share=np.sum(res)
        # print("lossA:===")
        # print(self.loss_overlap_share)


    def compute_loss(self,loss_share):
        # print("lossshare:===")
        # print(loss_share)
        # uA_overlap = - self.uA_overlap / self.feature_dim
        loss_overlap = self.loss_overlap_share+loss_share[1] 
        uB_phi=self.uB_phi_share+loss_share[0]
        loss_y = (-0.5 * np.sum(self.y_overlap * uB_phi) + 1.0 / 8 * np.sum(uB_phi * uB_phi)) + len(self.y_overlap) * np.log(2)
        loss2 = self.alpha * loss_y + loss_overlap
        self.loss=loss2
        # print("lossoverlap:===")
        # print(loss_overlap)
        # print("lossy:===")
        # print(loss_y)
        # print("loss:===")
        # print(loss2)
        

    def __compute_loss_y(self, uB_overlap, y_overlap, phi):
        # uB_phi has shape (len(overlap_indexes), 1)
        uB_phi = np.matmul(uB_overlap, phi.transpose())
        loss_y = (-0.5 * np.sum(y_overlap * uB_phi) + 1.0 / 8 * np.sum(uB_phi * uB_phi)) + len(y_overlap) * np.log(2)
        return loss_y

    def get_loss_grads(self):
        return self.loss_grads

    def predict(self, uB):
        if self.phi is None:
            self.uA = self.localModel.transform(self.X)
            self.phi = self.__compute_phi(self.uA, self.y)
        return sigmoid(np.matmul(uB, self.phi.transpose()))

    def restore_model(self, model_parameters):
        self.localModel.restore_model(model_parameters)

    def get_model_parameters(self):
        return self.localModel.get_model_parameters()


class OTEncryptedFTLHostModel(PartyModelInterface):

    def __init__(self, local_model, model_param, is_trace=False):
        super(OTEncryptedFTLHostModel, self).__init__()
        self.localModel = local_model
        self.feature_dim = local_model.get_encode_dim()
        self.alpha = model_param.alpha
        self.is_trace = is_trace
        self.logger = LOGGER

    def set_batch(self, X, overlap_indexes):
        self.X = X
        self.overlap_indexes = overlap_indexes


    #======= for gradient to host: we need uB_overlap_ex , to calculate matmul(uB_overlap_ex, y_overlap_2_phi_2) by OT matrix multiplication
    #======= for gradient to guset: we need uB_overlap, uB_overlap_2, mapping_comp_B
    def _compute_components(self):
        self.uB = self.localModel.transform(self.X)

        # following three parameters will be sent to guest
        # uB_overlap has shape (len(overlap_indexes), feature_dim)
        # uB_overlap_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # mapping_comp_B has shape (len(overlap_indexes), feature_dim)
        self.uB_overlap = self.uB[self.overlap_indexes]
        self.uB_overlap_2 = np.matmul(np.expand_dims(self.uB_overlap, axis=2), np.expand_dims(self.uB_overlap, axis=1))
        self.mapping_comp_B = - self.uB_overlap / self.feature_dim

        if self.is_trace:
            self.logger.debug("uB_overlap shape" + str(self.uB_overlap.shape))
            self.logger.debug("uB_overlap_2 shape" + str(self.uB_overlap_2.shape))
            self.logger.debug("mapping_comp_B shape" + str(self.mapping_comp_B.shape))

        self.uB_overlap_expands=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,self.uB_overlap)
        self.uB_overlap_2_expands=[]
        self.uB_overlap_ex_expands=[]
        for i in range(0,len(self.overlap_indexes)):
            expands=send_expand_matrix(self.feature_dim,self.feature_dim,self.uB_overlap_2[i])
            self.uB_overlap_2_expands.append(expands)
        uB_overlap_ex = np.expand_dims(self.uB_overlap, axis=1)
        for i in range(0,len(self.overlap_indexes)):
            expands=send_expand_matrix(1,self.feature_dim,uB_overlap_ex[i])
            self.uB_overlap_ex_expands.append(expands)
        
    def send_components(self):
        self._compute_components()
        return [self.uB_overlap, self.uB_overlap_2, self.uB_overlap_expands, self.uB_overlap_2_expands, self.uB_overlap_ex_expands,self.mapping_comp_B]

    def receive_components(self, components):
        self.y_overlap_2_phi_2 = components[0]
        self.y_overlap_phi = components[1]
        self.mapping_comp_A = components[2]
        # self.phi=components[3]#New!!!这里应该把component全部换掉，矩阵的信息等
        # self.y_overlap_2=components[4]#New!!!这里应该把component全部换掉，矩阵的信息等
        # self.uA_overlap=components[5]#New!!!这里应该把component全部换掉，矩阵的信息等,这个是特殊的uA
        # self.y_overlap=components[6]#New!!!这里应该把component全部换掉，矩阵的信息等
        # # print("--------------received uA------------")
        # # for i in range(0,6):
        # #     print("===================="+str(i)+"==================")
        # #     print(components[i])
        # #self._update_gradients()
        self.uA_overlap_expands = components[3]
        self.phit_expands= components[4]
        self.y_overlap_2_phi_expands= components[5]
        self.y_overlap_2_phi_2_expands= components[6]
        self.y_overlap_expands= components[7]

    def send_loss_shares(self):
        self.assist_compute_loss_part1()
        self.assist_compute_loss_part2()
        return [self.uB_phi_share,self.loss_overlap_share]
    #New!!!!
    def assist_compute_loss_part1(self):
        expandsA=self.phit_expands
        # begin_time=time.time()
        self.uB_phi_share=mulmatOT.mulmatOT(0,7891,self.uB_overlap,expandsA,len(self.overlap_indexes),self.feature_dim)
        # end_time=time.time()
        # print("2===========time===============")
        # print(end_time-begin_time)
        # print("------------------result------------------")
        # print(res1)
        # print("------------------phi-------------------")
        # print(phit)
        # print("------------------uB-------------------")
        # print(self.uB_overlap)
        
    def assist_compute_loss_part2(self):
        uB_overlap=self.uB_overlap
        # uA_overlap=self.uA_overlap
        # matrixA=np.zeros((1,1))#可优化
        # matrixB=np.zeros((1,1))
        # self.loss_overlap_share=0
        # for i in range(0,uA_overlap.shape[0]):
        #     for j in range(0,uA_overlap.shape[1]):
        #         matrixA[0][0]=uA_overlap[i][j]
        #         matrixB[0][0]=uB_overlap[i][j]
        #         msgfromA=send_expand_matrix(1,1,matrixA)
        #         expandsA=receive_expand_matrix(msgfromA)
        #         res=mulmatOT.mulmatOT(0,matrixB,expandsA,1,1)
        #         self.loss_overlap_share=self.loss_overlap_share+res[0][0]

        # msgfromA=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,uA_overlap)
        expandsA=self.uA_overlap_expands
        res=mulmatOT.mulmatOT_wise(0,7891,uB_overlap,expandsA,len(self.overlap_indexes),self.feature_dim)
        self.loss_overlap_share=np.sum(res)

        # print("lossB:===")
        # print(self.loss_overlap_share)
    #New!!!
    def send_gradient_shares(self):
        self.assist_gradient()
        return self.assist_const_share

    def assist_gradient(self):
        # y_overlap_2 have shape (len(overlap_indexes), 1),
        # phi has shape (1, feature_dim),
        # y_overlap_2_phi has shape (len(overlap_indexes), 1, feature_dim)
        # y_overlap_2_phi = np.expand_dims(self.y_overlap_2 * self.phi, axis=1)#本来获取不到的

        # uB_overlap_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # loss_grads_const_part1 has shape (len(overlap_indexes), feature_dim)
        #res1=y_overlap_2_phi
        res1=np.zeros((len(self.overlap_indexes), 1, self.feature_dim))
        processes=[]
        queues=[]
        for i in range(0,len(self.overlap_indexes)):
            # matrixB=y_overlap_2_phi[i]
            matrixA=self.uB_overlap_2[i]
            # msgfromB=send_expand_matrix(1,self.feature_dim,matrixB)
            expandsB=self.y_overlap_2_phi_expands[i]
            # temp=mulmatOT.mulmatOT(1,7891,matrixA,expandsB,self.feature_dim,self.feature_dim)
            # res1[i]=temp

            #speed up with process
            queue=multiprocessing.Queue()
            queues.append(queue)
            p=MulmatProcess(1,7891+i,matrixA,expandsB,self.feature_dim,self.feature_dim,queue)
            p.start()
            processes.append(p)
        for i,process in enumerate(processes):
            process.join()
            res1[i]=queues[i].get()

        #     p=Process(target=process_P,args=(1,7891+i,matrixA,expandsB,self.feature_dim,self.feature_dim,res1,i))
        #     p.start()
        #     processes.append(p)

        # for process in processes:
        #     process.join()

        # print("res1===")
        # print(res1)

        loss_grads_const_part1_share=0.25 * np.squeeze(res1, axis=1)

        # matrixA=np.zeros((1,1))#可优化
        # matrixB=np.zeros((1,1))
        # res2=np.zeros((len(self.overlap_indexes),self.feature_dim))
        # for i in range(0,self.uB_overlap.shape[0]): #倒过来
        #     for j in range(0,self.uB_overlap.shape[1]):
        #         matrixA[0][0]=self.y_overlap[i]
        #         matrixB[0][0]=self.uB_overlap[i][j]
        #         msgfromA=send_expand_matrix(1,1,matrixA)
        #         expandsA=receive_expand_matrix(msgfromA)
        #         temp=mulmatOT.mulmatOT(1,matrixB,expandsA,1,1)
        #         res2[i][j]=temp[0][0]

        # matrixA=np.zeros((len(self.overlap_indexes),self.feature_dim))
        # for i in range(0,self.uB_overlap.shape[0]): #倒过来
        #     for j in range(0,self.uB_overlap.shape[1]):
        #         matrixA[i][j]=self.y_overlap[i]
        matrixB=self.uB_overlap
        # msgfromA=send_expand_matrix(len(self.overlap_indexes),self.feature_dim,matrixA)
        expandsA=self.y_overlap_expands
        res2=mulmatOT.mulmatOT_wise(1,7891,matrixB,expandsA,len(self.overlap_indexes),self.feature_dim)

        loss_grads_const_part2_share=res2
        # self.assist_const_share = [loss_grads_const_part1_share,loss_grads_const_part2_share]
        self.assist_const_share=np.sum(loss_grads_const_part1_share, axis=0) - 0.5 * np.sum(loss_grads_const_part2_share, axis=0)


    def _update_gradients(self):
        # uB_overlap_ex has shape (len(overlap_indexes), 1, feature_dim)
        uB_overlap_ex = np.expand_dims(self.uB_overlap, axis=1)

        # y_overlap_2_phi_2 has shape (len(overlap_indexes), feature_dim, feature_dim)
        # uB_overlap_y_overlap_2_phi_2 has shape (len(overlap_indexes), 1, feature_dim)
        uB_overlap_y_overlap_2_phi_2 = np.matmul(uB_overlap_ex, self.y_overlap_2_phi_2)

        self.overlap_uB_y_overlap_2_phi_2 = np.squeeze(uB_overlap_y_overlap_2_phi_2, axis=1)
        # y_overlap_phi has shape (len(overlap_indexes), feature_dim)
        l1_grad_B = np.squeeze(uB_overlap_y_overlap_2_phi_2, axis=1) + self.y_overlap_phi
        loss_grad_B = self.alpha * l1_grad_B + self.mapping_comp_A
        loss_grads = loss_grad_B
        print("reallossgradB")
        print(loss_grads)
        # self.localModel.backpropogate(self.X[self.overlap_indexes], None, loss_grad_B)

    def prepare_gradient(self):
        # uB_overlap_ex has shape (len(overlap_indexes), 1, feature_dim)
        uB_overlap_ex = np.expand_dims(self.uB_overlap, axis=1)
        res1=uB_overlap_ex
        processes=[]
        queues=[]
        for i in range(0,uB_overlap_ex.shape[0]):
            matrixB=uB_overlap_ex[i]
            # matrixA=self.y_overlap_2_phi_2[i]
            # msgfromA=send_expand_matrix(self.feature_dim,self.feature_dim,matrixA)
            expandsA=self.y_overlap_2_phi_2_expands[i]
            # temp=mulmatOT.mulmatOT(0,7891,matrixB,expandsA,1,self.feature_dim)
            # res1[i]=temp

            #speed up with process
            queue=multiprocessing.Queue()
            queues.append(queue)
            p=MulmatProcess(0,7891+i,matrixB,expandsA,1,self.feature_dim,queue)
            p.start()
            processes.append(p)
        for i,process in enumerate(processes):
            process.join()
            res1[i]=queues[i].get()
        #     p=Process(target=process_P,args=(0,7891+i,matrixB,expandsA,1,self.feature_dim,res1,i))
        #     p.start()
        #     processes.append(p)

        # for process in processes:
        #     process.join()
        self.l1_grad_B_share=np.squeeze(res1, axis=1)

    def compute_gradients(self,guest_gradient_share_to_host):
        l1_grad_B=self.l1_grad_B_share+guest_gradient_share_to_host# + self.y_overlap_phi
        loss_grad_B = self.alpha * l1_grad_B + self.mapping_comp_A
        print("lossgradB===")
        print(loss_grad_B)
        self.loss_grads = loss_grad_B
        self.localModel.backpropogate(self.X[self.overlap_indexes], None, loss_grad_B)

    def get_loss_grads(self):
        return self.loss_grads

    def predict(self, X):
        return self.localModel.transform(X)

    def restore_model(self, model_parameters):
        self.localModel.restore_model(model_parameters)

    def get_model_parameters(self):
        return self.localModel.get_model_parameters()

def send_expand_matrix(ROW,COL,matrix):
    matrixInt,matrixIntTimes=mulmatOT.ChaneToInt(ROW,COL,matrix)
    return [ROW,COL,matrixIntTimes]

def receive_expand_matrix(msg):
    return msg


class LocalOTFederatedTransferLearning(object):

    def __init__(self, guest: OTEncryptedFTLGuestModel, host: OTEncryptedFTLHostModel, private_key=None):
        super(LocalOTFederatedTransferLearning, self).__init__()
        self.guest = guest
        self.host = host
        self.private_key = private_key

    def fit(self,role, X_A, X_B, y, overlap_indexes, non_overlap_indexes):
        self.guest.set_batch(X_A, y, non_overlap_indexes, overlap_indexes)
        self.host.set_batch(X_B, overlap_indexes)
        comp_B = self.host.send_components()
        comp_A = self.guest.send_components()
        self.guest.receive_components(comp_B)
        self.host.receive_components(comp_A)

        phit=self.guest.phi.transpose()
        msgfromA=send_expand_matrix(self.guest.feature_dim,1,phit)
        expandsA=receive_expand_matrix(msgfromA)
        msgfromB=send_expand_matrix(len(self.host.overlap_indexes),self.guest.feature_dim,self.host.uB_overlap)
        expandsB=receive_expand_matrix(msgfromB)
        
        # if role==0:
        #     res1=mulmatOT.mulmatOT(0,self.host.uB_overlap,expandsA,len(self.host.overlap_indexes),self.guest.feature_dim)
        #     print("------------------result------------------")
        #     print(res1)

        # else:
        #     res2=mulmatOT.mulmatOT(1,phit,expandsB,self.guest.feature_dim,1)
        #     print("------------------result------------------")
        #     print(res2)

        # outres=np.matmul(self.host.uB_overlap,phit)
        # print("-=======================================-")
        # print(outres)
        loss = self.guest.send_loss()
        return loss

    def predict(self, X_B):
        msg = self.host.predict(X_B)
        return self.guest.predict(msg)
