# -*- coding: utf-8 -*-
# A and B are matrixs, and A is in the front of B
import time
import sys,getopt
import threading
import numpy as np
import pickle
import ctypes
#这里
import random

BITS_LEN=16
NEW_BIT_LEN=0
so = ctypes.cdll.LoadLibrary   
lib = so("./federatedml/ftl/liboteIKNP.so")
#lib = so("./liboteIKNP.so")

# a=np.array([[100.1,2087,378.1,200.1,376.6,5801],[408.1,560.1,600.8,5007,108.6,250.8]])
# b=np.array([[145.5,27.71],[307.1,49.01],[550.1,652.1],[30.01,400.1],[370.1,405.1],[308.9,470.1]])
#c=np.array([[1,2],[3,4],[5,6]]) #random matrix c

def FloatToInt(f):
    i=f
    expanse=1
    ans=0
    if i>0:
        while(i*10<10000):
            i=i*10.0
            #print(i)
            expanse=expanse*10
        ans=int(round(i))
    else:
        i=-i
        while(i*10<10000):
            i=i*10
            expanse=expanse*10
        ans=int(round(i))
        expanse=-expanse
    #print(ans)
    return expanse, ans

def GetBTimes(BT):
    return BT

def CreateEij(i,j,m,n):
    #i，j is 1,others are 0
    Eij=[[0 for x in range(n)] for x in range(m)]
    Eij[i][j]=1
    return  Eij

def ChangeToBinary(n,m):
    #change n to binary with length of m
    binary=[0 for x in range(m)]
    index=0
    while(n>0):
        n, rem = divmod(n, 2)
        binary[index]=rem
        index+=1
    #print(binary)
    return binary
    #binary is the adverse order list of real one


def PreAddition():
    # A and B are matrixs
    pass

#change the float matrix to int matrix and note the times of expanse
def ChaneToInt(ROW,COL,matrix):
    matrixInt=np.zeros((ROW,COL),dtype=int)
    matrixIntTimes=np.zeros((ROW,COL))
    for i in range(0,ROW):
        for j in range(0,COL):
            matrixIntTimes[i][j],matrixInt[i][j]=FloatToInt(matrix[i][j])
    return matrixInt,matrixIntTimes

#Alice act as matrix breaker,and receiver of OT ,role 0
class Alice(object):
    def __init__(self,matrix,expands,ROW,COL):
        try:
            self.a=matrix
            self.ROW_B=expands[0]
            self.COL_B=expands[1]
            self.ROW_A=ROW
            self.COL_A=COL
            self.resArr=np.zeros((self.ROW_B,self.COL_A))
            self.zijkList=[0 for x in range(self.ROW_A*self.COL_A*BITS_LEN)]
            self.choice=[0 for x in range(self.ROW_A*self.COL_A*BITS_LEN)]
            self.aIntTimes=np.zeros((self.ROW_A,self.COL_A))
            self.aInt=np.zeros((self.ROW_A,self.COL_A),dtype=int)
            self.bIntTimes=expands[2]#tobechanged

            #这里
            self.resArr_wise=np.zeros((self.ROW_A,self.COL_A))
            #print ('Alice has been initialized\n')
        except:
            print('Fail')
            time.sleep(2)
            sys.exit()

    #change the float matrix to int matrix and note the times of expanse
    def ChaneToInt(self):
        for i in range(0,self.ROW_A):
            for j in range(0,self.COL_A):
                expanse,self.aInt[i][j]=FloatToInt(self.a[i][j])
                #print(self.aInt)
                self.aIntTimes[i][j]=expanse
        # for i in range(0,self.ROW_B):
        #     for j in range(0,self.COL_B):
        #         expanse,nothing=FloatToInt(b[i][j])
        #         self.bIntTimes[i][j]=expanse
    


    #choice with length of self.ROW_A*self.COL_A*BITS_LEN（在调用ot的时候每个choice都重复ROW_B*COL_B次）
    def MakeChoice(self,port):
        for i in range(0,self.ROW_A):
            for j in range(0,self.COL_A):
                jBinary=ChangeToBinary(self.aInt[i][j],BITS_LEN)
                for k in range(0,BITS_LEN):#TODO:try to use jBinary ,need 2-d vector
                    self.choice[self.COL_A*BITS_LEN*i+BITS_LEN*j+k]=jBinary[k]
                    #send to OT function
        self.CallOT(port)
        pass
    #这里
    def MakeChoice_wise(self):
        for i in range(0,self.ROW_A):
            for j in range(0,self.COL_A):
                jBinary=ChangeToBinary(self.aInt[i][j],BITS_LEN)
                for k in range(0,BITS_LEN):#TODO:try to use jBinary ,need 2-d vector
                    self.choice[self.COL_A*BITS_LEN*i+BITS_LEN*j+k]=jBinary[k]
                    #send to OT function
        
    def MultiMatrixByOT(self):
        # given a Matrix A from Alice,change it to sum(aijk*Eij)
        # A new matrix with element of aijk(binary list)
        for i in range(0,self.ROW_A):
            for j in range(0,self.COL_A):
                matrixSum=np.zeros((self.ROW_B,self.COL_B))
                indexBit=0
                for k in range(0,BITS_LEN):
                    zijk=self.GetZijk(i,j,k)
                    #print (zijk)
                    matrixSum=matrixSum+(zijk*(2**indexBit))
                    indexBit+=1
                matrixSum=matrixSum/self.aIntTimes[i][j]  #get back the times
                self.GetResult(matrixSum,i,j)
        pass
    
    #send jBinary as choice and get the zijklist with length of ROWA*COLA*8
    def CallOT(self,port):
        #fake c and cb
        # begin_time=time.time()
        int_arr2 = ctypes.c_int*2
        int_arr_c = int_arr2()
        int_arr_cb = int_arr2()

        int_arr_choice = ctypes.c_int*(self.ROW_A*self.COL_A*BITS_LEN)
        int_arr = int_arr_choice()
        for i in range(0,self.ROW_A*self.COL_A*BITS_LEN):
            int_arr[i]=self.choice[i]
        int_arr_res=ctypes.c_int*(self.ROW_A*self.COL_A*self.ROW_B*self.COL_B*BITS_LEN)
        int_res=int_arr_res()
        # end_time=time.time()
        # run_time=end_time-begin_time
        # print("=======OTpreparetime========")
        # print(run_time)

        begin_time=time.time()
        lib.pyotextension(1,port,int_arr_c,int_arr_cb,int_arr,self.ROW_A,self.COL_A,self.ROW_B,self.COL_B,int_res,NEW_BIT_LEN)
        # end_time=time.time()
        # run_time=end_time-begin_time
        # print("=======OTtime========")
        # print(run_time)

        # for i in range(0,self.ROW_A*self.COL_A*self.ROW_B*self.COL_B*BITS_LEN):
        #     print(int_res[i])

        # begin_time=time.time()

        #set int_res to zijklist per self.ROW_B*self.COL_B 这里容易出错，注意检查，尤其是行列的先后顺序，因为现在数字都设成一样
        count=0
        while count<self.ROW_A*self.COL_A*BITS_LEN:
            temp=np.zeros((self.ROW_B,self.COL_B))
            temp2=np.zeros((self.ROW_B,self.COL_B))
            for i in range(0,self.ROW_B):
                for j in range(0,self.COL_B):
                    temp[i][j]=int_res[count*self.ROW_B*self.COL_B+self.COL_B*i+j]/self.bIntTimes[i][j]
                    temp2[i][j]=int_res[count*self.ROW_B*self.COL_B+self.COL_B*i+j]
            # print("================received temp===============")
            # print(temp2)
            self.zijkList[count]=temp
            count+=1
        # count=0
        # while count<self.ROW_A*self.COL_A*BITS_LEN:
        #     temp=np.zeros((self.ROW_B,self.COL_B))
        #     for i in range(0,self.ROW_B):
        #         tempRes=int_res[count*self.ROW_B+i]
        #         for j in range(0,self.COL_B):
        #             #求余，右移
        #             temp[i][j]=tempRes%(2**(BITS_LEN))
        #             tempRes=tempRes>>BITS_LEN
        #     self.zijkList[count]=temp
        #     count+=1
        # end_time=time.time()
        # run_time=end_time-begin_time
        # print("=======resulttime========")
        # print(run_time)

    #这里
    def CallOT_wise(self,port):
        #fake c and cb
        int_arr2 = ctypes.c_int*2
        int_arr_c = int_arr2()
        int_arr_cb = int_arr2()

        int_arr_choice = ctypes.c_int*(self.ROW_A*self.COL_A*BITS_LEN)
        int_arr = int_arr_choice()
        for i in range(0,self.ROW_A*self.COL_A*BITS_LEN):
            int_arr[i]=self.choice[i]
        int_arr_res=ctypes.c_int*(self.ROW_A*self.COL_A*self.ROW_B*self.COL_B*BITS_LEN)
        int_res=int_arr_res()
        lib.pyotextension_wise(1,port,int_arr_c,int_arr_cb,int_arr,self.ROW_A,self.COL_A,self.ROW_B,self.COL_B,int_res,NEW_BIT_LEN)

        tempRes=np.zeros((self.ROW_B,self.COL_B))
        for i in range (0,self.ROW_A):
            for j in range (0,self.COL_A):
                sumij=0
                for k in range (0,BITS_LEN):
                    zk=int_res[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]
                    sumij=sumij+(zk*(2**k))
                sumij=sumij/self.bIntTimes[i][j]/self.aIntTimes[i][j]
                tempRes[i][j]=sumij
        self.resArr_wise=tempRes
    
    def GetZijk(self,i,j,k):
        return self.zijkList[self.COL_A*BITS_LEN*i+BITS_LEN*j+k]


    def GetResult(self,matrixS,I,J):
        #print("-------------2nd------------")
        #print(np.dot(matrixS,np.array(CreateEij(I,J,self.ROW_A,self.COL_A))))
        self.resArr=self.resArr+np.dot(matrixS,np.array(CreateEij(I,J,self.ROW_A,self.COL_A)))

    def ShowResult(self):
        return self.resArr
    
    #这里
    def ShowResult_wise(self):
        # print("------------------result------------------")
        # print(self.resArr_wise)
        return self.resArr_wise

#Bob act as server of OT , role 1
class Bob(object):
    def __init__(self,matrix,expands,ROW,COL):
        try:
            self.b=matrix
            self.ROW_A=expands[0]
            self.COL_A=expands[1]
            self.ROW_B=ROW
            self.COL_B=COL
            self.resArr=np.zeros((self.ROW_B,self.COL_A))
            self.randomMatrixList=[]
            self.bIntTimes=np.zeros((self.ROW_B,self.COL_B))
            self.bInt=np.zeros((self.ROW_B,self.COL_B),dtype=int)
            self.aIntTimes=expands[2]
            #这里
            self.randomMatrix_wise=[]
            self.resArr_wise=np.zeros((self.ROW_B,self.COL_B))
            # print ('Bob has been initialized\n')
        except:
            print('Fail')
            time.sleep(2)
            sys.exit()

    #change the float matrix to int matrix and note the times of expanse
    def ChaneToInt(self):
        for i in range(0,self.ROW_B):
            for j in range(0,self.COL_B):
                expanse,self.bInt[i][j]=FloatToInt(self.b[i][j])
                self.bIntTimes[i][j]=expanse
        # for i in range(0,self.ROW_A):
        #     for j in range(0,self.COL_A):
        #         expanse,nothing=FloatToInt(a[i][j])
        #         self.aIntTimes[i][j]=expanse

    #这里可以再优化，生成随机矩阵C
    def GetRandomC(self):
        for i in range(0,self.ROW_A*self.COL_A*BITS_LEN):
            temp=np.random.randint(1000,10000,(self.ROW_B,self.COL_B))#要改
            self.randomMatrixList.append(temp)
        #print(self.randomMatrixList)

    #这里
    def GetRandomC_wise(self):
        for i in range(0,self.ROW_B*self.COL_B*BITS_LEN):
            temp=random.randint(1000,10000)
            self.randomMatrix_wise.append(temp)


    def GetResult(self):
    # in the client of Bob, Bob input Matrix B and constant matrix C,
    # sending C,B+C as input of OT
    # setting -c as z1
        matrixSum=np.zeros((self.ROW_B,self.COL_B))
        for i in range (0,self.ROW_A):
            for j in range (0,self.COL_A):
                for k in range (0,BITS_LEN):
                    zijk=-1*self.randomMatrixList[i*self.COL_A*BITS_LEN+j*BITS_LEN+k]
                    matrixSum=matrixSum+(zijk*(2**k))
                for m in range(0,self.ROW_B):
                    for n in range(0,self.COL_B):
                        matrixSum[m][n]=matrixSum[m][n]/self.bIntTimes[m][n]  #get back of the times
                self.resArr=self.resArr+np.dot(matrixSum,np.array(CreateEij(i,j,self.ROW_A,self.COL_A)))/self.aIntTimes[i][j]  #this devide get back of the times is important!!!!!
                matrixSum=np.zeros((self.ROW_B,self.COL_B))
    #这里
    def GetResult_wise(self):
        tempRes=np.zeros((self.ROW_B,self.COL_B))
        for i in range (0,self.ROW_B):
            for j in range (0,self.COL_B):
                sumij=0
                for k in range (0,BITS_LEN):
                    zk=-1*self.randomMatrix_wise[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]
                    sumij=sumij+(zk*(2**k))
                sumij=sumij/self.bIntTimes[i][j]/self.aIntTimes[i][j]
                tempRes[i][j]=sumij
        self.resArr_wise=tempRes

    #send c and b+c to OT protocol
    def CallOT(self,port):
        #choice=[0 for x in range(2)] 
        #fake ,have no meaning ,to be improved for param
        int_arr2 = ctypes.c_int*2
        int_arr = int_arr2()
        #for i in range(0,2):
            #int_arr[i]=choice[i]#fake ,have no meaning ,to be improved for param

        #旧版本
        
        
        int_arr_matrixb = ctypes.c_int*(self.ROW_B*self.COL_B*self.ROW_A*self.COL_A*BITS_LEN)
        int_arr_c = int_arr_matrixb()
        int_arr_cb = int_arr_matrixb()
        #print(self.randomMatrixList)
        # begin_time=time.time()
        for k in range(0,self.ROW_A*self.COL_A*BITS_LEN):
            for i in range(0,self.ROW_B):
                for j in range(0,self.COL_B):
                    int_arr_c[k*self.ROW_B*self.COL_B+i*self.COL_B+j]=self.randomMatrixList[k][i][j]
                    int_arr_cb[k*self.ROW_B*self.COL_B+i*self.COL_B+j]=self.randomMatrixList[k][i][j]+self.bInt[i][j]
                    # print("!!!!!!!!!!!!!!!!!!!!!!!!")
                    # print(int_arr_c[k*self.ROW_B*self.COL_B+i*self.COL_B+j])
                    # print(int_arr_cb[k*self.ROW_B*self.COL_B+i*self.COL_B+j])
                    
        # end_time=time.time()
        # run_time=end_time-begin_time
        # print("=======preparetime========")
        # print(run_time)
        # #将B和B+C的m每一行统一成一个数，减少OT的数目
        # int_arr_matrixb = ctypes.c_uint64*(self.ROW_B*self.ROW_A*self.COL_A*BITS_LEN)
        # int_arr_c = int_arr_matrixb()
        # int_arr_cb = int_arr_matrixb()
        # for k in range(0,self.ROW_A*self.COL_A*BITS_LEN):
        #     for i in range(0,self.ROW_B):
        #         new_c=0
        #         new_bc=0
        #         for j in range(0,self.COL_B):
        #             #注意顺序
        #             if j!=0:
        #                 new_c=new_c<<BITS_LEN
        #                 new_bc=new_bc<<BITS_LEN
        #             new_c=new_c+self.randomMatrixList[k][i][self.COL_B-1-j]
        #             #print("c:")
        #             #print(self.randomMatrixList[k][i][self.COL_B-1-j])
        #             new_bc=new_bc+self.randomMatrixList[k][i][self.COL_B-1-j]+b[i][self.COL_B-1-j]
        #         int_arr_c[k*self.ROW_B+i]=new_c
        #         int_arr_cb[k*self.ROW_B+i]=new_bc
        #         #print(new_c)
        # # for i in range(0,self.ROW_B*self.ROW_A*self.COL_A*BITS_LEN):
        # #     print(int_arr_c[i])
        # #     #print(int_arr_cb[i])
        
        int_arrres=ctypes.c_uint64*(self.ROW_A*self.COL_A*self.ROW_B*self.COL_B*BITS_LEN)
        int_res=int_arrres()#fake ,have no meaning ,to be improved for param
        # begin_time=time.time()
        lib.pyotextension(0,port,int_arr_c,int_arr_cb,int_arr,self.ROW_A,self.COL_A,self.ROW_B,self.COL_B,int_res,NEW_BIT_LEN)
        # end_time=time.time()
        # run_time=end_time-begin_time
        # print("=======libtime========")
        # print(run_time)

    #这里
    def CallOT_wise(self,port):
        #choice=[0 for x in range(2)] 
        #fake ,have no meaning ,to be improved for param
        int_arr2 = ctypes.c_int*2
        int_arr = int_arr2()

        int_arr_matrixb = ctypes.c_int*(self.ROW_B*self.COL_B*BITS_LEN)
        int_arr_c = int_arr_matrixb()
        int_arr_cb = int_arr_matrixb()
        for i in range(0,self.ROW_B):
            for j in range(0,self.COL_B):
                for k in range(0,BITS_LEN):
                    int_arr_c[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]=self.randomMatrix_wise[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]
                    int_arr_cb[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]=self.randomMatrix_wise[i*self.COL_B*BITS_LEN+j*BITS_LEN+k]+self.bInt[i][j]
                    # print("!!!!!!!!!!!!!!!!!!!!!!!!")
                    # print(int_arr_c[k*self.ROW_B*self.COL_B+i*self.COL_B+j])
                    # print(int_arr_cb[k*self.ROW_B*self.COL_B+i*self.COL_B+j])
        
        int_arrres=ctypes.c_uint64*(self.ROW_A*self.COL_A*self.ROW_B*self.COL_B*BITS_LEN)
        int_res=int_arrres()#fake ,have no meaning ,to be improved for param
        lib.pyotextension_wise(0,port,int_arr_c,int_arr_cb,int_arr,self.ROW_A,self.COL_A,self.ROW_B,self.COL_B,int_res,NEW_BIT_LEN)   
        

    def ShowResult(self):   
        return self.resArr
    
    #这里
    def ShowResult_wise(self):
        # print("------------------result------------------")
        # print(self.resArr_wise)
        return self.resArr_wise


def mulmatOT(role,port,matrix,expands,ROW,COL):
    if role==0: #act as Bob.server in OT
        bob=Bob(matrix,expands,ROW,COL)
        bob.ChaneToInt()
        bob.GetRandomC()
        bob.CallOT(port)
        bob.GetResult()
        res=bob.ShowResult()
    else:       #act as Alice:receiver in OT
        alice=Alice(matrix,expands,ROW,COL)
        alice.ChaneToInt()
        alice.MakeChoice(port)
        alice.MultiMatrixByOT()
        res=alice.ShowResult()
    return res

#这里
def mulmatOT_wise(role,port,matrix,expands,ROW,COL):
    if role==0: #act as Bob.server in OT
        bob=Bob(matrix,expands,ROW,COL)
        bob.ChaneToInt()
        bob.GetRandomC_wise()
        bob.CallOT_wise(port)
        bob.GetResult_wise()
        res=bob.ShowResult_wise()
    else:       #act as Alice:receiver in OT
        alice=Alice(matrix,expands,ROW,COL)
        alice.ChaneToInt()
        alice.MakeChoice_wise()
        alice.CallOT_wise(port)
        res=alice.ShowResult_wise()
    return res

# def main(argv):
#     try:
#         opts, args = getopt.getopt(argv,"hr:")
#     except getopt.GetoptError:
#         print ('xxx.py -r (0 or 1)')
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt == '-h':
#             print ('xxx.py -r (0 or 1)')
#             sys.exit()
#         elif opt in ("-r"):
#             role = int(arg)
    
#     # hhh=np.dot(b,a)
#     # print(hhh)
      
#     if role==0: #act as Bob.server in OT
#         bob=Bob()
#         bob.ChaneToInt()
#         bob.GetRandomC()
#         bob.CallOT()
#         bob.GetResult()
#         bob.ShowResult()
#     else:       #act as Alice:receiver in OT
#         alice=Alice()
#         alice.ChaneToInt()
#         alice.MakeChoice()
#         alice.MultiMatrixByOT()
#         alice.ShowResult()

# if __name__ == "__main__":
#     main(sys.argv[1:])




