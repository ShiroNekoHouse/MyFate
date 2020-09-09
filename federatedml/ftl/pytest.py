import time
import numpy as np
import random
from multiprocessing import Process
def run(name,i,res):
    print('%s runing %d' %(name,i))
    sum=0
    for j in range(0,10):
        sum+=10*j+j**2
    print('%s running end' %name)
    res[i]=sum

res=np.zeros((1,10))
processes = []

start_time=time.time()
for i in range(0,10):
    p=Process(target=run,args=('anne',i,res)) #必须加,号 
    p.start()
    processes.append(p)

for process in processes:
    process.join()

# for i in range(0,10):
#     run("anne",i)

end_time=time.time()
run_time=end_time - start_time
print('主线程')
print(run_time)

# for i in range(0,res.shape[0]):
#     run('anne',i,res)
#     print(res[i])