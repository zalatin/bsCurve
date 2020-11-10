from math import *
import logging

logging.basicConfig(level=logging.INFO)

class BSplineCurSurf(object):
    def __init__(self, p):
        self.p = p # 基函数阶数
    
    # 求取节点向量U (已知型值点值Q)
    def paraAndNodVect(self, Q):
        p = self.p
        n = len(Q)
        U = [0 for i in range(n+p+3)]
        # 第i个型值点 Q[i]=[x,y,z] & 第i行j列的型值点 Q[i][j]=[x,y,z]
        temp = [0 for i in range(n-1)] # n个型值点，n-1段弦长
        for i in range(n-1):
            temp[i] = sqrt(pow(Q[i+1][1]-Q[i][1],2)+pow(Q[i+1][2]-Q[i][2],2)+pow(Q[i+1][3]-Q[i][3],2))
        sumtemp = sum(temp) # 顺序相邻两型值点之间距离的和，即总弦长
        for i in range(p+1): # 前p+1个节点为0
            U[i] = 0
        for i in range((n+2)+1, n+p+3): # 后p+1个节点为1
            U[i] = 1
        for i in range(p+1, n+p-2):
            U[i+1] = U[i] + temp[i-p] / sumtemp
        return U
    
# 测试
if __name__ == '__main__':
    a = 5
    logging.info('hello %s' % a)


        
        
        




