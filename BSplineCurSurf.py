from math import *
import logging
import numpy as np

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
            #logging.info('i now is %s' % (i+1))
            temp[i] = sqrt(pow(Q[i+1][0]-Q[i][0],2)+pow(Q[i+1][1]-Q[i][1],2)+pow(Q[i+1][2]-Q[i][2],2))
        sumtemp = sum(temp) # 顺序相邻两型值点之间距离的和，即总弦长
        for i in range(p+1): # 前p+1个节点为0
            U[i] = 0
        for i in range(n+2, n+p+3): # 后p+1个节点为1
            U[i] = 1
        for i in range(p+1, n+p-1):
            U[i] = U[i-1] + temp[i-p-1] / sumtemp # 内节点
        return U
    
    # 找到基函数自变量u的相对于节点集合U所在区间的左端索引
    def findSpan(self, n, u, U): # n为u能在的最右区间的左端索引
        if (u == U[n+1]):
            return n
        low = self.p
        high = n + 1
        mid = floor((low+high)/2)
        while (u < U[mid] or u >= U[mid+1]):
            if (u < U[mid]):
                high = mid
            else:
                low = mid
            mid = floor((low+high)/2)
        return mid
    
    # 计算自变量u所在区间内的所有不为零的基函数N
    def basisFuns(self, i, u, U): # i为findSpan()的返回值
        N = [0 for i in range(self.p+1)]
        left = [0 for i in range(self.p+1)]
        right = [0 for i in range(self.p+1)]
        N[0] = 1.0
        for j in range(1,self.p+1):
            left[j] = u - U[i+1-j]
            right[j] = U[i+j] - u
            saved = 0.0
            for r in range(j):
                temp = N[r] / (right[r+1]+left[j-r])
                N[r] = saved + right[r+1] * temp
                saved = left[j-r] * temp
            N[j] = saved
        return N
    
    # 使用三对角矩阵算法解决反求B样条控制点P
    def solveTridiagonal(self, Q, U, P): # Q为型值点集，P为控制点集，假设P[0],P[1],P[n+1],P[n+2]已知
        n = len(Q) - 1 # 型值点个数 = n+1
        R = [[0 for i in range(3)] for i in range(n+1)]
        dd = [0 for i in range(n+1)]
        for i in range(3,n):
            R[i] = Q[i-1]
        abc = self.basisFuns(4, U[4], U)
        den = abc[1]
        P[2] = ((np.array(Q[1])-abc[0]*np.array(P[1])) / den).tolist()
        for i in range(3,n):
            dd[i] = abc[2] / den
            abc = self.basisFuns(i+2,U[i+2],U)
            den = abc[1] - abc[0] * dd[i]
            P[i] = ((np.array(R[i])-abc[0]*np.array(P[i-1])) / den).tolist()
        dd[n] = abc[2] / den
        abc = self.basisFuns(n+2,U[n+2],U)
        den = abc[1] - abc[0] * dd[n]
        P[n] = ((np.array(Q[n-1])-abc[2]*np.array(P[n+1])-abc[0]*np.array(P[n-1])) / den).tolist()
        for i in range(n-1,1,-1):
            P[i] = (np.array(P[i]) - dd[i+1] * np.array(P[i+1])).tolist()
        return P

# 测试
if __name__ == '__main__':
    pt = np.loadtxt('data2.txt') # 下载型值点数据
    Q = pt.tolist() # 把numpy变为list
    logging.info(Q[2])
    logging.info('---测试paraAndNodVect---')
    BS = BSplineCurSurf(3) # 基函数阶数p=3
    U = BS.paraAndNodVect(Q) # 节点集合
    logging.info(U)
    logging.info('---测试findSpan---')
    n = len(Q) + 1 # n为u能在的最右区间的左端索引
    u = 0.7
    i = BS.findSpan(n, u, U) # u 所在区间的左端索引值
    logging.info(i)
    logging.info('---测试basisFuns---')
    N = BS.basisFuns(i, u, U) # 计算i所在区间内所有不为零的基函数N
    logging.info(N)
    logging.info('---测试solveTridiagonal---')
    P = [[0 for i in range(3)] for i in range(len(Q)+2)] # 初始化控制点集
    P[0] = Q[0] 
    P[len(Q)+1] = Q[-1]
    P[1] = P[0] # 使用自由端点边界条件
    P[len(Q)] = P[len(Q)+1] # 使用自由端点边界条件
    P = BS.solveTridiagonal(Q, U, P)
    logging.info(P)
    




        
        
        




