from math import *
import logging
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

logging.basicConfig(level=logging.INFO)

class BSplineCurSurf(object):
    def __init__(self, p=3, q=3):
        self.p = p # 基函数阶数,暂时固定为3
        self.q = q # 涉及曲面问题时才用到，为纵向基函数阶数，本项目暂时固定为3
    
    # 求取节点向量U (已知型值点值Q)
    def paraAndNodVect(self, Q):
        p = self.p
        n = len(Q)
        U = [0 for i in range(n+p+3)]
        # 第i个型值点 Q[i]=[x,y,z] & 第i行j列的型值点 Q[i][j]=[x,y,z]
        temp = [0 for i in range(n-1)] # n个型值点，n-1段弦长
        for i in range(n-1):
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
    def findSpan(self, n, u, U):  # n为u能在的最右区间的左端索引，n=len(U)-5
        if (u == U[n+1]): #特殊情况
            return n
        # 二分搜索
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

    # 计算非有理B样条曲线上的对应某一自变量u的点值
    def curvePoint(self, n, U, P, u): # P为控制点集, n为u能在的最右区间的左端索引, n=len(U)-5
        p = self.p # 基函数阶数
        span = self.findSpan(n, u, U)
        N = self.basisFuns(span, u, U)
        C = [0 for i in range(3)]
        for i in range(p+1):
            C = (np.array(C) + N[i] * np.array(P[span-p+i])).tolist()
        return C

    # 已知型值点，得出要求密度点的曲线
    def curvePlot(self, n, Q, numpoi): # numpoi为要求的密度点数量，Q为型值点集, n为u能在的最右区间的左端索引, n=len(Q)+1=len(U)-5
        # 第一步：计算节点矢量
        U = self.paraAndNodVect(Q) 
        # 第二步：反算n+2个控制点，采用自由端点边界条件
        P = [[0 for i in range(3)] for i in range(len(Q)+2)] # 初始化控制点集
        P[0] = Q[0] 
        P[len(Q)+1] = Q[-1]
        P[1] = P[0] # 使用自由端点边界条件
        P[len(Q)] = P[len(Q)+1] # 使用自由端点边界条件
        P = self.solveTridiagonal(Q, U, P)
        # 第三步：求出要求密度点的曲线
        C = [[0 for i in range(3)] for i in range(numpoi)]
        for i in range(numpoi-1):
            C[i] = self.curvePoint(n, U, P, i/(numpoi-1))
        C[numpoi-1] = self.curvePoint(n, U, P, 1)
        return C

    # 已知节点向量和控制点向量，求曲面上的点
    def surfacePoint(self, n, U, m, V, P, u, v): # U、V分别为横向和纵向节点集合，P为控制点<P[i][j]=[x,y,z]>，u、v分别为横向和纵向基函数自变量
        #n = len(U) - 5 # n为横向维度上u能在的最右区间的左端索引
        #m = len(V) - 5 # m为纵向维度上v能在的最右区间的左端索引
        p = self.p # 横向基函数阶数
        q = self.q # 纵向基函数阶数，本项目q=p
        uspan = self.findSpan(n, u, U)
        Nu = self.basisFuns(uspan, u, U)
        vspan = self.findSpan(m, v, V)
        Nv = self.basisFuns(vspan, v, V)
        uind = uspan - p
        S = [0 for i in range(3)]
        for h in range(q+1):
            temp = [0 for i in range(3)]
            vind = vspan - q + h
            for k in range(p+1):
                temp = (np.array(temp) + Nu[k] * np.array(P[uind+k][vind])).tolist()
            S = (np.array(S) + Nv[h] * np.array(temp)).tolist()
        return S



# 测试
if __name__ == '__main__':
    pt = np.loadtxt('data2.txt') # 下载型值点数据
    Q = pt.tolist() # 把numpy变为list
    logging.info(Q[2])
    logging.info('---测试paraAndNodVect---')
    BS = BSplineCurSurf() # 基函数阶数p=3
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
    logging.info('---测试curvePoint---')
    C = BS.curvePoint(n, U, P, u)
    logging.info(C)
    logging.info('---测试curvePlot---')
    numpoi = 101
    Cur = BS.curvePlot(n, Q, numpoi) # 求得曲线上的点集
    # 画图
    fig = plt.figure()
    ax = Axes3D(fig)
    CurA = np.array(Cur) # 把list变为array
    np.savetxt('data2gen', CurA) # 保存曲线点集数据到.txt文件
    ax.scatter((CurA[:,0]).tolist(), (CurA[:,1]).tolist(), (CurA[:,2]).tolist(), c='r')
    # 添加坐标轴标记及坐标标题
    ax.set_xlabel('X',fontdict={'size':15,'color':'black'})
    ax.set_ylabel('Y',fontdict={'size':15,'color':'black'})
    ax.set_zlabel('Z',fontdict={'size':15,'color':'black'})
    plt.show() # 显示图像
    #plt.pause(1) # 暂停3秒
    #plt.close() # 关闭当前显示的图像
    logging.info('---测试surfacePoint---')
    n = 4
    m = 3
    U = [0,0,0,0,1/2,1,1,1,1]
    V = [0,0,0,0,1,1,1,1]
    P = [[[0 for i in range(3)] for i in range(4)] for i in range(5)]
    P[0][0] = [0,0,0];P[1][0] = [3,0,3];P[2][0] = [6,0,3];P[3][0] = [9,0,0];P[4][0] = [12,0,3]
    P[0][1] = [0,2,2];P[1][1] = [3,2,5];P[2][1] = [6,2,5];P[3][1] = [9,2,2];P[4][1] = [12,2,5]
    P[0][2] = [0,4,0];P[1][2] = [3,4,3];P[2][2] = [6,4,3];P[3][2] = [9,4,0];P[4][2] = [12,4,3]
    P[0][3] = [0,6,2];P[1][3] = [3,6,5];P[2][3] = [6,6,5];P[3][3] = [9,6,2];P[4][3] = [12,6,5]
    num = 100
    S = [[0 for i in range(num+1)] for i in range(num+1)]
    for i in range(num+1):
        logging.info(i)
        u = i / num
        for j in range(num+1):
            v = j / num
            S[i][j] = BS.surfacePoint(n, U, m, V, P, u, v)
    # 画图
    fig = plt.figure()
    ax = Axes3D(fig)
    SA = np.array(S) # 把list变为array
    #np.savetxt('data2gen', CurA) # 保存曲线点集数据到.txt文件
    ax.scatter((SA[:,:,0]).tolist(), (SA[:,:,1]).tolist(), (SA[:,:,2]).tolist(), c='r')
    # 添加坐标轴标记及坐标标题
    ax.set_xlabel('X',fontdict={'size':15,'color':'black'})
    ax.set_ylabel('Y',fontdict={'size':15,'color':'black'})
    ax.set_zlabel('Z',fontdict={'size':15,'color':'black'})
    plt.show() # 显示图像

    
