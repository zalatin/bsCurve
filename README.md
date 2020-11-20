# bsCurve

## 一、环境搭建
### 1 建立项目BS_Curve
### 2 配置虚拟环境venv
        py -3 -m venv .venv #安装虚拟环境到工作目录
        Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process #临时更改PowerShell执行策略
        .venv\scripts\activate #启动虚拟环境
### 3 安装软件包
        python -m pip install matplotlib #使用pip工具安装python软件工具包matplotlib


## 二、尝试编写Class BSplineCurSurf()：
        属性
             1 p ---基函数阶数
        方法
            1 findSpan(n, p, u, U) 
            2 basisFuns(i, u, p, U) 
            3 paraAndNodVect(Q) 
            4 solveTridiagonal(n, Q, U, P) 
                        需要注意的是list不能直接进行矩阵数组式的加减乘除，需要借助numpy
                        a = [1,2,3]    ---    b = 2.3 * np.array(a)    ----    c = b.tolist()
            5 curvePoint(n, p, U, P, u) 
            6 curvePlot() 
                        保存为txt文件     np.savetxt('文件名'，curve)


## 三、尝试解决曲面重构问题
          方法
                1 surfacePoint(n, p, U, m, q, V, P, u, v) 
                2 globalSurfInterp(Q) 
                            貌似BS在对于某些特别的三点采样的时候会产生很大误差(表现在BS曲线不通过中间的型值点)，而NURBS产生误差的几率和程度相对小很多---未找到原因(11/16)
            貌似python中的numpy.array()有一个坑：
          Q = [2.34222,5.143522,8.5254235,1.23]
    Q4 = [0 for i in range(4)] 
    Q4A = np.array(Q4)
    QA = np.array(Q)
    logging.info('Q4n is %s' % QA[0:4])
    Q4A = QA[0:4] # 第6行
    #Q4A[0:4] = QA # 第7行
    #Q4A[3] = 1
    logging.info('Q4n is %s' % Q4A)
            当如第6行这种赋值，不会出现什么奇怪的现象。
            但是，当如第7行这种赋值，则所有的内容都会被截断小数部分从而变为整数。
            在实际中应该尽量避免这种坑，实际上Q4A可以直接使用列表，不必化为数组(把等式右边的数组使用xx.tolist()重新化为列表即可)，这样既可以克服列表不能参与加减乘除的缺点，又避免了上述的数组的坑。

需要时再去尝试解决曲面上任意点的法向量问题
