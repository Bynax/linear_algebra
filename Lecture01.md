## <p align='center'>方程组的几何解释</p>

### 概述

本节的主要目的是从**行图像(row picture)**和**列图像(column picture)**的角度来理解方程。

### 方程组的几何解释

在本节中以方程
$$
\left\{\begin{array}{ccc}{2 x-y} & {=} & {0} \\ {-x+2 y} & {=} & {3}\end{array}\right.
$$
为例，通过行图像和列图像来阐释此方程组的含义。

#### 行图像角度

首先我们可以按照行将方程写成如下的矩阵形式：
$$
\left[\begin{array}{cc}{2} & {-1} \\ {-1} & {2}\end{array}\right]\left[\begin{array}{l}{x} \\ {y}\end{array}\right]=\left[\begin{array}{l}{0} \\ {3}\end{array}\right]
$$
约定将前面的2X2矩阵记为A，称作系数矩阵。将后面的系数向量记为**X**，称作系数向量，最后的结果向量记为b。

```
系数矩阵(A): 将方程组系数按行提取出来，构造完成的一个矩阵。
未知向量(x): 将方程组的未知数提取出来，按列构成一个向量。
向量(b): 将等号右侧结果按列提取，构成一个向量。
```

**行图像**则是将改方程组中的方程绘制出来，交点即是改方程组的解。

#### 列图像角度

首先同行图像方法相同，我们将上面的方程组按照列提取出来构成如下矩阵：
$$
x\left[\begin{array}{c}{2} \\ {-1}\end{array}\right]+y\left[\begin{array}{c}{-1} \\ {2}\end{array}\right]=\left[\begin{array}{l}{0} \\ {3}\end{array}\right]
$$
此时我们就将方程组的问题转化成寻找一组合适的参数，将向量$\left[\begin{array}{c}{2} \\ {-1}\end{array}\right]$和$\left[\begin{array}{c}{-1} \\ {2}\end{array}\right]$线性组合构成$\left[\begin{array}{l}{0} \\ {3}\end{array}\right]$。因此$AX$可以看作是A各列的线性组合。

