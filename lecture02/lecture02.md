## <p align='center'>矩阵消元</p>

### 概述

高斯消元法（gaussian elimination）是一种方程组求解的方法，能够使用消元法要求系数矩阵A是一个**足够好**的矩阵，所谓的足够好其实就是后面要介绍的可逆矩阵。

在消元法中还有一个关键的概念称为**主元(pivot)**，主元的意思是用来旋转和变换的，在消元法中的意思是用主元来对其他行进行消元。因为消元的过程是将系数矩阵变为上三角矩阵(upper triangle)的过程，可以笼统将主元理解为对角线上的元素。

矩阵操作中遵循**左行右列**的思想，即在矩阵左边乘矩阵表示对该矩阵进行的是行变换，在矩阵的右边乘矩阵表示对该矩阵的列进行变换。

在此部分中以求解方程组
$$
\left\{\begin{array}{c}{x+2 y+z=2} \\ {3 x+8 y+z=12} \\ {4 y+z=2}\end{array}\right.
$$
为例。

### 消元过程

因为消元法是进行行变换，因此将方程组改写成行矩阵的形式，如下：
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {3} & {8} & {1} \\ {0} & {4} & {1}\end{array}\right]\left[\begin{array}{l}{x} \\ {y} \\ {z}\end{array}\right]=\left[\begin{array}{c}{2} \\ {12} \\ {2}\end{array}\right]
$$
消元法针对的是系数矩阵，其余需要处理的则是将每主元下面的元素通过行化简的方法变为0（因为要将系数矩阵变为上三角）。

首先是对第二行的处理，如下：
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {3} & {8} & {1} \\ {0} & {4} & {1}\end{array}\right] \underset{(2,1)}{\longrightarrow}\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {4} & {1}\end{array}\right]
$$
因为主元1下面的元素是3和0，因此只需要处理第二行即可，将3变为0，则需要将第二行减去第一行的3倍，从而得到右面的矩阵。

接着将对角线上的第二个元素看作主元，或者是将原矩阵中第一行第一列去掉得到新的矩阵，其实实质都是一样的，则接下来的变换根据主元及主元下面的元素确定要将第三行减去第二行的二倍。因此得到下一步的处理为：
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {4} & {1}\end{array}\right] \underset{(3,2)}{\longrightarrow}\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {0} & {5}\end{array}\right]
$$
此时得到的矩阵即是我们想要的上三角矩阵。

```
注：
在高斯消元法中主元不能为0，若是出现有0的情况，应该做换行处理：
首先看它的下一行对应位置是不是 0，如果不是，就将这两行位置互换，将非零数视为主元。
如果是，就再看下下行，以此类推。若其下面每一行都看到了，仍然没有非零数的话，那就意味着这个矩阵不可逆，消元法求出的解不唯一。
```

### 回带求解

回带求解应该是和消元法同步进行的，这里引入了**增广矩阵**的概念，即是指在系数矩阵的右边添上线性方程组等号右边的常数列得到的矩阵。因此针对方程组：
$$
\left\{\begin{array}{c}{x+2 y+z=2} \\ {3 x+8 y+z=12} \\ {4 y+z=2}\end{array}\right.
$$
对应的增广矩阵为：
$$
\left[\begin{array}{ccc}{1} & {2} & {1} & {2} \\ {3} & {8} & {1} & {12} \\ {0} & {4} & {1} & {2}\end{array}\right]
$$
之后的过程即是上述中高斯消元的过程，只是这个时候需要和等号右边的常数向量一起进行。过程如下：
$$
\left[\begin{array}{ccc}{1} & {2} & {1} & {2} \\ {3} & {8} & {1} & {12} \\ {0} & {4} & {1} & {2}\end{array}\right] \underset{(2,1)}{\longrightarrow}\left[\begin{array}{ccc}{1} & {2} & {1} &{2} \\ {0} & {2} & {-2} & {6} \\ {0} & {4} & {1} & {6}\end{array}\right] \underset{(3,2)}{\longrightarrow}\left[\begin{array}{cccc}{1} & {2} & {1} & {2} \\ {0} & {2} & {-2} & {6} \\ {0} & {0} & {5} & {-10}\end{array}\right] \rightarrow\left[\begin{array}{cccc}{1} & {2} & {1} & {2} \\ {0} & {2} & {-2} & {6} \\ {0} & {0} & {5} & {-10}\end{array}\right]
$$
将最后结果带入方程AX=b，则方程为：
$$
\left\{\begin{aligned} 1 x+2 y+z &=2 \\ 2 y-2 z &=6 \\ 5 z &=-10 \end{aligned}\right.
$$
可以很轻易求出解。

### 消元矩阵

#### 向量与矩阵乘法

在上一节中我们知道了
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {3} & {8} & {1} \\ {0} & {4} & {1}\end{array}\right]\left[\begin{array}{l}{x} \\ {y} \\ {z}\end{array}\right]
$$
可以看作是对A矩阵的三个列向量的线性组合，即上述矩阵可以看作是矩阵列的线性组合
$$
x\left[\begin{array}{l}{1} \\ {3} \\ {0}\end{array}\right]+y\left[\begin{array}{l}{2} \\ {8} \\ {4}\end{array}\right]+z\left[\begin{array}{l}{1} \\ {1} \\ {1}\end{array}\right]
$$


但是在高斯消元法中我们用到的是行向量的变换，那么行向量与矩阵相乘的结果是什么呢？还是以上述方程组为例，假设有行向量与矩阵相乘
$$
\left[\begin{array}{l}{x} & {y} & {z}\end{array}\right]\left[\begin{array}{lll}{1} & {2} & {1} \\ {3} & {8} & {1} \\ {0} & {4} & {1}\end{array}\right]
$$
结果可以看作是矩阵行的线性组合。
$$
x\left[\begin{array}{l}{1} & {2} & {1}\end{array}\right]+y\left[\begin{array}{l}{3} & {8} & {1}\end{array}\right]+z\left[\begin{array}{l}{0} & {4} & {1}\end{array}\right]
$$

### 消元矩阵

因为在消元法中我们主要就是针对的是矩阵的行变换，而所谓的**消元矩阵，就是将行变化的形式写成矩阵的乘法形式**。

首先根据上面的知识我们可以知道：
$$
\left[\begin{array}{l}{1} & {0} & {0}\end{array}\right]\left[\begin{array}{lll}{1} & {2} & {1} \\ {？} & {？} & {？} \\ {？} & {？} & {？}\end{array}\right]=\left[\begin{array}{1}{1}&{2}&{1}\end{array}\right]
$$
因为行向量中的每个元素控制矩阵对应的行数，若第一个元素为1其余元素为0，则最后的结果就是矩阵的第一行元素。

同理有
$$
\left[\begin{array}{l}{0} & {1} & {0}\end{array}\right]\left[\begin{array}{lll}{？} & {？} & {？}\\{1} & {2} & {1}\\ {？} & {？} & {？}\end{array}\right]=\left[\begin{array}{1}{1}&{2}&{1}\end{array}\right]
$$

$$
\left[\begin{array}{l}{0} & {0} & {1}\end{array}\right]\left[\begin{array}{lll}{？} & {？} & {？}\\{？} & {？} & {？}\\ {1} & {2} & {1}\end{array}\right]=\left[\begin{array}{1}{1}&{2}&{1}\end{array}\right]
$$

因此我们将$\left[\begin{array}{l}{1} & {0} & {0}\end{array}\right]$、$\left[\begin{array}{l}{0} & {1} & {0}\end{array}\right]$、$\left[\begin{array}{l}{0} & {0} & {1}\end{array}\right]$构成矩阵
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {0} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]
$$
称为单位矩阵记做$I$。**单位矩阵与任何矩阵相乘都得到该矩阵本身**（保证乘法有效，即矩阵的形状符合要求）。

接着我们使用矩阵的形式来阐述矩阵行变换的过程：
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {3} & {8} & {1} \\ {0} & {4} & {1}\end{array}\right] \underset{(2,1)}{\longrightarrow}\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {4} & {1}\end{array}\right]
$$
上面矩阵的变化过程是将第二行减去第一行的3倍，我们从单位矩阵入手，将单位矩阵的第二行减去第一行的三倍，得到矩阵
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {-3} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]
$$
而我们又由上面的知识可知$\left[\begin{array}{l}{0} & {1} & {0}\end{array}\right]$和$\left[\begin{array}{l}{0} & {0} & {1}\end{array}\right]$与原矩阵相乘得到矩阵对应的行数，因此将变化后的第二行元素$\left[\begin{array}{l}{-3} & {1} & {0}\end{array}\right]$与矩阵相乘得到矩阵的线性组合，为$\left[\begin{array}{l}{0} & {2} & {-2}\end{array}\right]$。因此得到该步骤对应的变换矩阵为
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {-3} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]
$$
我们称此类矩阵为**Elementary Matrix**，因为此步是变换矩阵的第二行第一个元素的，记为$E_{21}$。同理计算
$$
\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {4} & {1}\end{array}\right] \underset{(3,2)}{\longrightarrow}\left[\begin{array}{lll}{1} & {2} & {1} \\ {0} & {2} & {-2} \\ {0} & {0} & {5}\end{array}\right]
$$
此步骤对应的矩阵为
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {0} & {1} & {0}  \\ {0} & {-2} & {1} \end{array}\right]
$$
记为$E_{32}$。因此我们由系数矩阵A得到上三角矩阵U的过程可以表示为：
$$
E_{32}E_{21}A(系数矩阵)=U(上三角矩阵)
$$
使用结合律，先计算$E_{32}\cdot E_{21}$得到矩阵记做E，E即为此消元过程的消元矩阵。

注：**矩阵乘法不满足交换律但是满足结合律。**

因此求消元法对应的消元矩阵就是**从单位矩阵I入手，按照A每次变换的步骤对I进行相应的操作得到Elementary matrix E，最后累积得到E即可。**

#### 行变换和列变换

这里以2*2的矩阵为例简单阐述矩阵中的行变换和列变换。

假设2*2的矩阵为：
$$
\left[\begin{array}{cc}{a}  & {b}  \\ {c} & {d}  \end{array}\right]
$$
则交换两行的矩阵为
$$
\left[\begin{array}{cc}{0}  & {1}  \\ {1} & {0}  \end{array}\right]\left[\begin{array}{cc}{a}  & {b}  \\ {c} & {d}  \end{array}\right]=\left[\begin{array}{cc}{c}  & {d}  \\ {b} & {a}  \end{array}\right]
$$
交换两列的矩阵为：
$$
\left[\begin{array}{cc}{a}  & {b}  \\ {c} & {d}  \end{array}\right]\left[\begin{array}{cc}{0}  & {1}  \\ {1} & {0}  \end{array}\right]=\left[\begin{array}{cc}{b}  & {a}  \\ {d} & {c}  \end{array}\right]
$$
**左乘矩阵相当于对矩阵进行行变换，右乘表示对矩阵进行列变换。**

### 逆矩阵初探

在上面我们知道可以将一个矩阵A通过乘以一个矩阵E变为上三角矩阵U，那么如何将这个上三角矩阵U重新变为矩阵A呢？这个变化称为reversing，如上面提到的$E_{21}$，是将第二行减去第一行的3倍，那么我们将第二行再加上第一行的三倍即可以复原这一运算过程，即：
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {3} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {-3} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]=\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {0} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]
$$
此时矩阵
$$
\left[\begin{array}{ccc}{1} & {0} & {0}  \\ {3} & {1} & {0}  \\ {0} & {0} & {1} \end{array}\right]
$$
称为矩阵$E_{21}$的逆矩阵，记做$E_{21}^{-1}$。