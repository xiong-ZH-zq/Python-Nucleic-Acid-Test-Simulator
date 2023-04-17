# Python-Nucleic-Acid-Test-Simulator
本文中的数学建模问题来源于NKU的数学建模第二次实战演练，由于本次是我来进行程序的编写，故将代码与笔记记录在这里。

# 问题提要
现有800万市民报名参与核酸检测，如果对每人逐一进行检测，所需时间和检测能力都超过现实情况，所以拟采用混样检测(group testing)方式进行。先考虑混样规模为$k$人份一组（一份的咽拭子样本不宜过多，限定$k\leq10$）。$k$个人的样本混在一起进行检测，如果是阴性，那就全部通过，如果是阳性，可以将该组再进一步分组，或者逐一进行检测。假设人群的阳性率为$\alpha$  (不超过10\%)，每个人最多采集$n$ 个样本（限定$n$不超过3，太多，对于被测试者不友好，也不宜保存和管理）。在满足上述条件的情况下，请设计一个方案（可以是理论证明和公式，也可以是模拟程序），提供给核酸检测团队，帮助核酸检测人员尽量高效的完成大规模（假设固定为800万）核酸检测。

# 问题分析
根据问题的实际情况，考虑利用Phython程序来进行模拟。为了计算结果的速度，我们考虑将样本压缩到80000人，当然接下来的思路实际上可以套用到任意总数上，只不过到了800000程序运行时间就肉眼可见的慢了。
下面我们考虑进行多次分组。
## 第一次分组
根据问题条件，首次混样分组以$k$人份为一组，且限定$k\leq10$.则此时可以知道第一次分组混检次数为$N_1=\lceil \frac{8\times10^7}{k} \rceil$.第一次分组混检次数与分组在$k$确定时已经确定。

## 第二次分组
下面考虑第二次分组，第一次混检结果为阴性的组此时已无需考虑，结果为阳性的组仍需进一步检测以确定阳性病例。为使问题研究简便，下面我们考虑单独一个人数为$k$的阳性小组在不同分组情形下的检测次数。

下面我们将第二次分组的问题简化抽象为以下问题：

1. 对于4人小组，我们可将该小组拆分为$1+1+1+1,1+2+1,2+2,1+3$的正整数和的形式，针对本问题不考虑其本身，对于任意的正整数$k$,如何罗列出$k$人的小组拆分成多个正整数和的方案？

2. 我们将1.中拆分之后的小组称为**细分小组**，如何计算或模拟不同细分方案下核酸检测所需次数？

针对以上问题，考虑到实际细分小组中的阳性病例的随机性，考虑利用Python程序模拟，通过产生随机样本，经过多次模拟试验，估测出各方案的可能取值范围，从而帮助我们做出决策。

## 第三次核酸检测
在经历以上步骤后，第三次核酸检测由于$n\leq3$的限制只能进行每人的单检（在仍然不能确定阳性与阴性的前提下）。由此我们仅需考虑两次分组问题。

# 模型假设
为方便有效地对问题进行研究, 我们不妨作如下假设:

1. 假设在新冠病毒核酸检测中, 所用的试剂判断准确率为 100%.

2. 假设在新冠病毒核酸检测中, 所用的试剂同样可以准确地检测出混合样本, 且判断准确率为 100%. 在混合样本均为阴性时, 显阴性特性, 否则显阳性特性。

3. 假设每个人的取样数足够多, 在样本分组后, 各组之间相互独立, 而且样本分组或稀释不会影响最终检测的结果。

# 程序建构
## 生成样本
下面是生成总样本的代码：
```python
def generate_mass(n):
    mass = []
    #alpha = random.uniform(0.01,0.1)
    alpha = 0.05
    for i in range(1,n+1):
        if i <= n*alpha:
            mass.append(1)
        else:
            mass.append(0)
    random.shuffle(mass)    # 洗牌模拟真实情形
    return mass
```
这个函数在传入总人数$n$后生成并返回一个长度为$n$的列表，其中内容为整数0，1.0代表实际阴性人员，1表示实际阳性病例。我们利用`random.shuffle()`方法将生成的列表洗牌打乱以模拟真实情形。如果想要让$\alpha$也随机生成，可以选用注释掉的那一句代码。

## 第一次分组
在生成样本之后我们考虑将样本列表进行第一次分组，即生成一个包含$N_1$个长度为$k$的列表：
```python
def group1(mass,k,n):
    global N,N1
    N=0
    groups1 = []
    group = []
    for i in range(1,n+1,k):
        for j in range(i,i+k):
            if j<=n:
                group.append(mass[j-1])
            else:
                pass
        groups1.append(group)
        N +=1
        group = []
    #print(N)
    N1=N
    return groups1
```
注意这里的$N$与$N1$为程序开始时定义的全局变量，用于存储总检测次数和第一次分组混检的检测次数。返回的`groups1`即为我们需要的列表。

## 第二次分组
### 分组方案生成
第二次分组由于既可单检，也可混检，因此我们需要考虑$k$的加法分解。
例如$k=4$，我们用列表表示可能的正整数分组就有$[1,1,1,1],[1,1,2],[2,2],[1,3]$，我们不考虑$[4]$本身，因为这对第二次检测是没有意义的。
因此我们需要将任意$k$的方案输出：
```python
# 考虑找出k的所有正整数加法分解，且每个数最大为k-1

def find_all_partitions(k):
    # 初始化结果列表
    partitions = []
    # 递归函数，用于生成所有的分解
    def find_partitions(n, lst):
        if n == 0:
            # 如果n为0，说明已经找到了一组分解
            lst.sort() # 对分解进行排序
            if lst not in partitions: # 如果分解不在结果列表中，就添加到结果列表
                partitions.append(lst)
        else:
            # 否则从1开始尝试，递归查找所有的分解
            for i in range(1, n + 1):
                find_partitions(n - i, lst + [i])

    # 调用递归函数，生成所有的分解
    find_partitions(k, [])
    # 返回结果列表
    return partitions

```
利用递归函数我们可以将任意正整数$k$的加法分解求出。注意代码中的排序作用是除去仅仅只是顺序不同的方案。我们利用循环从$k=2$至$k=10$可以输出我们需要的全部方案。

### 分组方案实施
```python
def group2(groups1):
    
        global N,solutions,N1
        N=N1
        positives = []
        positives2=[]   #需要第三次检测的组
        groups2 = []
        # 先踢掉全阴性的组
        for group in groups1:
            for i in range(0,len(group)):
                if group[i] == 1:
                    positives.append(group)
        #print(positives)

    # 再考虑第二次分组检测，从所有可能分组中进行循环
    
        for positive in positives: 
            num=0
            for i in range(0,len(solution)):
                    # 利用切片取出细分小组
                    groups2.append(positive[num:num+solution[i]])
                    N+=1
                    num += solution[i]
        print('在{0}的情况下第二次分组为'.format(solution))
        #print(groups2)
        for group in groups2:
            for i in range(0,len(group)):
                if group[i] == 1 and len(group)>=2:
                    positives2.append(group)
        #print('这里面需要参加第三次的组为：')
        #print(positives2)
        for i in range(0,len(positives2)):
            for j in range(0,len(positives2[i])):
                N+=1
        print(str(N))
```
上面的代码将会输出对应分组方案的总检测次数，如果需要多次试验可以将第一个`print`注释掉（在明确结果对应的分组之后），从而方便我们将数据复制到表格当中进行分析。
# 结果分析
在运行程序20次试验之后我们将结果汇总并将检测次数最少的八个方案列表如下（总数为$80000$，$\alpha=0.05$）：
| 分组方案            | 平均次数  | 标准差         |
|-----------------|-------|-------------|
| [2, 2, 2]       | 35339 | 78.56818695 |
| [1, 2, 2]       | 35661 | 80.75215167 |
| [1, 1, 1, 2]    | 35833 | 64.58018272 |
| [1, 1, 1, 1]    | 36000 | 0           |
| [1, 1, 1, 1, 1] | 36000 | 0           |
| [1, 1, 2, 2]    | 36007 | 84.90959899 |
| [2, 2, 3]       | 36038 | 114.1926442 |
| [1, 2, 3]       | 36184 | 106.2002236 |

根据表中的结果来看分组方案中平均次数较少的方案有如下特点：


+ 第二次分组中少量均分方案检测次数相对较少，其中单检和最大组为2人的分组方案尤其较少。

+ 第一次分组时$k=5,6$的情况检测次数相对较少。

+ 第二次分组单检的情况次数稳定，而混检次数可能根据实际情况出现接近上百次的浮动。

我们因此可以得到如下的结论：
+ 第一次分组令$k=5,6$对于减少检测次数是有效的，考虑到整除问题选用$k=5$较为合理且便于操作。

+ 最具有实操性且检测次数在较少的同时还能保持稳定的方案是第二次分组时直接采用单检方案。

综上我们可以认为$k=5$且第二次对所有阳性小组进行单检的方案是最适宜于实际操作的。

# 注意事项
模型还能从以下角度进行改进：
+ 在试验时可以在$\alpha = 0.01,0.02,\cdots,0.1$的条件下分别进行20次试验，可以找到$\alpha$不同时分别的最优方案。

+ 代码中仍有$\frac{n}{k}$不整除时带来的误差，这源于将分组方案直接套用到不足$k$人小组造成的误差。

# 完整代码
完整代码如下：
```python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 22:54:04 2023
核酸检测模拟程序

@author: xzqbear
"""
import random
N,N1=0,0
solution = []
solutions = []
# 生成总人数列表
def generate_mass(n):
    mass = []
    #alpha = random.uniform(0.01,0.1)
    alpha = 0.05
    for i in range(1,n+1):
        if i <= n*alpha:
            mass.append(1)
        else:
            mass.append(0)
    random.shuffle(mass)    # 洗牌模拟真实情形
    return mass

# 第一次分组
def group1(mass,k,n):
    global N,N1
    N=0
    groups1 = []
    group = []
    for i in range(1,n+1,k):
        for j in range(i,i+k):
            if j<=n:
                group.append(mass[j-1])
            else:
                pass
        groups1.append(group)
        N +=1
        group = []
    #print(N)
    N1=N
    return groups1

'''
在进行第二次分组之前我们需要考虑k人小组可能的分组情形
我们考虑先寻找k这个正整数可能的加法分解情况，例如对于4，可能的分解有：
1+1+1+1，1+2+1，1+3，2+2，4
根据问题，我们可以认为其本身作为加法分解是没有意义的，所以上述情况
符合提议的分解需要去掉4.
下面考虑找出所有可能的分解情况
'''
# 考虑找出k的所有正整数加法分解，且每个数最大为k-1

def find_all_partitions(k):
    # 初始化结果列表
    partitions = []
    
    # 递归函数，用于生成所有的分解
    def find_partitions(n, lst):
        if n == 0:
            # 如果n为0，说明已经找到了一组分解
            lst.sort() # 对分解进行排序
            if lst not in partitions: # 如果分解不在结果列表中，就添加到结果列表
                partitions.append(lst)
        else:
            # 否则从1开始尝试，递归查找所有的分解
            for i in range(1, n + 1):
                find_partitions(n - i, lst + [i])

    # 调用递归函数，生成所有的分解
    find_partitions(k, [])
    
    # 返回结果列表
    return partitions


# 第二次分组

def group2(groups1):
    
        global N,solutions,N1
        N=N1
        positives = []
        positives2=[]   #需要第三次检测的组
        groups2 = []
        # 先踢掉全阴性的组
        for group in groups1:
            for i in range(0,len(group)):
                if group[i] == 1:
                    positives.append(group)
        #print(positives)
    

        
    # 再考虑第二次分组检测，从所有可能分组中进行循环
    
        for positive in positives: 
            num=0
            for i in range(0,len(solution)):
                    # 利用切片取出细分小组
                    groups2.append(positive[num:num+solution[i]])
                    N+=1
                    num += solution[i]
        #print('在{0}的情况下第二次分组为'.format(solution))
        #print(groups2)
        for group in groups2:
            for i in range(0,len(group)):
                if group[i] == 1 and len(group)>=2:
                    positives2.append(group)
        #print('这里面需要参加第三次的组为：')
        #print(positives2)
        for i in range(0,len(positives2)):
            for j in range(0,len(positives2[i])):
                N+=1
        print(str(N))
    

if __name__ == '__main__':

    n = 80000    # 人数
    mass = generate_mass(n)
    #print(mass)
    for k in range(2,11):
        groups1 = group1(mass, k, n)
        #print(groups1)
        solutions = find_all_partitions(k)
        del solutions[-1]
        for solution in solutions:
            group2(groups1)
```
