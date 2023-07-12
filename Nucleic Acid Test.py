# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 22:54:04 2023
核酸检测模拟程序
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
    


    
