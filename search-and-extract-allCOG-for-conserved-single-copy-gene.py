# coding:utf-8
#使用之前一定要备份一下*annotations文件和faa文件，处理过程中会修改faa文件
#以eggnog注释获得的annotation文件和蛋白faa文件为输入，提取所有COG，标准是：检测到多个基因就比较evalue,选择数值最小的一个。如果数值相同就取检测到的第一个。
import re
import os
from collections import Counter

'''
faaEdit函数说明：
功能：删除输入.faa文件中的 以# 开头，以0.xxx结尾的没用的数据
输入 path: 输入.faa文件的保存路径
'''


def faaEdit(path):

    suffix = ".faa"  # 输入文件的后缀名，这里为.faa

    listdir = os.listdir(path)  # 查找path中所有文件

    fileList = []  # 输入文件名的保存列表

    for name in listdir:

        if name.endswith(suffix):
            fileList.append(os.path.join(path, name))  # 如果文件以后缀 suffix结尾，就将文件全路径放到 fileList

    fileList.sort()  # 文件名排序

    # 开始循环
    for file in fileList:

        strOrg = ""

        # 读取 xxx.faa中的所有字符
        with open(file, "r") as f:

            lines = f.readlines()
            for line in lines:

                if ">" in line:
                    strOrg +=line.split(" ")[0]+"\n"
                else:
                    strOrg += line

        # 重新写入文件
        with open(file, "w") as f:
            f.write(strOrg)
        print(file + "  replace finished")  # 某文件处理完成


'''
getArcogGene函数说明：
功能：处理 .annotations文件
输入 path: 输入.annotations文件的保存路径
返回： 列表，列表中每个元素为一个字典，该字典保存了一个输入文件的内容，其中 key为基因名，value为arcog名

[{},{},{},...]---->{gene1:arcog1,gene2:arcog2,gene3:arcog3,gene4:arcog4,....}

'''


def getArcogGene(path):
    suffix = ".annotations"  # 以 .annotations结尾
    listdir = os.listdir(path)

    fileList = []

    for name in listdir:

        if name.endswith(suffix):
            fileList.append(os.path.join(path, name))

    fileList.sort()  # 文件名排序

    dataList = []  # 最终返回的结果列表

    # 开始循环读入文件
    for file in fileList:

        dictData = {}  # 存放单个文件的结果
        # tempDelList = []  # 临时存放所有arcog的列表

        # 打开文件
        with open(file, "r") as f:

            lines = f.readlines()
            for line in lines:

                temp = line.strip("\n").split("\t")  # 以 tab键 \t分割

                # 删除空格
                while ('' in temp):
                    temp.remove('')

                #跳过开头
                if (len(temp) > 2):

                    #跳过表头
                    if ("#query" in temp):
                        continue

                    key1 = ">" + temp[0]  # 第一个[0]就是基因名，给基因名前面加上 >

                    valueTempList = temp[4].split(",")  # 第5列 就是 arcog所在的列， 该列用","隔开，
                    
                    cog_val=float(temp[2])

                    value1 = ""
                    for i in valueTempList:

                        # if "arCOG" in i:
                        if "arCOG" not in i and "COG" in i:
                            value1 = i.split("@")[0]  # # 如果包含arcog，将arcog用"@"隔开，只取 "@" 前面的

                            resValue = re.sub('"', '', value1)  # 有些 arcog 开头带有 "，需要去除掉
                            dictData[key1] = (resValue,cog_val)  # key就是 基因名，value就是 arcog，存到字典中
                            # tempDelList.append(resValue)  # 同时把所有的 arcog都放到 tempDelList列表中，以作后续处理
                            break


        saveList=[]

        for key,val in dictData.items():

            flag=0
            for d in saveList:
                if val[0]  in d:
                    flag = 1

                    if val[1]>d[1]:
                        break
                    else:
                        saveList.remove(d)
                        saveList.append(val)
                        break
            if flag==0:
                saveList.append(val)

        dictData1={}

        dict_tmp={}

        for key, val in dictData.items():

            flag=0
            for _,val1 in dict_tmp.items():
                if val1==val:
                    flag=1
                    break

            if flag==0:
                dict_tmp[key]=val



        for key,val in dict_tmp.items():

            if val in saveList:

                dictData1[key]=val[0]

        # delDict = dict(Counter(tempDelList))  # Counter  可以查出 列表中，每个元素出现的次数， 并返回如 ：[(a,1),(b,2)....]的列表，将其转化为字典
        #
        # delList = []  # 存放出现超过两次的arcog名称，用于删除
        #
        # for key, val in delDict.items():
        #     if val > 1:
        #         delList.append(key)  # 如果出现 arcog在 tempDellist中出现超过两次，就将该arcog放到 待删除列表中
        #
        # for key, val in list(dictData.items()):
        #
        #     if val[0] in delList:
        #
        #         dictData.pop(key)  # 循环该文件的字典，如果arcog出现在delDict列表中，就从字典中删除这一条信息

        dataList.append(dictData1)  # 处理完的的字典放到返回的总列表中

    return dataList  # 返回列表


'''
getFaaData函数说明：
功能：处理 .faa文件
输入 path: 输入.faa文件的保存路径
返回： 列表，列表中每个元素为一个字典，该字典保存了一个输入文件的内容，其中 key为基因名，value为基因序列

[{},{},{},...]---->{gene1:序列1,gene2:序列2,gene3:序列3,gene4:序列4,....}

'''


def getFaaData(path):


    faaEdit(path)  # 先调用函数，删除不需要的字符串

    comp = re.compile("(^>.*_\d+)([A-Za-z].*)")  # 用于查找基因名，名称以 >开头 ，以 _xxx结尾,xxx为数字  如：>liuronghua_222

    suffix = ".faa"
    listdir = os.listdir(path)

    fileList = []

    for name in listdir:

        if name.endswith(suffix):
            fileList.append(os.path.join(path, name))

    fileList.sort()

    dataList = []

    for file in fileList:

        dictData = {}  # 存放单个文件的结果

        strTemp = ""

        with open(file, "r") as f:

            lines = f.readlines()
            for line in lines:
                strTemp += line.strip("\n")  # 读取每一行，删除 换行符 \n

        tempDataList = [">" + i for i in strTemp.split(">") if i]  # 以 >分割，得到每个基因，但split分割的同时会删除掉 > 所有，再给每个基因开头补上 >

        for data in tempDataList:

            res = comp.search(data)  # 查找每个基因的名字

            if (res.groups()):
                key = res.groups()[0]  # 如果找到了，那么第一个[0]就是基因名
                val = res.groups()[1].replace(" ", "")  # 第二个[1]就是基因序列，  同时去掉基因序列里面的空格，可能存在空格
                dictData[key] = val

        dataList.append(dictData)

    return dataList


'''
getArcogFile函数说明：

功能：得到 arcog_xxx.faa文件

输入 arcogPath: 输入.annotations文件的保存路径
输入 faaPath: 输入.faa文件的保存路径
输入 savePath: arcog_xxx.faa文件的保存路径

'''


def getArcogFile(arcogPath, faaPath, savePath):
    arcogDataList = getArcogGene(arcogPath)  # 调用函数，得到 annotations文件数据
    faaDataList = getFaaData(faaPath)  # 调用函数，得到 faa文件数据

    # 对每个文件信息开始循环
    for i in range(len(arcogDataList)):

        arcogDict = arcogDataList[i]  # 单个annotations文件的字典
        faaDict = faaDataList[i]  # 单个faa文件的字典

        for key, val in arcogDict.items():

            listdir = os.listdir(savePath)  # 先判断结果保存路径中有没有要输出到的 arcog_xxx文件
            if val + ".faa" in listdir:
                # 如果有，那么就以追加的方式打开
                with open(os.path.join(savePath, val + ".faa"), "a") as f:
                    if key in faaDict.keys():
                        f.write(key + "\n" + faaDict[key] + "\n")  # 判断该annotations文件中的该arcog对应基因名是否出现在该输入的faa文件中
                        # 如果存在，就将 对应的基因序列写入到该arcog文件中
                        # 写入的格式为  基因名 + \t + 序列 +\n
                    else:
                        f.write(key + "\n")  # 如果不存在，就只写 基因名 + \t  +\n
            else:
                # 如果没有该 arcog_xxx.faa文件，那么就以写入的方式打开，相当于新建，其余跟上面操作一致
                with open(os.path.join(savePath, val + ".faa"), "w") as f:
                    if key in faaDict.keys():
                        f.write(key + "\n" + faaDict[key] + "\n")
                    else:

                        f.write(key +"\n")

        print(i + 1)  # 完成一个 输入文件，就输出一个序号，方便确定进度，当然也可以给你搞个进度条，但是觉得没必要


if __name__ == '__main__':
    faaPath = "/lustre/home/zhangxiaohua/liuronghua/MT5-0904/20230713-amended-paper/making-new-tree/COG_markergene/231-faa-168-41-22/"  # 输入 faa保存路径
    arcogPath = "/lustre/home/zhangxiaohua/liuronghua/MT5-0904/20230713-amended-paper/making-new-tree/COG_markergene/231-emapper-annotations/"  # 输入annoataions保存路径
    savePath = "/lustre/home/zhangxiaohua/liuronghua/MT5-0904/20230713-amended-paper/making-new-tree/COG_markergene/COG_files/"  # 输出结果保存路径

    getArcogFile(arcogPath, faaPath, savePath)  # just run it ！！！！！

