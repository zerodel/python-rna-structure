python-rna-structure
====

这里面包括了一些python代码，用于

- Rnasnp运行结果中抓取数据
- 建立模型矩阵
- 自动编写HYPHY 批处理文件


文件简述
====

- .pyc 文件是 同名.py 文件被调用后产生的。

- global_constants.py : 一些常量，包括常见的密码子等等

- gu_model.mdl ： Gu2010 文献中使用的模型矩阵
- rebuild_model.mdl  ： 用 gu_model.mdl 重建出来的GTR模型
- rna_structure.mdl : 新建立的模型文件。


pyRNAsnp 文件夹
----

里面包括了一些从.rnasnp 文件中抓取数据的函数。
- pyRNAsnp/pyRNAsnp.py ： 抓取.rnasnp 数据


pyHYPHY 文件夹
----

里面是和HYPHY打交道的一些脚本和文件

- pyHYPHY/BfHYPHY.py: 用于建立HYPHY简单的ML分析所需批处理文件（batch file）的一个类
- pyHYPHY/DataHYPHY.py: 用于建立HYPHY分析所需的序列文件（.input)
- pyHYPHY/MatrixHYPHY.py: 建立HYPHY矩阵
- pyHYPHY/ModelHYPHY.py : 之前用于建矩阵的一个类
- pyHYPHY/SpeciesSpecificTree.py ： 进化树的new-ick string 相关的一些函数。
— pyHYPHY/*.bf : batch file 的模板

test 文件夹
----
里面是每天为了一些功能所写的代码片段