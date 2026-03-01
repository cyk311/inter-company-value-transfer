# 说明文档
## 核心计算代码
| 核心计算部分分为五个步骤，拆解放在三个文件中
### Step_1.m （Matlab代码）
  - 预处理代码，实现步骤1
    - 用来将人大IO表转化为面板数据，方便使用stata处理
### Step_23.do（Stata代码）
  - 步骤2
    - 为每个(i,j)对应的技术时间序列做DFT，计算主周期c_ij
  - 步骤3
    - 分别带入3个线性回归模型 + 模型选择 + 导出 Ahat_panel.csv
### Step_45.m（Matlab代码）
- 步骤4
  - 计算价值和生产价格，并使用A_hat矩阵构造p_ij
- 步骤5
  - 计算UE_ij矩阵
  
## 数据文件
| 数据文件共四个，分别为计算部分代码的输入和输出
### 原始数据：IO(1981-2018).xlsx
### Step_1.m输出文件：A_panel.csv
### Step_23.m输出文件：Ahat_panel.csv
### Step45.m输出文件：value_transfer_long.csv

## 绘图代码
| 除了DFT示例图之外，均基于value_transfer_long.csv
### DFT示意图：DFT_exg
### 不平衡交换程度热力图：Heatmap
- Heatmap.m 绘制了全部时间的热力图，存放在ue_heatmaps文件夹中
- Heatmap_2.m 挑选了九张间隔四年的热力图，用于Slides展示
### 显著性热力图：Significance
- Significance_1.m 计算了显著性矩阵
- Significance_2.m 绘制热力图
### 上下游相对位置对价值转移程度的影响：Up-downstream
- Up_downstream.m 计算了指标数据，存放在Nt_summary.csv中，并绘制了指标的时序图