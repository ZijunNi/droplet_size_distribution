<!-- <style>
  body {
    counter-reset: section; /* 初始化 section 计数器 */
  }
  h2 {
    counter-reset: subsection; /* 在每个 h2 中重置 subsection 计数器 */
  }
  h2::before {
    counter-increment: section; /* 递增 section 计数器 */
    content: counter(section) ". "; /* 显示 section 编号 */
  }
  h3::before {
    counter-increment: subsection; /* 递增 subsection 计数器 */
    content: counter(section) "." counter(subsection) ". "; /* 显示 subsection 编号 */
  }
</style> -->

# 流场后处理MATLAB代码使用说明

For English version of this readme file, please refer to ``readme_english.md``.

## 1. 简介
![alt text](<./doc/fig 1.png>)

 - 本流程新增第四步（Step 4）用于解决仅考虑最大持续时间时难以收敛的问题。通过计算所有剪切应力超过阈值的持续时间的概率密度函数（PDF），从中提取特定百分位数（对应图4中的"Given Ratio"，代码中对应`data.ratio`）对应的修正后时间$t_m$，从而建立新的阈值-持续时间关系（如Pic.5所示）。

## 2. 数据输入部分
以下是需要修改的变量详细说明：

### 2.1 输入数据文件
- **`data_set`**: 该变量指定输入数据文件名，文件需与`post_processing.mat`位于同一目录。文件应为`.mat`格式且包含以下变量：
  - `U`: 流向速度分量（使用摩擦速度$u_\tau$无量纲化，三维数组）
  - `V`: 展向速度分量（使用摩擦速度$u_\tau$无量纲化，三维数组）
  - `zpos_delta`: 壁面法向网格位置（使用通道半宽$\delta$无量纲化，一维数组）
  - `xpos_delta`: 流向网格位置（使用通道半宽$\delta$无量纲化，一维数组）

  示例：
  ```matlab
  data_set = "example_data.mat"; % 替换为实际数据文件名
  ```

### 2.2 雷诺数
- **`data.Reynolds_number`**: 该变量指定流动的雷诺数，**必须与所提供数据的雷诺数保持一致**。

  示例：
  ```matlab
  data.Reynolds_number = 3200; % 替换为实际雷诺数
  ```
### 2.3 二分法初始上下界
- **`left_bound`与`right_bound`**: 这两个变量分别代表程序使用二分法逼近颗粒临界直径时选择的初始上下界值。**需注意这两个变量使用通道半宽进行无量纲化**


## 3. 脚本运行流程
完成数据文件及雷诺数设置后，可直接运行脚本。脚本将自动执行以下步骤：

### 3.1 数据加载
脚本将加载输入数据文件并提取所需变量（`U`, `V`, `zpos_delta`, `xpos_delta`）。

### 3.2 创建数据存储目录
脚本将在当前工作目录下创建名为`data_N`的文件夹（其中`N`代表具体雷诺数`data.Reynolds_number`）用于存储结果和图表。若文件夹已存在则不会重复创建。

### 3.3 计算平均速度剖面
通过对流向速度分量（`U`）进行展向和流向平均，计算得到平均速度剖面，结果保存于`data.mean_U`。

生成平均速度剖面图并保存为PDF文件至`data_N`目录。

### 3.4 计算临界液滴尺寸
基于输入数据计算临界液滴尺寸，主要步骤包含：
- 定义液滴尺寸初始上下界（`left_bound`和`right_bound`）
- 遍历设定比例参数（`data.ratio`）计算临界值、物理持续时间、物理阈值和物理tau

生成各比例参数对应的临界液滴尺寸图表并保存为PDF文件至`data_N`目录。

### 3.5 结果保存
脚本将处理后的数据（`data_N`）保存为`.mat`文件存储于`data_N`目录，文件名基于雷诺数命名。

## 4. 输出文件
运行完成后将在`data_N`目录下生成：
- **平均速度剖面图**：PDF格式的平均速度分布图
- **临界液滴尺寸图表**：各比例参数对应的临界值PDF图表
- **处理后数据文件**：包含处理数据的`.mat`文件

## 5. 错误与警告提示
- **<font color=red>Initial Bound Error, [left_bound, right_bound] need to be adjusted.</font>**
该提示表明当前设置的上下界值对于当前数据集不适用，可通过调整**2.3节**所述的两个变量解决。实际使用中，颗粒临界直径通常略小于$y^+\approx 1$，可作为调整参考。
若上述方法无效，可修改文件`extract_2d_slice_x_interp.m`中第8行`if`语句后括号内的参数（默认值为1）。该参数控制是否使用对数插值，在低摩擦雷诺数且壁面网格解析度足够时（近壁网格点位于线性区），应将其设为0使用线性插值。注意修改后需恢复默认值1。

- **<font color=red>The mean velocity profile deviates from the logarithmic law. Please verify that the input velocity values are properly non-dimensionalized.</font>**
该提示表明输入数据的平均速度剖面与标准对数律存在显著偏差。**特别注意：输入速度场应使用摩擦速度$u_\tau$进行无量纲化，输入坐标应使用通道半宽$\delta$进行无量纲化**。本文所采用对数律为$u^+ = \ln(y^+)/\kappa + B$，其中$\kappa = 0.41$，$B = 5$。