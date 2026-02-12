# 等离子体色散函数（Plasma Dispersion Function）
本项目提供了一个基于 MATLAB 编写的函数 std_PDF，用于计算等离子体色散函数及其高阶导数的数值，并内置误差估计。该函数接受复自变量 zeta 进行高效、高精度的计算，适用于等离子体物理、空间物理及相关领域的数值模拟。

## 引用要求 🔖
若您在本研究或任何公开发表的成果中使用了此函数，必须引用以下文献：

Zhe Gao and Zhihao Su. "The Plasma Dispersion Function and its Numerical Calculation," [待发表：期刊/会议名称], 卷(期): 页码, 年份. DOI: [待填写]


## 函数说明
```matlab
function [Value, Error] = std_PDF(Zeta, Diff_Order, Error_Abs, Error_Rel)
```
### 输入参数
|参数名|类型|描述|
|:--------:|:--------:|:-------|
|Zeta|复数|标量|
|Diff_Order|整数|导数阶数，0 （默认）表示函数本身，1 表示一阶导数，依此类推。|
|Error_Abs|实数|绝对误差容限，控制数值计算的精度，默认1e-15。|
|Error_Rel|实数|相对误差容限，控制数值计算的精度，默认1e-15。|

### 输出参数
|参数名|类型|描述|
|:--------:|:--------:|:-------|
|Value|复数|函数值或导数值。|
|Error|实数|绝对误差估计值。|


## 快速开始
### 示例 1：计算单点处的函数值（0阶导数）
```matlab
[val, err] = std_PDF(1 + 0.5i, 0, 1e-10, 1e-10)
```

```text
val =
  -0.6077 + 0.6290i
err =
   9.8632e-12
```

### 示例 2：计算一组复数点的一阶导数
```matlab
[val, err] = std_PDF(3 + 3i, 1, 1e-12, 1e-12)

```

```text
val =
  -0.0045 - 0.0549i
err =
   2.1036e-13
```


## 依赖环境
MATLAB R2016b 或更高版本（推荐 R2020a 以上）

无需额外工具箱

## 作者信息
|项目|内容|
|:--------:|:-------|
|作者|Zhihao Su, Zhe Gao|
|机构|Department of Engineering Physics, Tsinghua University, Beijing 100084, CHINA|
|邮箱|gaozhe@tsinghua.edu.cn|
|创建日期|2022-06|
|最后更新|2023-09|



**致使用者**：请务必在"引用要求"部分填入您需要引用的具体文献。若您修改或扩展了此函数，建议注明原代码来源。
