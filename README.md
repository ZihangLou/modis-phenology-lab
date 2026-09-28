# MODIS植被物候实验 · 学生离线包



面向植被遥感课程的独立实验：从单像元曲线和物候提取，逐步扩展到当年空间制图、多年变化与产品一致性比较。



## Windows 下载

**[前往下载学生实验包](https://github.com/ZihangLou/modis-phenology-lab/releases/latest)**

在发布页的 **Assets** 中下载 **modis-phenology-student.zip**（约1.03 GiB）。请下载这个完整附件，不要下载GitHub自动生成的“Source code (zip)”。

## Mac 离线包

要求macOS 14或更新版本，请按芯片选择：
- [Apple M芯片版](https://github.com/ZihangLou/modis-phenology-lab/releases/download/v1.1.1-macos/modis-phenology-macos-arm64-v1.1.1.zip)（约945 MiB，Numba后端）。

完整解压后双击 **02_启动实验.command**，首次自动完成离线安装和自检并打开Notebook。无需另装Python、Conda或Homebrew。

## Windows 开始实验

1. 将ZIP完整解压到Windows 64位电脑的可写目录，建议额外预留至少5 GB空间。

2. 双击 **01_首次准备.cmd**，等待离线安装和环境检查通过。

3. 双击 **02_启动实验.cmd**，在JupyterLab打开 `notebooks/modis_phenology_methods.ipynb`，从第一格依次运行。

4. 环境异常时运行 **03_环境检查.cmd**。实验结果写入包内 `outputs/`。



随包提供Python及锁定的离线依赖，无需另装Python、Conda或登录GEE。初始化后不要移动目录；换位置时请重新解压原始ZIP并准备环境。

## 实验内容

- MODIS植被指数读取与QA质量控制。

- DL拟合、参数物候日期与50%动态阈值。

- 迭代SG、Whittaker和DL的单像元对比。

- 分段Logistic几何特征提取。

- 当年SG50%的SOS、POS、EOS制图，独立计算LOS，并比较SG50%与DL参数法的SOS/EOS差值。

- 固定像元多年变化、可选全区多年制图与趋势。

- MCD12Q2参考产品一致性比较。

## 实验03、04补充资料

已在仓库按原学生包目录结构提供 `notebooks/`、`data/`、`src/`、`config/` 和 `docs/`。

- [下载补充ZIP](https://github.com/ZihangLou/modis-phenology-lab/releases/download/v1.0.1-supplement/experiments-03-04-supplement.zip)
- [合并到已有学生包的操作说明](docs/实验03-04_补充资料使用说明.md)

下载后将五个同名文件夹**合并**到原Windows学生包根目录，保留已有文件。无需重装环境。补充包不含原MODIS数据和运行环境，需配合Windows v1.0.0使用。实验04的完整物候空间分析需要此前生成的多年DL结果；缺失时显示待计算。
