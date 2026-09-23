# MODIS植被物候实验 · 学生离线包

面向植被遥感课程的独立实验：从单像元曲线和物候提取，逐步扩展到当年空间制图、多年变化与产品一致性比较。

## 下载

**[前往下载学生实验包](https://github.com/ZihangLou/modis-phenology-lab/releases/latest)**

在发布页的 **Assets** 中下载 **第八章_学生实验包.zip**（约1.03 GiB）。请下载这个完整附件，不要下载GitHub自动生成的“Source code (zip)”。也可下载 `.sha256` 文件核对完整性。

## 开始实验

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

## 数据与运行范围

MODIS覆盖2001—2025年，MCD12Q2参考产品覆盖2001—2024年。当年制图默认2024年、349×466格网。可选全区25年DL任务默认关闭（`RUN_FULL_SPATIAL=False`），不影响当年制图。

包内Notebook已清空运行输出，未包含教师机环境缓存和历史结果。已完成本机独立解压、离线安装、Jupyter内核及真实数据抽样检查；尚未在所有Windows版本和学生电脑上验证。

## 方法来源

两篇主要方法论文集中列在Notebook最后的参考文献单元。第三方依赖和字体的许可证随包保留。
