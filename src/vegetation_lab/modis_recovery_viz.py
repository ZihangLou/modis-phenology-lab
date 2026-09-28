"""Plotting only for the third MODIS application notebook."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from .plot_style import configure_plot_fonts


def setup():
    configure_plot_fonts()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def finish(fig, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def annual_lines(table, columns, title, path, units):
    setup()
    fig, axes = plt.subplots(len(columns), 1, figsize=(11, 2.4*len(columns)), layout="constrained", squeeze=False)
    for ax, column, unit in zip(axes[:, 0], columns, units):
        ax.plot(table.year, table[column], "o-", color="#27805B", ms=4)
        ax.set(title=column, xlabel="年份", ylabel=unit)
        ax.grid(alpha=.2)
    fig.suptitle(title)
    finish(fig, path)


def recovery_maps(greenness, phenology, alpha, index, path):
    setup()
    panels = [(greenness["annual_mean"].slope, f"年均{index.upper()}线性趋势", f"{index.upper()}/年", False),
              (greenness["annual_max"].slope, f"年峰值{index.upper()}线性趋势", f"{index.upper()}/年", False)]
    slope, p = greenness["annual_mean"].slope, greenness["annual_mean"].p
    valid = np.isfinite(slope) & np.isfinite(p)
    classes = np.where(valid, np.where((p < alpha) & (slope > 0), 1,
                       np.where((p < alpha) & (slope < 0), -1, 0)), np.nan)
    panels.append((slope.copy(data=classes), "年均绿度方向：逐像元p未作多重检验校正", "", True))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), layout="constrained")
    for ax, (array, title, unit, categorical) in zip(axes.flat, panels):
        draw_map(ax, array, title, unit, categorical)
    for ax, metric in zip(axes[1], ["SOS", "EOS", "LOS"]):
        if phenology is None:
            ax.text(.5, .5, f"{metric}空间趋势\n待前一实验生成完整多年DL产品\n未用3×3验证或Zhang结果替代",
                    ha="center", va="center", transform=ax.transAxes, fontsize=11)
            ax.set_axis_off()
        else:
            draw_map(ax, phenology[metric].slope, metric+"空间趋势", "天/年", False)
    fig.suptitle("原下载矩形ROI：绿度变化不是生态工程因果证据")
    finish(fig, path)


def draw_map(ax, array, title, unit, categorical):
    lat, lon = array.lat.values, array.lon.values
    dx = abs(lon[1]-lon[0])/2 if len(lon) > 1 else .0025
    dy = abs(lat[1]-lat[0])/2 if len(lat) > 1 else .0025
    extent = [lon.min()-dx, lon.max()+dx, lat.min()-dy, lat.max()+dy]
    if categorical:
        image = ax.imshow(array, extent=extent, origin="upper", interpolation="nearest",
                          cmap=ListedColormap(["#be4248", "#d6d6d6", "#238257"]),
                          norm=BoundaryNorm([-1.5, -.5, .5, 1.5], 3))
        cb = ax.figure.colorbar(image, ax=ax, ticks=[-1, 0, 1], shrink=.8)
        cb.ax.set_yticklabels(["显著减少", "不显著", "显著增加"])
    else:
        finite = np.asarray(array)[np.isfinite(array)]
        limit = max(float(np.max(np.abs(finite))), 1e-9) if len(finite) else 1
        image = ax.imshow(array, extent=extent, origin="upper", interpolation="nearest",
                          cmap="RdYlGn", vmin=-limit, vmax=limit)
        ax.figure.colorbar(image, ax=ax, label=unit, shrink=.8)
    ax.set(title=title, xlabel="经度", ylabel="纬度")


def climate_scatter(analysis, pairs, path):
    setup()
    fig, axes = plt.subplots(1, len(pairs), figsize=(14, 4), layout="constrained", squeeze=False)
    for ax, (x, y) in zip(axes[0], pairs):
        data = analysis[[x, y]].dropna()
        ax.scatter(data[x], data[y], color="#27805B")
        ax.set(xlabel=x, ylabel=y, title=f"n={len(data)}")
        ax.grid(alpha=.2)
    fig.suptitle("气候相关：共同趋势、时间自相关与其他因素尚未控制")
    finish(fig, path)


def climate_correlation_maps(results, alpha, path):
    setup()
    computed = [(name, result) for name, result in results.items() if result is not None]
    if not computed:
        return
    fig, axes = plt.subplots(len(computed), 3, figsize=(14, 4.3*len(computed)),
                             layout="constrained", squeeze=False)
    for row, (relationship, result) in enumerate(computed):
        lat, lon = result.lat.values, result.lon.values
        dx, dy = abs(lon[1]-lon[0])/2, abs(lat[1]-lat[0])/2
        extent = [lon.min()-dx, lon.max()+dx, lat.min()-dy, lat.max()+dy]
        for ax, variable, title, cmap, low, high in zip(
                axes[row], ["r", "p", "n"], ["Pearson r", "双侧p（未校正）", "有效共同年份n"],
                ["RdBu_r", "viridis_r", "YlGnBu"], [-1, 0, 0], [1, 1, float(result.n.max())]):
            artist = ax.imshow(result[variable], extent=extent, origin="upper", interpolation="nearest",
                               cmap=cmap, vmin=low, vmax=max(high, 1))
            fig.colorbar(artist, ax=ax, shrink=.8, label=variable)
            ax.set(title=title, xlabel="经度", ylabel="纬度")
    fig.suptitle("气候双线性对齐MODIS格网：年降水与年均EVI")
    finish(fig, path)


def xr_like(reference, values):
    return reference.copy(data=np.asarray(values))


def stl_plot(table, path):
    setup()
    fig, axes = plt.subplots(4, 1, figsize=(12, 8), layout="constrained", sharex=True)
    for ax, name in zip(axes, table.columns):
        ax.plot(table.index, table[name], lw=1)
        ax.set_ylabel(name)
        ax.grid(alpha=.2)
    finish(fig, path)
