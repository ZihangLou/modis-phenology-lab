"""HLS v3 IO and figures; calculations use the notebook's supplied functions."""
from pathlib import Path
from datetime import datetime
import hashlib
import inspect
import json
import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from rasterio.features import bounds
from affine import Affine
import matplotlib.pyplot as plt
from . import hls_phenology_io as hio
from .plot_style import configure_plot_fonts

LABELS = {"DB_FOREST": "落叶阔叶林", "EB_FOREST": "常绿阔叶林", "STEPPE": "草地",
          "NC_WHEAT": "华北冬小麦春季", "NC_MAIZE": "华北夏玉米",
          "NE_MAIZE": "东北玉米", "NE_RICE": "东北水稻", "NE_SOY": "东北大豆", "JS_RICE": "江苏水稻"}
METRICS = ["SOS", "POS", "EOS", "LOS"]


def new_run(root):
    path = root / "outputs/hls_teaching_v3" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    (path / "figures").mkdir(parents=True)
    return path


def load_region(root, roi, full=False):
    config = json.loads((root / "config/phenology_display_samples.json").read_text(encoding="utf8"))
    selected = config.get(roi)
    center = (selected["lon"], selected["lat"]) if selected else None
    source = "existing_display_config_NOT_verified_crop" if selected else "downloaded_ROI_center"
    if center is not None:
        folder = root / "data/hls" / ("Ch08_HLS_NDVI_v1_"+roi)
        geometry = json.loads(sorted(folder.glob("*roi.geojson"))[-1].read_text(encoding="utf8"))["features"][0]["geometry"]
        west, south, east, north = bounds(geometry)
        if not (west <= center[0] <= east and south <= center[1] <= north):
            print(f"{roi}: existing display coordinate {center} is outside the current downloaded ROI; using its recorded center, no pixel search.")
            center = None
            source = "downloaded_ROI_center_old_config_outside_current_ROI"
    ds = hio.load_patch(root / "data/hls", roi,
                        patch_size=None if full else 1, center=center)
    ds.attrs["sample_source"] = source
    return ds


def point_record(ds, qa_function):
    row, col = int(ds.attrs["sample_row"]), int(ds.attrs["sample_col"])
    point = ds.isel(y=row, x=col)
    valid = qa_function(point.ndvi_raw.values, point.fmask.values)
    return pd.DataFrame({"date": point.time.values, "raw": point.ndvi_raw.values,
                         "fmask": point.fmask.values, "qc": np.where(valid, point.ndvi_raw.values, np.nan)})


def dates_row(name, values, start):
    origin = pd.Timestamp(start)
    row = {"method": name}
    for metric, value in zip(METRICS, values):
        if metric == "LOS":
            row["LOS_days"] = value
        else:
            date = origin + pd.to_timedelta(value, unit="D") if np.isfinite(value) else pd.NaT
            row[metric] = date
            row[metric+"_doy"] = (date-pd.Timestamp(date.year, 1, 1)).total_seconds()/86400+1 if pd.notna(date) else np.nan
    return row


def save_configuration(root, folder, settings, functions):
    script = root / "src/vegetation_lab/chen_sg_filter.py"
    metadata = {**settings, "SG_source_task": "01a0b563-b2b0-7fa2-94be-fe009821819e",
                "SG_script": str(script.relative_to(root)), "SG_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
                "algorithms": {f.__name__: inspect.getsource(f) for f in functions}}
    (folder / "configuration.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf8")


def _figure():
    configure_plot_fonts()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def _finish(fig, destination):
    fig.savefig(destination, dpi=135, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def show_series(points, results, title, path, windows=None):
    _figure()
    fig, ax = plt.subplots(figsize=(12, 4.8), layout="constrained")
    for key, point in points.items():
        label = LABELS.get(key, key)
        if key in results:
            grid, curve = results[key]
            point = point[(point.date >= grid[0]) & (point.date < grid[-1]+pd.Timedelta(days=1))]
        ax.scatter(point.date, point.qc, s=10, alpha=.45, label=label+"观测")
        if key in results:
            grid, curve = results[key]
            ax.plot(grid, curve, lw=1.6, label=label+" Chen SG")
    for name, (start, end) in (windows or {}).items():
        ax.axvspan(pd.Timestamp(start), pd.Timestamp(end), alpha=.08, label=name)
    ax.set(title=title, ylabel="NDVI（不分别归一化）", xlabel="日期")
    ax.legend(frameon=False, ncol=3, fontsize=9)
    ax.grid(alpha=.2)
    fig.autofmt_xdate()
    _finish(fig, path)


def show_methods(point, grid, curves, path):
    _figure()
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
    ax.scatter(point.date, point.raw, s=12, color=".75", label="原始")
    ax.scatter(point.date, point.qc, s=14, color="black", label="QA后真实观测")
    for name, curve in curves.items():
        ax.plot(grid, curve, label=name, lw=1.5)
    ax.set(xlabel="日期", ylabel="NDVI", title="落叶阔叶林：同一观测，不同曲线方法")
    ax.legend(frameon=False, ncol=3)
    ax.grid(alpha=.2)
    fig.autofmt_xdate()
    _finish(fig, path)


def show_events(grid, curve, events, path, title, ratio=None):
    _figure()
    fig, ax = plt.subplots(figsize=(12, 4.5), layout="constrained")
    ax.plot(grid, curve, color="#27805B", label="拟合/重建曲线")
    for name, offsets in events.items():
        for metric, offset in zip(METRICS[:3], offsets[:3]):
            if np.isfinite(offset):
                date = grid[0] + pd.to_timedelta(offset, unit="D")
                ax.axvline(date, lw=1, ls="--", label=f"{name} {metric}: {date:%m-%d}")
    if ratio is not None and np.isfinite(curve).any():
        peak = int(np.nanargmax(curve))
        for base in [np.nanmin(curve[:peak+1]), np.nanmin(curve[peak:])]:
            ax.axhline(base+ratio*(curve[peak]-base), color=".5", ls=":")
    ax.set(title=title, xlabel="日期", ylabel="NDVI")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.grid(alpha=.2)
    fig.autofmt_xdate()
    _finish(fig, path)


def show_derivatives(grid, curve, path):
    _figure()
    valid = np.flatnonzero(np.isfinite(curve))
    fig, axes = plt.subplots(3, 1, figsize=(12, 7), layout="constrained", sharex=True)
    if len(valid) >= 3:
        y = curve[valid]
        d1 = np.gradient(y, valid)
        d2 = np.gradient(d1, valid)
        for ax, data, label in zip(axes, [y, d1, d2], ["NDVI", "一阶导数", "二阶导数"]):
            ax.plot(grid[valid], data)
            ax.set_ylabel(label)
            ax.grid(alpha=.2)
    fig.suptitle("速率/加速度特征不是同一个阈值事件")
    _finish(fig, path)


def run_regions(root, windows, sg_settings, functions, folder):
    qa, daily_fn, interpolate_fn, sg_fn, threshold_fn, stats_fn = functions
    maps, tables, checks = {}, [], []
    for key, window in windows.items():
        roi = "NC_ROT" if key.startswith("NC_") else key
        start, end = window
        ds = load_region(root, roi, full=True)
        valid = qa(ds.ndvi_raw.values, ds.fmask.values)
        values = np.where(valid & ds.inside_roi.values[None, :, :], ds.ndvi_raw.values, np.nan)
        grid, observed = daily_fn(ds.time.values, values, start, end)
        interpolated = interpolate_fn(observed)
        smooth = sg_fn(interpolated, sg_settings)
        cube = np.full((4, ds.sizes["y"], ds.sizes["x"]), np.nan)
        for row, col in np.argwhere(ds.inside_roi.values):
            cube[:, row, col] = threshold_fn(smooth[:, row, col], ratio=.2)
        # Same algorithm and support as the representative pixel; verify before DOY conversion.
        row, col = int(ds.attrs["sample_row"]), int(ds.attrs["sample_col"])
        standalone = sg_fn(interpolate_fn(observed[:, row, col]), sg_settings)
        np.testing.assert_allclose(smooth[:, row, col], standalone, equal_nan=True, atol=1e-12)
        np.testing.assert_allclose(cube[:, row, col], threshold_fn(standalone, .2), equal_nan=True, atol=1e-8)
        checks.append(dict(ROI=key, shape=[ds.sizes["y"], ds.sizes["x"]], representative_matches=True))
        origin_doy = (grid[0]-pd.Timestamp(grid[0].year, 1, 1)).days+1
        cube[:3] += origin_doy
        out = xr.Dataset({metric: (("y", "x"), cube[i].astype("float32")) for i, metric in enumerate(METRICS)},
                          coords={"y": ds.y, "x": ds.x},
                          attrs={"roi": roi, "season": key, "crs": ds.attrs["crs"], "transform": ds.attrs["transform"],
                                 "window_start": start, "window_end": end, "date_origin_year": int(grid[0].year),
                                 "method": "provided_Chen_SG_then_20_percent", "scope": ds.attrs["spatial_scope"],
                                 "SG_settings": json.dumps(sg_settings), "sample_source": ds.attrs["sample_source"]})
        out["inside_roi"] = ds.inside_roi.astype("uint8")
        for metric in METRICS:
            out[metric].attrs["units"] = "days" if metric == "LOS" else f"DOY relative to {grid[0].year}-01-01"
        out.to_netcdf(folder / f"{key}_full_chen20.nc", engine="h5netcdf")
        with rasterio.open(folder / f"{key}_full_chen20.tif", "w", driver="GTiff", height=ds.sizes["y"],
                            width=ds.sizes["x"], count=4, dtype="float32", crs=ds.attrs["crs"],
                            transform=Affine(*json.loads(ds.attrs["transform"])), nodata=np.nan, compress="deflate") as dst:
            dst.write(cube.astype("float32"))
            for i, name in enumerate(METRICS, 1):
                dst.set_band_description(i, name)
        maps[key] = out
        tables.append(stats_fn(cube, key, start, end))
        print(f"{key}: complete {ds.sizes['y']} x {ds.sizes['x']} native grid", flush=True)
    stats = pd.concat(tables, ignore_index=True)
    stats.to_csv(folder / "roi_statistics.csv", index=False, encoding="utf-8-sig")
    (folder / "spatial_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf8")
    return maps, stats


def show_distributions(maps, keys, path, title):
    _figure()
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), layout="constrained")
    for ax, metric in zip(axes.flat, METRICS):
        values = [maps[key][metric].values[np.isfinite(maps[key][metric].values)] for key in keys]
        ax.boxplot([a if len(a) else [np.nan] for a in values], tick_labels=[LABELS[k] for k in keys],
                   showfliers=False)
        ax.set(title=metric, ylabel="天" if metric == "LOS" else "DOY")
        ax.grid(axis="y", alpha=.2)
    fig.suptitle(title+"（箱体为Q1—Q3；散点未绘制，不从统计中剔除）")
    _finish(fig, path)


def show_maps(ds, path, title):
    _figure()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for ax, metric in zip(axes.flat, METRICS):
        args = {} if metric == "LOS" else dict(vmin=1, vmax=366)
        im = ax.imshow(ds[metric], interpolation="nearest", cmap="viridis", **args)
        ax.set(title=metric, xlabel="原生格网列", ylabel="原生格网行")
        fig.colorbar(im, ax=ax, label="天" if metric == "LOS" else "DOY", shrink=.8)
    fig.suptitle(title+"；完整ROI，空白为未得到该指标或区外")
    _finish(fig, path)


def comparison_table(stats, pairs):
    rows = []
    for a, b in pairs:
        for metric in METRICS:
            left = stats.loc[(stats.ROI == a) & (stats.metric == metric)].iloc[0]
            right = stats.loc[(stats.ROI == b) & (stats.metric == metric)].iloc[0]
            rows.append(dict(A=a, B=b, metric=metric, A_median=left["median"], B_median=right["median"],
                             difference_B_minus_A=right["median"]-left["median"],
                             A_IQR=left.IQR, B_IQR=right.IQR, A_N=left.N, B_N=right.N))
    return pd.DataFrame(rows)
