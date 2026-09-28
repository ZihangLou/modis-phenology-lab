"""Annual greenness and recovery teaching support; no changes to previous labs."""
from pathlib import Path
from datetime import datetime
import calendar
import json
import hashlib
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
import pymannkendall as mk
from . import modis_teaching_support as prior
from .modis_dl import teaching_mask, extract_dl

QA = dict(summary_max=1, vi_quality_max=1, usefulness_max=4,
          adjacent_cloud=True, mixed_cloud=True, possible_snow=True, possible_shadow=True)
METRICS = ["annual_mean", "annual_max", "annual_min", "annual_amplitude", "annual_integral"]
CLIMATE_COLUMNS = ["year", "annual_precip", "mean_temp", "spring_temp", "autumn_temp"]
TERRACLIMATE_VARIABLES = ["tmin", "tmax", "tmean", "ppt", "soil"]


def new_run(root):
    folder = root / "outputs/modis_recovery" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    (folder / "figures").mkdir(parents=True)
    return folder


def open_region(root, region_id, start_year, end_year, index):
    if region_id != "beijing_download_bbox":
        raise ValueError("Only beijing_download_bbox is configured; register actual input paths for another region.")
    if index not in {"evi", "ndvi"}:
        raise ValueError("INDEX must be evi or ndvi")
    catalog = prior.inventory(root)
    catalog = catalog[catalog.year.between(start_year, end_year)].copy()
    if catalog.year.tolist() != list(range(start_year, end_year+1)):
        raise ValueError("Requested years are not all present in the MODIS source catalog")
    first = prior.load_year(catalog, start_year)
    first[index+"_analysis"] = first[index].where(teaching_mask(first, index, QA))
    return catalog, first


def annual_metrics(t, values):
    """Equal-composite mean; trapezoids join valid samples without endpoint extrapolation."""
    t, values = np.asarray(t, float), np.asarray(values, float)
    if values.shape[0] != len(t) or np.any(np.diff(t) <= 0):
        raise ValueError("Time must be strictly increasing and match values")
    valid = np.isfinite(values)
    count = valid.sum(axis=0)
    mean = np.divide(np.where(valid, values, 0).sum(axis=0), count,
                     out=np.full(count.shape, np.nan), where=count > 0)
    maximum = np.max(np.where(valid, values, -np.inf), axis=0)
    minimum = np.min(np.where(valid, values, np.inf), axis=0)
    maximum = np.where(count > 0, maximum, np.nan)
    minimum = np.where(count > 0, minimum, np.nan)
    integral = np.zeros(count.shape)
    duration = np.zeros(count.shape)
    previous_t = np.full(count.shape, np.nan)
    previous_y = np.full(count.shape, np.nan)
    for day, y in zip(t, values):
        present = np.isfinite(y)
        pair = present & np.isfinite(previous_y)
        integral += np.where(pair, .5*(y+previous_y)*(day-previous_t), 0)
        duration += np.where(pair, day-previous_t, 0)
        previous_t = np.where(present, day, previous_t)
        previous_y = np.where(present, y, previous_y)
    return {"annual_mean": mean, "annual_max": maximum, "annual_min": minimum,
            "annual_amplitude": maximum-minimum, "annual_integral": np.where(count >= 2, integral, np.nan),
            "integral_span_days": np.where(count >= 2, duration, np.nan)}


def finite_mean(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    return float(finite.mean()) if len(finite) else np.nan


def prepare_annual(catalog, first, index, region_id, run_dir, metric_function=annual_metrics):
    maps, records, regional_series, coverage = [], [], [], []
    for row in catalog.itertuples():
        ds = first if row.year == int(first.attrs["year"]) else prior.load_year(catalog, row.year)
        if not np.array_equal(ds.lat, first.lat) or not np.array_equal(ds.lon, first.lon):
            raise ValueError("Native grids differ between years")
        masked = ds[index].where(teaching_mask(ds, index, QA))
        values = masked.values
        t = pd.DatetimeIndex(ds.time.values).dayofyear.to_numpy(float)
        fields = metric_function(t, values)
        maps.append(xr.Dataset({key: (("lat", "lon"), np.asarray(array, np.float32)) for key, array in fields.items()},
                               coords={"lat": ds.lat, "lon": ds.lon}).expand_dims(year=[row.year]))
        regional = masked.mean(dim=("lat", "lon"), skipna=True).values
        regional_series.append(pd.DataFrame({"date": ds.time.values, "year": row.year, "doy": t, "vi": regional}))
        records.append({"year": row.year, **{key: finite_mean(array) for key, array in fields.items()}})
        coverage.append({"year": row.year, "composites": len(t), "height": ds.sizes["lat"],
                         "width": ds.sizes["lon"], "valid_observation_fraction": float(np.isfinite(values).mean())})
        print(f"{row.year}: annual metrics on {ds.sizes['lat']} x {ds.sizes['lon']} native grid", flush=True)
    annual_cube = xr.concat(maps, dim="year")
    annual_cube.attrs.update(index=index, region_id=region_id, QA=json.dumps(QA),
                             scope="full_downloaded_grid_equal_pixel_summary_NOT_administrative_boundary",
                             integral="trapezoidal between first and last valid composite; gaps bridged; no extrapolation")
    for key in METRICS:
        annual_cube[key].attrs["units"] = f"{index.upper()} days" if key == "annual_integral" else index.upper()
    annual_cube.integral_span_days.attrs["units"] = "days"
    annual_cube.to_netcdf(run_dir / "annual_greenness.nc", engine="h5netcdf")
    table, series, cov = pd.DataFrame(records), pd.concat(regional_series, ignore_index=True), pd.DataFrame(coverage)
    for name, frame in [("annual_greenness", table), ("regional_composites", series), ("coverage", cov)]:
        frame.to_csv(run_dir / f"{name}.csv", index=False, encoding="utf-8-sig")
    return annual_cube, table, series, cov


def trend_summary(years, values, alpha=.05):
    x, y = np.asarray(years, float), np.asarray(values, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    order = np.argsort(x)
    x, y = x[order], y[order]
    if len(x) < 3 or np.any(np.diff(x) <= 0):
        return dict(linear_slope=np.nan, linear_p=np.nan, sen_slope=np.nan, sen_low=np.nan,
                    sen_high=np.nan, mk_trend="not computable", mk_p=np.nan, tau=np.nan, n=len(x))
    linear = stats.linregress(x, y)
    sen = stats.theilslopes(y, x, alpha=1-alpha)
    result = mk.original_test(y, alpha=alpha)
    # scipy.linregress returns undefined r/p for constant input in some versions.
    p = 1. if np.ptp(y) == 0 else float(linear.pvalue)
    return dict(linear_slope=float(linear.slope), linear_p=p, sen_slope=float(sen.slope),
                sen_low=float(sen.low_slope), sen_high=float(sen.high_slope),
                mk_trend=result.trend, mk_p=float(result.p), tau=float(result.Tau), n=len(x))


def pixelwise_trend(years, field):
    """Vectorized OLS and two-sided Student-t p; actual years, pairwise finite samples."""
    values = np.asarray(field, float)
    years = np.asarray(years, float)
    if len(years) != len(values) or np.any(np.diff(years) <= 0):
        raise ValueError("Years must be unique, increasing and match the field")
    x = (years-years[0]).reshape((-1,)+(1,)*(values.ndim-1))
    valid = np.isfinite(values)
    n = valid.sum(axis=0)
    mean_x = np.divide(np.where(valid, x, 0).sum(axis=0), n, out=np.zeros(n.shape), where=n > 0)
    mean_y = np.divide(np.where(valid, values, 0).sum(axis=0), n, out=np.zeros(n.shape), where=n > 0)
    dx, dy = np.where(valid, x-mean_x, 0), np.where(valid, values-mean_y, 0)
    sxx, syy, sxy = (dx*dx).sum(axis=0), (dy*dy).sum(axis=0), (dx*dy).sum(axis=0)
    slope = np.divide(sxy, sxx, out=np.full(n.shape, np.nan), where=sxx > 0)
    intercept = mean_y-slope*mean_x
    residual = np.where(valid, values-(intercept+x*slope), 0)
    sse = (residual*residual).sum(axis=0)
    variance = np.divide(sse, n-2, out=np.full(n.shape, np.nan), where=n > 2)
    se = np.sqrt(np.divide(variance, sxx, out=np.full(n.shape, np.nan), where=sxx > 0))
    statistic = np.divide(np.abs(slope), se, out=np.full(n.shape, np.nan), where=se > 0)
    p = 2*stats.t.sf(statistic, np.maximum(n-2, 1))
    p = np.where((n >= 3) & (syy == 0), 1., p)
    p = np.where((n >= 3) & (se == 0) & (slope != 0), 0., p)
    return xr.Dataset({"slope": (field.dims[1:], slope), "p": (field.dims[1:], p)},
                      coords={d: field.coords[d] for d in field.dims[1:]})


def regional_phenology(series, region_id, index, run_dir):
    """Reuse the exact DL extractor, not a mismatched sample-point annual table."""
    source_hash = hashlib.sha256(pd.util.hash_pandas_object(series, index=False).values.tobytes()).hexdigest()
    rows = []
    for year, part in series.groupby("year", sort=True):
        daily = np.arange(1., pd.Timestamp(int(year), 12, 31).dayofyear+1)
        values = extract_dl(part.doy.to_numpy(float), part.vi.to_numpy(float), daily)
        rows.append({"year": int(year), **dict(zip(["SOS", "POS", "EOS", "LOS"], values))})
    table = pd.DataFrame(rows)
    table.to_csv(run_dir / "regional_dl_phenology.csv", index=False, encoding="utf-8-sig")
    metadata = dict(region_id=region_id, index=index, QA=QA, source_series_sha256=source_hash,
                    method="DL_parameter", support="regional_mean_curve_NOT_mean_pixel_dates",
                    definition="SOS=rise; EOS=rise+duration; POS=daily maximum; LOS=duration")
    (run_dir / "regional_dl_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf8")
    return table


def load_spatial_dl(root, years, index, reference, directory=None):
    """Only accept full-grid, matching-QA DL runs; never verification cubes or Zhang products."""
    candidates = [Path(directory)] if directory else sorted((root / "outputs/modis_dl_teaching").glob("*"), reverse=True)
    rejected = []
    for folder in candidates:
        metadata_path = folder / "configuration.json"
        if not metadata_path.is_file():
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf8"))
        if (metadata.get("scope") != "full_downloaded_grid" or metadata.get("method") != "DL_parameter"
                or metadata.get("index") != index or metadata.get("QA") != QA
                or metadata.get("years") != list(years)):
            rejected.append(str(folder))
            continue
        files = [folder / f"DL_{index}_{int(y)}.nc" for y in years]
        if not all(p.is_file() for p in files):
            rejected.append(str(folder))
            continue
        arrays = []
        try:
            for year, path in zip(years, files):
                with xr.open_dataset(path, engine="h5netcdf") as ds:
                    a = ds.phenology.load()
                if (a.attrs.get("scope") != "full_downloaded_grid" or a.attrs.get("method") != "DL_parameter"
                        or a.attrs.get("index") != index or int(a.attrs.get("year", -1)) != year
                        or not np.array_equal(a.lat, reference.lat) or not np.array_equal(a.lon, reference.lon)
                        or a.metric.values.tolist() != ["SOS", "POS", "EOS", "LOS"]):
                    raise ValueError("DL grid, year, method or metric mismatch")
                arrays.append(a.expand_dims(year=[year]))
        except ValueError:
            rejected.append(str(folder))
            continue
        return xr.concat(arrays, dim="year"), {"state": "loaded", "source": str(folder)}
    return None, {"state": "pending_full_grid_DL", "rejected_runs": rejected,
                  "message": "尚无匹配年份、指数、QA与完整格网的DL产品；3×3验证和Zhang结果不用于本图。"}


def load_climate_grid_dl(root, years, index, climate_reference, directory=None):
    """Load DL phenology fitted after MODIS was area-averaged to the climate grid."""
    candidates = ([Path(directory)] if directory else
                  sorted((root / "outputs/modis_dl_climate_grid").glob("*"), reverse=True))
    expected_years = [int(year) for year in years]
    rejected = []
    for folder in candidates:
        config = folder / "configuration.json"
        if not config.is_file():
            continue
        metadata = json.loads(config.read_text(encoding="utf8"))
        if (metadata.get("scope") != "climate_grid_from_modis_area_average"
                or metadata.get("method") != "DL_parameter" or metadata.get("index") != index
                or metadata.get("QA") != QA or metadata.get("years") != expected_years):
            rejected.append(str(folder))
            continue
        files = [folder / f"DL_{index}_{year}.nc" for year in expected_years]
        if not all(path.is_file() for path in files):
            rejected.append(str(folder))
            continue
        arrays = []
        try:
            for year, path in zip(expected_years, files):
                with xr.open_dataset(path, engine="h5netcdf") as ds:
                    array = ds.phenology.load()
                if (array.metric.values.tolist() != ["SOS", "POS", "EOS", "LOS"]
                        or not np.array_equal(array.lat, climate_reference.lat)
                        or not np.array_equal(array.lon, climate_reference.lon)
                        or int(array.attrs.get("year", -1)) != year):
                    raise ValueError("Climate-grid DL metadata or coordinates mismatch")
                arrays.append(array.expand_dims(year=[year]))
        except ValueError:
            rejected.append(str(folder))
            continue
        return xr.concat(arrays, dim="year"), {"state": "loaded", "source": str(folder)}
    return None, {"state": "pending_climate_grid_DL", "rejected_runs": rejected,
                  "message": "尚无完整的MODIS聚合到TerraClimate格网后的多年DL产品。"}


def read_climate(path, years, region_id):
    if not Path(path).is_file():
        return pd.DataFrame(columns=CLIMATE_COLUMNS), "待补充真实气候CSV；本次不计算气候相关。"
    table = pd.read_csv(path)
    missing = set(CLIMATE_COLUMNS)-set(table.columns)
    if missing:
        raise ValueError(f"Climate columns missing: {sorted(missing)}")
    if "region_id" in table and not table.region_id.eq(region_id).all():
        raise ValueError("Climate region_id differs from this experiment")
    for column in CLIMATE_COLUMNS:
        table[column] = pd.to_numeric(table[column], errors="raise")
    if table.year.isna().any() or (table.year % 1 != 0).any() or table.year.duplicated().any():
        raise ValueError("Climate years must be unique nonmissing integers")
    table.year = table.year.astype(int)
    if np.isinf(table[CLIMATE_COLUMNS].to_numpy(float)).any():
        raise ValueError("Climate CSV contains infinite values")
    return table[table.year.isin(years)].sort_values("year"), "真实CSV已读取；请依据数据说明核对区域、单位和季节月份。"


def build_terraclimate_annual_csv(source_dir, output_path, region_id,
                                  start_year=2001, end_year=2025):
    """Validate annual TerraClimate stacks and derive the climate CSV used by lab 04."""
    import rasterio

    source_dir, output_path = Path(source_dir), Path(output_path)
    files = sorted(source_dir.glob("Beijing_TerraClimate_*.tif"))
    if not files:
        raise FileNotFoundError(f"No TerraClimate GeoTIFFs under {source_dir}")
    expected_descriptions = tuple(
        f"m{month:02d}_{variable}"
        for month in range(1, 13) for variable in TERRACLIMATE_VARIABLES
    )
    rows, reference_grid, used_files = [], None, []
    for path in files:
        try:
            year = int(path.name.split("_")[2])
        except ValueError as exc:
            raise ValueError(f"TerraClimate filename has no valid year: {path}") from exc
        if not start_year <= year <= end_year:
            continue
        with rasterio.open(path) as src:
            grid = (str(src.crs), src.width, src.height, tuple(src.transform), tuple(src.bounds))
            if reference_grid is None:
                reference_grid = grid
            elif grid != reference_grid:
                raise ValueError(f"TerraClimate grid mismatch: {path}")
            if src.count != 60 or tuple(src.descriptions) != expected_descriptions:
                raise ValueError(f"TerraClimate bands are not the expected 12 x 5 monthly stack: {path}")
            if src.dtypes != ("float32",) * 60 or src.nodata != -9999:
                raise ValueError(f"TerraClimate dtype or NoData mismatch: {path}")
            data = src.read(masked=True)

        monthly_temp = np.array([finite_mean(data[(month-1)*5+2]) for month in range(1, 13)])
        monthly_precip = np.array([finite_mean(data[(month-1)*5+3]) for month in range(1, 13)])
        if not np.isfinite(monthly_temp).all() or not np.isfinite(monthly_precip).all():
            raise ValueError(f"TerraClimate has a month with no valid temperature or precipitation: {path}")
        days = np.array([calendar.monthrange(year, month)[1] for month in range(1, 13)], float)
        spring = np.array([2, 3, 4])
        autumn = np.array([8, 9, 10])
        rows.append({
            "year": year,
            "annual_precip": float(monthly_precip.sum()),
            "mean_temp": float(np.average(monthly_temp, weights=days)),
            "spring_temp": float(np.average(monthly_temp[spring], weights=days[spring])),
            "autumn_temp": float(np.average(monthly_temp[autumn], weights=days[autumn])),
            "region_id": region_id,
        })
        used_files.append(str(path))

    table = pd.DataFrame(rows).sort_values("year")
    if table.empty or table.year.duplicated().any() or not table.year.is_monotonic_increasing:
        raise ValueError("TerraClimate years must be nonempty, unique and increasing")
    expected_available = sorted(
        int(folder.name) for folder in source_dir.iterdir()
        if folder.is_dir() and folder.name.isdigit() and start_year <= int(folder.name) <= end_year
    )
    if table.year.tolist() != expected_available:
        raise ValueError("TerraClimate year folders and readable annual stacks differ")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False, encoding="utf-8-sig")
    metadata = {
        "dataset_id": "IDAHO_EPSCOR/TERRACLIMATE",
        "region_id": region_id,
        "roi": "Beijing boundary buffered 10 km, then rectangular bounds",
        "years_available": table.year.tolist(),
        "years_missing_from_requested_range": sorted(set(range(start_year, end_year+1))-set(table.year)),
        "variables": {"annual_precip": "sum of 12 monthly ppt; mm/year",
                      "mean_temp": "day-weighted mean of monthly tmean; degree Celsius",
                      "spring_temp": "day-weighted March-May tmean; degree Celsius",
                      "autumn_temp": "day-weighted September-November tmean; degree Celsius"},
        "spatial_aggregation": "equal-weight mean of valid 1/24 degree cells in the downloaded rectangular ROI",
        "crs": reference_grid[0],
        "width": reference_grid[1],
        "height": reference_grid[2],
        "transform": list(reference_grid[3]),
        "source_files": used_files,
    }
    metadata_path = output_path.with_suffix(".metadata.json")
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf8")
    return table, metadata


def load_terraclimate_spatial(source_dir, start_year=2001, end_year=2025):
    """Load monthly stacks and derive annual/seasonal fields on the native climate grid."""
    import rasterio

    source_dir = Path(source_dir)
    files = sorted(source_dir.glob("Beijing_TerraClimate_*.tif"))
    expected = tuple(f"m{m:02d}_{v}" for m in range(1, 13) for v in TERRACLIMATE_VARIABLES)
    years, records, reference = [], [], None
    for path in files:
        year = int(path.name.split("_")[2])
        if not start_year <= year <= end_year:
            continue
        with rasterio.open(path) as src:
            signature = (str(src.crs), src.width, src.height, tuple(src.transform), tuple(src.bounds))
            if reference is None:
                reference = signature
                transform = src.transform
            elif signature != reference:
                raise ValueError(f"TerraClimate grid mismatch: {path}")
            if src.count != 60 or tuple(src.descriptions) != expected or src.nodata != -9999:
                raise ValueError(f"Unexpected TerraClimate stack structure: {path}")
            data = src.read(masked=True).filled(np.nan).astype(float)
        tmean = data[2::5]
        ppt = data[3::5]
        days = np.array([calendar.monthrange(year, month)[1] for month in range(1, 13)], float)

        def weighted_mean(indices):
            values = tmean[indices]
            weights = days[indices, None, None]
            valid = np.isfinite(values)
            return np.divide(np.where(valid, values*weights, 0).sum(axis=0),
                             np.where(valid, weights, 0).sum(axis=0),
                             out=np.full(values.shape[1:], np.nan),
                             where=np.where(valid, weights, 0).sum(axis=0) > 0)

        records.append({"mean_temp": weighted_mean(np.arange(12)),
                        "spring_temp": weighted_mean(np.array([2, 3, 4])),
                        "autumn_temp": weighted_mean(np.array([8, 9, 10])),
                        "annual_precip": np.nansum(ppt, axis=0)})
        years.append(year)
    if not records:
        raise FileNotFoundError(f"No TerraClimate annual stacks under {source_dir}")
    if years != sorted(set(years)):
        raise ValueError("TerraClimate years must be unique and increasing")
    lon = transform.c + (np.arange(reference[1])+.5)*transform.a
    lat = transform.f + (np.arange(reference[2])+.5)*transform.e
    ds = xr.Dataset({name: (("year", "lat", "lon"),
                            np.stack([record[name] for record in records]).astype(np.float32))
                     for name in ["mean_temp", "spring_temp", "autumn_temp", "annual_precip"]},
                    coords={"year": years, "lat": lat, "lon": lon})
    ds.attrs.update(dataset_id="IDAHO_EPSCOR/TERRACLIMATE", crs=reference[0],
                    transform=json.dumps(list(transform)), spatial_resolution="1/24 degree",
                    support="native climate grid; not upsampled to MODIS")
    return ds


def aggregate_to_climate_grid(field, climate_reference, return_coverage=False):
    """Overlap-area weighted mean on a regular geographic grid (spherical area).

    Latitude overlaps use sin(lat); longitude overlaps use radians. The common
    Earth-radius factor cancels in the weighted mean and coverage fraction.
    Missing source pixels and areas outside the downloaded ROI reduce coverage.
    """
    if field.dims[-2:] != ("lat", "lon"):
        raise ValueError("Field must end with lat, lon dimensions")

    def intervals(coordinates):
        centers = np.asarray(coordinates, float)
        steps = np.diff(centers)
        if len(centers) < 2 or not np.allclose(steps, steps[0]) or steps[0] == 0:
            raise ValueError("Regular, monotonic coordinates are required")
        half = abs(steps[0]) / 2
        return centers-half, centers+half

    sy0, sy1 = intervals(field.lat)
    sx0, sx1 = intervals(field.lon)
    dy0, dy1 = intervals(climate_reference.lat)
    dx0, dx1 = intervals(climate_reference.lon)
    ylo = np.maximum(dy0[:, None], sy0[None, :])
    yhi = np.minimum(dy1[:, None], sy1[None, :])
    wy = np.maximum(np.sin(np.deg2rad(yhi))-np.sin(np.deg2rad(ylo)), 0)
    wx = np.deg2rad(np.maximum(np.minimum(dx1[:, None], sx1[None, :])
                              - np.maximum(dx0[:, None], sx0[None, :]), 0))
    cell_area = np.outer(np.sin(np.deg2rad(dy1))-np.sin(np.deg2rad(dy0)),
                         np.deg2rad(dx1-dx0))
    values = np.asarray(field, float)
    destination = np.full(values.shape[:-2]+cell_area.shape, np.nan)
    coverage = np.zeros_like(destination)
    for index in np.ndindex(values.shape[:-2]):
        valid = np.isfinite(values[index])
        support = wy @ valid.astype(float) @ wx.T
        total = wy @ np.where(valid, values[index], 0) @ wx.T
        destination[index] = np.divide(total, support, out=np.full_like(total, np.nan), where=support > 0)
        coverage[index] = np.clip(support/cell_area, 0, 1)
    coords = {d: field.coords[d] for d in field.dims[:-2]}
    coords.update(lat=climate_reference.lat, lon=climate_reference.lon)
    result = xr.DataArray(destination, dims=field.dims, coords=coords, name=field.name, attrs=field.attrs)
    result.attrs["spatial_handoff"] = "spherical overlap-area mean to native TerraClimate grid"
    fractions = result.copy(data=coverage).rename("coverage_fraction")
    fractions.attrs = {"units": "1", "definition": "valid MODIS area / full climate cell area"}
    return (result, fractions) if return_coverage else result


def interpolate_climate_to_modis_grid(climate, modis_reference, method="bilinear"):
    """Interpolate continuous climate fields onto the MODIS grid for aligned analysis."""
    from rasterio.transform import from_origin
    from rasterio.warp import reproject, Resampling

    methods = {"bilinear": Resampling.bilinear, "cubic": Resampling.cubic}
    if method not in methods:
        raise ValueError("Climate interpolation method must be bilinear or cubic")
    src_lat, src_lon = np.asarray(climate.lat), np.asarray(climate.lon)
    dst_lat, dst_lon = np.asarray(modis_reference.lat), np.asarray(modis_reference.lon)
    src_dx, src_dy = abs(np.diff(src_lon).mean()), abs(np.diff(src_lat).mean())
    dst_dx, dst_dy = abs(np.diff(dst_lon).mean()), abs(np.diff(dst_lat).mean())
    src_transform = from_origin(src_lon.min()-src_dx/2, src_lat.max()+src_dy/2, src_dx, src_dy)
    dst_transform = from_origin(dst_lon.min()-dst_dx/2, dst_lat.max()+dst_dy/2, dst_dx, dst_dy)
    output = {}
    for name, variable in climate.data_vars.items():
        values = np.asarray(variable, np.float32)
        if src_lat[0] < src_lat[-1]:
            values = values[..., ::-1, :]
        destination = np.full(values.shape[:-2]+(len(dst_lat), len(dst_lon)), np.nan, np.float32)
        for index in np.ndindex(values.shape[:-2]):
            reproject(values[index], destination[index], src_transform=src_transform, src_crs="EPSG:4326",
                      dst_transform=dst_transform, dst_crs="EPSG:4326", src_nodata=np.nan,
                      dst_nodata=np.nan, resampling=methods[method])
        output[name] = ((variable.dims[:-2]+("lat", "lon")), destination)
    result = xr.Dataset(output, coords={"year": climate.year, "lat": dst_lat, "lon": dst_lon})
    result.attrs.update(climate.attrs)
    result.attrs.update(interpolation=method, target_grid="MODIS 349x466",
                        interpretation="grid alignment only; no new sub-climate-pixel information")
    return result


def spatial_correlation(x, y, minimum_n=3):
    """Pairwise Pearson correlation through time for every shared grid cell."""
    x, y = xr.align(x, y, join="inner")
    if x.dims != y.dims or x.dims[0] != "year":
        raise ValueError("Inputs must share year, lat, lon dimensions")
    xv, yv = np.asarray(x, float), np.asarray(y, float)
    valid = np.isfinite(xv) & np.isfinite(yv)
    n = valid.sum(axis=0)
    mx = np.divide(np.where(valid, xv, 0).sum(axis=0), n, out=np.zeros(n.shape), where=n > 0)
    my = np.divide(np.where(valid, yv, 0).sum(axis=0), n, out=np.zeros(n.shape), where=n > 0)
    dx, dy = np.where(valid, xv-mx, 0), np.where(valid, yv-my, 0)
    ssx, ssy = (dx*dx).sum(axis=0), (dy*dy).sum(axis=0)
    r = np.divide((dx*dy).sum(axis=0), np.sqrt(ssx*ssy),
                  out=np.full(n.shape, np.nan), where=(ssx > 0) & (ssy > 0) & (n >= minimum_n))
    r = np.clip(r, -1, 1)
    statistic = np.abs(r)*np.sqrt(np.divide(n-2, np.maximum(1-r*r, np.finfo(float).eps)))
    p = np.where(n >= minimum_n, 2*stats.t.sf(statistic, n-2), np.nan)
    return xr.Dataset({"r": (("lat", "lon"), r.astype(np.float32)),
                       "p": (("lat", "lon"), p.astype(np.float32)),
                       "n": (("lat", "lon"), n.astype(np.int16))},
                      coords={"lat": x.lat, "lon": x.lon},
                      attrs={"years": f"{int(x.year.min())}-{int(x.year.max())}",
                             "test": "pairwise Pearson correlation; two-sided unadjusted p"})


def spatial_correlation_summary(results, alpha=.05):
    rows = []
    for relationship, result in results.items():
        if result is None:
            rows.append({"relationship": relationship, "state": "pending_spatial_phenology"})
            continue
        valid = np.isfinite(result.r) & np.isfinite(result.p)
        positive = valid & (result.p < alpha) & (result.r > 0)
        negative = valid & (result.p < alpha) & (result.r < 0)
        total = int(valid.sum())
        rows.append({"relationship": relationship, "state": "computed", "valid_cells": total,
                     "median_r": float(result.r.where(valid).median()),
                     "significant_positive_percent": 100*int(positive.sum())/total if total else np.nan,
                     "significant_negative_percent": 100*int(negative.sum())/total if total else np.nan})
    return pd.DataFrame(rows)


def correlation_row(frame, xname, yname):
    pair = frame[[xname, yname]].dropna()
    row = dict(relationship=f"{xname} vs {yname}", n=len(pair), pearson_r=np.nan,
               pearson_p=np.nan, spearman_rho=np.nan, spearman_p=np.nan)
    if len(pair) >= 3 and all(pair[c].nunique() > 1 for c in pair):
        p = stats.pearsonr(pair[xname], pair[yname])
        s = stats.spearmanr(pair[xname], pair[yname])
        row.update(pearson_r=p.statistic, pearson_p=p.pvalue, spearman_rho=s.statistic, spearman_p=s.pvalue)
    return row


def pettitt(years, values):
    years, values = np.asarray(years), np.asarray(values, float)
    present = np.isfinite(values)
    years, values = years[present], values[present]
    n = len(values)
    if n < 3:
        return dict(change_after_year=None, approximate_p=np.nan)
    ranks = stats.rankdata(values)
    u = 2*np.cumsum(ranks)-np.arange(1, n+1)*(n+1)
    split = int(np.argmax(np.abs(u[:-1])))
    k = abs(u[split])
    return dict(change_after_year=int(years[split]), next_observed_year=int(years[split+1]),
                statistic=float(k), approximate_p=float(min(1., 2*np.exp(-6*k*k/(n**3+n*n)))))


def monthly_stl(series):
    from statsmodels.tsa.seasonal import STL
    series = series.set_index("date").vi.sort_index()
    # Explicit temporal regularization, not period=23 on irregular composite labels.
    daily = series.reindex(pd.date_range(series.index.min(), series.index.max(), freq="D"))
    daily = daily.interpolate(method="time", limit_area="inside")
    monthly = daily.resample("MS").mean()
    if monthly.isna().any() or len(monthly) < 24:
        raise ValueError("STL requires at least two complete monthly cycles")
    fit = STL(monthly, period=12, robust=True).fit()
    return pd.DataFrame({"observed": monthly, "trend": fit.trend, "seasonal": fit.seasonal, "residual": fit.resid})
