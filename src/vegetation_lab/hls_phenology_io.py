"""Read HLS NDVI/FMASK scene stacks without modifying or resampling input files."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from rasterio.warp import transform as transform_xy, transform_geom
from rasterio.features import geometry_mask
import xarray as xr

REASONS = {1:'无效NDVI', 2:'无效FMASK', 4:'云', 8:'邻云/阴影', 16:'云影', 32:'雪冰', 64:'高气溶胶', 128:'水体'}


def quality_mask(ndvi, fmask, *, exclude_adjacent=True, exclude_high_aerosol=False, exclude_water=False):
    """HLS v2 bit0 unused, bits1/2/3/4 cloud/adjacent/shadow/snow, bit5 water, bits6-7 aerosol.

    Reject NaN/inf/-9999, Fmask 255 fill, non-integer QA BEFORE uint conversion.
    Range filtering does not clip invalid NDVI into plausible vegetation values.
    """
    ndvi, fmask = np.asarray(ndvi), np.asarray(fmask)
    nvalid = np.isfinite(ndvi) & (ndvi >= -1) & (ndvi <= 1)
    qvalid = np.isfinite(fmask) & (fmask >= 0) & (fmask < 255) & (fmask == np.rint(fmask))
    q = np.where(qvalid, fmask, 0).astype('uint16')
    reason = (~nvalid).astype('uint16') | ((~qvalid).astype('uint16') << 1)
    for bit, code, enabled in [(1,4,True),(2,8,exclude_adjacent),(3,16,True),(4,32,True),(5,128,exclude_water)]:
        if enabled: reason |= np.where(qvalid & ((q & (1<<bit)) != 0), code, 0).astype('uint16')
    if exclude_high_aerosol: reason |= np.where(qvalid & ((q>>6)==3),64,0).astype('uint16')
    return reason == 0, reason


def inventory(root):
    records=[]
    for folder in sorted(Path(root).glob('Ch08_HLS_NDVI_v1_*')):
        manifests=sorted(folder.glob('*band_manifest.csv'))
        if not manifests: continue
        table=pd.concat([pd.read_csv(p) for p in manifests],ignore_index=True)
        table=table.drop_duplicates(['output_file','band_index_1based'])
        dates=pd.to_datetime(table.acquisition_time_utc,utc=True)
        records.append(dict(roi=folder.name.removeprefix('Ch08_HLS_NDVI_v1_'),
                            first=str(dates.min().date()),last=str(dates.max().date()),
                            files=table.output_file.nunique(),scenes=int((table.variable=='NDVI').sum()),
                            grids=table[['grid_crs','grid_transform']].drop_duplicates().shape[0]))
    return pd.DataFrame(records)


def load_patch(root, roi, *, patch_size=15, grid_crs=None, grid_transform=None, center=None):
    """Read a center patch on ONE physical grid. No cross-grid mosaicking or reprojection.

    Prefer the center's UTM zone, then the grid with most scenes. Excluded-grid
    scene counts are explicit. patch_size=None reads the full exported grid for explicit spatial screening.
    This is a sample, not a class-verified plot.
    """
    if patch_size is not None and (patch_size < 1 or patch_size % 2 != 1): raise ValueError('Use odd patch_size >= 1')
    folder=Path(root)/('Ch08_HLS_NDVI_v1_'+roi)
    manifests=sorted(folder.glob('*band_manifest.csv'))
    if not manifests:raise FileNotFoundError(f'Manifest missing: {folder}')
    table=pd.concat([pd.read_csv(p) for p in manifests],ignore_index=True)
    # Repeated copies may be identical, but conflicting provenance must not be silently dropped.
    for _, group in table.groupby(['output_file','band_index_1based']):
        if len(group[['band_name','source_id','acquisition_time_utc']].drop_duplicates())>1:
            raise ValueError('Conflicting manifest rows')
    table=table.drop_duplicates(['output_file','band_index_1based'])
    geo=json.loads(sorted(folder.glob('*roi.geojson'))[-1].read_text(encoding='utf8'))
    if center is None:
        coords=np.asarray(geo['features'][0]['geometry']['coordinates'][0],float)
        center=tuple((coords.min(0)+coords.max(0))/2)
    lon,lat=center
    desired=f'EPSG:{32600+int((lon+180)//6)+1}'
    ndvi=table[table.variable.eq('NDVI')]
    counts=ndvi.groupby(['grid_crs','grid_transform']).size().reset_index(name='scenes')
    if grid_crs is not None: counts=counts[counts.grid_crs.eq(grid_crs)]
    if grid_transform is not None: counts=counts[counts.grid_transform.eq(grid_transform)]
    if counts.empty:raise ValueError('Requested native grid not found')
    counts['preferred']=counts.grid_crs.eq(desired)
    chosen=counts.sort_values(['preferred','scenes','grid_crs'],ascending=[False,False,True]).iloc[0]
    ndvi=ndvi[ndvi.grid_crs.eq(chosen.grid_crs)&ndvi.grid_transform.eq(chosen.grid_transform)].copy()
    ndvi=ndvi.drop_duplicates('source_id').sort_values(['time_start_ms','source_id'])
    qa=table[table.variable.eq('FMASK')].set_index(['output_file','source_id'])
    with rasterio.open(folder/ndvi.iloc[0].output_file) as src:
        x,y=transform_xy('EPSG:4326',src.crs,[lon],[lat]); row,col=src.index(x[0],y[0])
        if not 0 <= row < src.height or not 0 <= col < src.width:raise ValueError('Sample center not inside file')
        if patch_size is None:
            r0=c0=0
            window=Window(0,0,src.width,src.height)
        else:
            r0=max(0,row-patch_size//2);c0=max(0,col-patch_size//2)
            window=Window(c0,r0,min(patch_size,src.width-c0),min(patch_size,src.height-r0))
        shape=(int(window.height),int(window.width));tr=src.window_transform(window)
        reference=(src.crs,src.transform,src.width,src.height)
    n=len(ndvi);values=np.empty((n,*shape),np.float32);flags=values.copy()
    for file, rows in ndvi.groupby('output_file'):
        with rasterio.open(folder/file) as src:
            if (src.crs,src.transform,src.width,src.height)!=reference:raise ValueError('Export grids differ: '+file)
            for index, r in rows.iterrows():
                i=ndvi.index.get_loc(index); q=qa.loc[(file,r.source_id)]
                if isinstance(q,pd.DataFrame):raise ValueError('Duplicate QA source')
                for dest,band,description in [(values,int(r.band_index_1based),r.band_name),(flags,int(q.band_index_1based),q.band_name)]:
                    if src.descriptions[band-1]!=description:raise ValueError('Manifest/band description mismatch')
                    a=src.read(band,window=window,masked=True).filled(np.nan).astype('float32')
                    a[~np.isfinite(a)]=np.nan
                    dest[i]=a
    times=pd.to_datetime(ndvi.acquisition_time_utc,utc=True).dt.tz_localize(None).to_numpy()
    ds=xr.Dataset({'ndvi_raw':(('scene','y','x'),values),'fmask':(('scene','y','x'),flags)},
       coords={'scene':np.arange(n),'time':('scene',times),'source_id':('scene',ndvi.source_id.astype(str).to_numpy()),
               'branch':('scene',ndvi.branch.astype(str).to_numpy()),
               'y':tr.f+(np.arange(shape[0])+.5)*tr.e,'x':tr.c+(np.arange(shape[1])+.5)*tr.a},
       attrs={'roi_id':roi,'crs':str(reference[0]),'transform':json.dumps(list(tr)[:6]),
              'center_lon':float(lon),'center_lat':float(lat),'sample_row':row-r0,'sample_col':col-c0,
              'source_grids':int(len(table[['grid_crs','grid_transform']].drop_duplicates())),
              'excluded_other_grid_scenes':int((table.variable=='NDVI').sum()-n),
              'spatial_support':'single native grid center patch; no independently verified crop labels'})
    projected_roi = transform_geom('EPSG:4326', ds.attrs['crs'], geo['features'][0]['geometry'])
    inside = geometry_mask([projected_roi], out_shape=shape, transform=tr, invert=True)
    ds['inside_roi'] = (('y','x'), inside)
    ds.attrs['spatial_scope'] = 'full_downloaded_5km_roi_on_one_native_grid' if patch_size is None else 'center_patch'
    ds.attrs['roi_mask_rule'] = 'native pixel center inside downloaded ROI; no resampling'
    return ds


def apply_qc(ds, **kwargs):
    ds=ds.copy()
    valid,reason=quality_mask(ds.ndvi_raw.values,ds.fmask.values,**kwargs)
    ds['valid_qc']=(('scene','y','x'),valid)
    ds['qc_reason']=(('scene','y','x'),reason)
    ds['ndvi_qc']=ds.ndvi_raw.where(ds.valid_qc)
    ds.attrs['qa_policy']=json.dumps(kwargs,sort_keys=True)
    return ds


def pixel_series(ds, row=None, col=None):
    if row is None:row=int(ds.attrs['sample_row'])
    if col is None:col=int(ds.attrs['sample_col'])
    s=ds.isel(y=row,x=col)
    return pd.DataFrame({'date':s.time.values,'raw':s.ndvi_raw.values,'qa':s.fmask.values,
                         'qc':s.ndvi_qc.values,'branch':s.branch.values,'reason':s.qc_reason.values})
