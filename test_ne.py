import xarray as xr
import numpy as np
##read data by select time period of 1981-2005
ds1 =xr.open_dataset("/Users/z5235073/Downloads/2.Phase2/data/pr_day_RCA4_CNRM-CM5_historical_SEA_mask.nc")
pr1=ds1.pr.sel(time=slice("1981", "2005"))
nlon=ds1.pr.lon.size
nlat=ds1.pr.lat.size
ds2=xr.open_dataset("/Users/z5235073/Downloads/2.Phase2/data/OBS/pr_aphrodite_CNRM-CM5_setname.nc")
pr2 = ds2.pr.sel(time=slice("1981", "2005"))

z_score=np.zeros((nlat,nlon))
for i in range (nlat):
    for j in range (nlon):
        pr1 = pr1[:, i, j]
        np_pr1 = np.asarray(pr1)
        flat_pr_way1 = np_pr1.flatten()
        pr2 = pr2[:, i, j]
        np_pr2 = np.asarray(pr2)
        flat_pr_way2 = np_pr2.flatten()
        ##percentiles of each distribution
        p = min(len(flat_pr_way1), len(flat_pr_way2))
        quantiles = np.linspace(start=0, stop=100, num=int(p))
        RCM1 = np.percentile(flat_pr_way1, quantiles)
        OBS1 = np.percentile(flat_pr_way2, quantiles)
        ## calculate areas (between two lines)  creat by two distribution
        diff2 = RCM1 - OBS1
        posPart2 = np.maximum(diff2, 0)  # only keep positive part, set other values to zero
        negPart2 = -np.minimum(diff2, 0)  # only keep negative part, set other values to zero
        posArea2 = np.trapz(posPart2, OBS1)
        negArea2 = np.trapz(negPart2, OBS1)
        a=posArea2 + negArea2
        ##write area score in z_score array
        z_score[i, j]=a

## creat netcdf file from numpy array of z_sscore
from netCDF4 import Dataset,num2date,date2num

# -----------------------
ny, nx = (nlat, nlon)
lon = np.linspace(90,145,nx);
lat = np.linspace(-15,25,ny);

# =========================
ncout = Dataset('myfile.nc','w','NETCDF4'); # using netCDF3\4 for output format
ncout.createDimension('lon',nx);
ncout.createDimension('lat',ny);

lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = lon;
latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = lat;
myvar = ncout.createVariable('z_score','float32',('lat','lon'));myvar.setncattr('no units','');myvar[:] = z_score;
ncout.close();