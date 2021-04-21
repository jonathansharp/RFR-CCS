% importETOPO2

path = '/Volumes/2TB Hard Drive/ETOPO2/'; % this can be changed as necessary

ETOPO2 = netcdfreader(strcat(path,'ETOPO2v2c_f4.nc'));

ETOPO2.lat = double(ETOPO2.y);
ETOPO2.lon = double(ETOPO2.x);
ETOPO2.bottomdepth = -ETOPO2.z;

ETOPO2 = rmfield(ETOPO2,{'x' 'y' 'z'});

clear path