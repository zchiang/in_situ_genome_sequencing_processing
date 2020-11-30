function [dr,ds] = getDistancesForEachSegment(data)

data = sortrows(data,'hg38_pos');

dx = abs(data.x_um(1:end-1)-data.x_um(2:end));
dy = abs(data.y_um(1:end-1)-data.y_um(2:end));
dz = abs(data.z_um(1:end-1)-data.z_um(2:end));

ds = abs(data.hg38_pos(1:end-1)-data.hg38_pos(2:end));

dr = sqrt(dx.^2+dy.^2+dz.^2);

end