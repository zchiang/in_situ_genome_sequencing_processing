function dist = doPariwiseDistP(x,y,z)

xm = x*ones(1,numel(x));
ym = y*ones(1,numel(y));
zm = z*ones(1,numel(z));

dx = (xm-xm').^2;
dy = (ym-ym').^2;
dz = (zm-zm').^2;

dist = sqrt(dx+dy+dz);

end