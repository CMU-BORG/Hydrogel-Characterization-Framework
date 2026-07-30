function err = ErrEllipse(x,y,rx,ry,x0,y0)

dx = x - x0;
dy = y - y0;

err = sum( ( 1 - ( (dx.^2 / rx^2)  +  (dy.^2 / ry^2) ) ).^2 ) + 1000*(abs(rx-ry) > 0);
end