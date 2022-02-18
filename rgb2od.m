function od = rgb2od(rgb)
%asummes double rgb
od=real(-log((rgb+1)/256));

end