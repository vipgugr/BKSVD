function rgb = od2rgb(od)
    rgb=256*exp(-od)-1;
end