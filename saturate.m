function im = saturate(im,range)

if nargin == 1
    range = [0 255];
end

a = range(1);
b = range(2);

im(im<a)=a;
im(im>b)=b;