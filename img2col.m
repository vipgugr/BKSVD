%
%Convierte una imagen a una matriz fila. Cada capa RGB se incluye en una fila
%diferente.
%

function y = img2col(img)
    layer= img(:,:,1);
    y=transpose(layer(:));
    for i=2:size(img,3)
        layer=img(:,:,i);
        y=cat(1,y,transpose(layer(:)));
    end
end