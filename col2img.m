%
%Convierte una imagen expresada como fila a una matriz tipica. Cada fila se
%interpreta como una capa de color
%

function y = col2img(imgcol)
    %Supongo img cuadrada!!!!! NECESARIO
    imsize=sqrt(size(imgcol,2));
    y=transpose(vec2mat(imgcol(1,:),imsize));
    for i=2:size(imgcol,1)
        aux=transpose(vec2mat(imgcol(i,:),imsize));
        y=cat(3,y,aux);
    end
end