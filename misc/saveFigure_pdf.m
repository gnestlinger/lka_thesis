function saveFigure_pdf(fig,filename,delOnOff)
% saveFigure_pdf(fig,filename,delOnOff)
% prints figure fig to filename.eps and converts to filename.pdf
% if delOnOff is set to 1, filename.eps is deleted

if nargin < 3
    delOnOff = 0; 
end%if

if (delOnOff ~= 0) && (delOnOff ~= 1)
    disp('ERROR: delOnOff can only have values 0 or 1');
    return;
end%if

print(fig, '-depsc2', filename);
system(['epstopdf ', filename,'.eps']);

if delOnOff == 1
    delete([filename,'.eps']);
end%if

end