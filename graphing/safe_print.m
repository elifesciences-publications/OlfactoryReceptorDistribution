function safe_print(fname)
% safe_print Print image to PDF, avoiding overwrite.

if exist(fname, 'file')
    warning([mfilename ':fexist'], 'File already exists.');
    return;
end
print('-dpdf', fname);

end