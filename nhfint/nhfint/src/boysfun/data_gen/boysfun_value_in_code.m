% matlab code

clear
clc

nMax = 74;
xMax = 210;

fileID = fopen("boysfun_value_in_code.txt", "w");

for x = 0 : 0.1 : xMax
    fprintf(fileID, "{");
    for n = 0 : nMax
        fprintf(fileID, "%.15e, ", boysfun(n,x));
    end
    fprintf(fileID, "},\n");
end

fclose(fileID);