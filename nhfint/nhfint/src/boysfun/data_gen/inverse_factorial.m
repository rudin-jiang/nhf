clc;
clear;

for i = 0 : 9
    fprintf("%.15e, ", 1.0 / factorial(i));
    if (mod(i+1, 5) == 0)
        fprintf("\n");
    end
end