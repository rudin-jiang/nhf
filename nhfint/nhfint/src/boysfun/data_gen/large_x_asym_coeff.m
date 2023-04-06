clc;
clear;
nMax = 64;

for n = 0 : nMax
    val = sqrt(pi) * double_factorial(2*n-1) / (2.^(n+1));
    fprintf("%.15e, ", val);
    if (mod(n+1, 5) == 0)
        fprintf("\n");
    end
end

function result = double_factorial(n) 
    if (n == 1 || n==0) 
        result = 1; 
    else 
        if (mod(n,2) == 0) % the number is even 
            result = prod(2:2:n); 
        else % the number is odd 
            result = prod(1:2:n); 
        end
    end
end