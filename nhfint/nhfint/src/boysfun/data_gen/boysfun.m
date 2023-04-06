function ret = boysfun(n, x)
    f = @(t) (t .^(2*n) .* exp(-x .* t.^2));
    ret = integral(f, 0, 1, 'RelTol', 0.0, 'AbsTol', 0.0);
end