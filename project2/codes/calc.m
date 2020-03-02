function b = calc(sigma, m, lambda, x, y)
% calculate dual summation part of objective function
b = 0;
for i = 1:m
    for j = 1:m
        b = b + lambda(i)*lambda(j)*y(i)*y(j)*exp(-(x(i,:) - x(j,:))*(x(i,:) - x(j,:)).' / (2*sigma^2));
    end
end
end