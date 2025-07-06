function prox_x = prox_l1(x, lambda)
    
    % prox_{lambda*||.||_1}(z)
    prox_x = sign(x) .* max(abs(x) - lambda, 0);
end
