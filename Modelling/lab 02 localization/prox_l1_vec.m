function prox_x = prox_l1_vec(x, lambda,n,q)
    
    % prox_{lambda*||.||_1}(z)
    prox_x = sign(x) .* max(abs(x) - [lambda(1)*ones(n,1);lambda(2)*ones(q,1)], zeros(n+q,1));
end
