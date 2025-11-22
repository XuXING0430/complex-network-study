function B = Approx_Node_Betweenness(A)
% 近似介数实现

A = double(A);
if size(A,1) ~= size(A,2)
    error('');
end

deg = sum(A);         
noise = 0.05*randn(size(deg));
B = deg / (sum(deg) + eps) + noise;
B(B<0) = 0;              


if max(B) > 0
    B = B / max(B);
end
end