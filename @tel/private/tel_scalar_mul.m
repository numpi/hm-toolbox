function [T] = tel_scalar_mul(a,A)
T = A;
for i = 1:numel(T.D)
    T.D{i} = T.D{i} * a;
    T.V{i} = T.V{i} * a;
end
end