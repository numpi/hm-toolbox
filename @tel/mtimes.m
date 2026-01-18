function T = mtimes(T1, T2)
%MTIMES Matrix multiplication for tel objects

if isa(T1, 'tel') && isa(T2, 'tel')
    % Both arguments are tel objects
    T = tel_mat_mat_mul(T1, T2); %FARE come???
elseif isa(T2, 'tel')
    % T2 is a tel object
    if isscalar(T1)
        T = tel_scalar_mul(T1, T2);
    else
        T = tel_vec_mat_mul(T1, T2); 
    end
else
    % T1 is a tel object
    if isscalar(T2)
        T = tel_scalar_mul(T2, T1);
    else
        T = tel_mat_vec_mul(T1, T2);
    end
end