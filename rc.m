function d = rc(A, B, m0)
    len = length(A) + length(B) - 1;
    d = zeros(2,len);
    d(1,:) = (A-m0(1)) ;
    d(2,:) = (B-m0(2));
end

