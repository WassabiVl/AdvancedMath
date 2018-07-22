function [x] = gauss2(A,b)
    if size(b,2)>1
        b = b';
    end
    n = length(b);
    A1 = zeros(size(A));
    b1 = zeros(size(b));
    for k=1:(n-1)
        for i=(k+1):n
            if A(k,k) ~= 0
                lambda = A(i,k)/A(k,k);
                A(i,(k+1:n)) = A(i,(k+1:n)) - lambda*A(k,(k+1:n));
                b(i) = b(i) - lambda*b(k);
            else
                %error('Division by zero');
                A1((2:n-1),:) = A((2:n-1),:);
                A1(1,:) = A(n,:);
                A1(n,:) = A(1,:);
                A = A1;
                b1((2:n-1)) = b((2:n-1));
                b1(1) = b(n);
                b1(n) = b(1);
                b = b1;
                lambda = A(i,k)/A(k,k);
                A(i,(k+1:n)) = A(i,(k+1:n)) - lambda*A(k,(k+1:n));
                b(i) = b(i) - lambda*b(k);
            end
        end
    end
    for k=n:-1:1
        b(k) = (b(k) - A(k,(k+1):n)*b(k+1:n))/A(k,k);
    end
    x = b;
end