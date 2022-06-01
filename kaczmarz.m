function [ x ] = kaczmarz( A,b,iterations,lambd,shuff,enforceReal,enforcePositive )
% 20220520 lilei
% Ax=b,
% shuff：whether to use the randomized Kaczmarz

% initialization of the variable
% A has n rows and m columns（n * m）
% b must be 1 * m, m * 1 will cause a wrong answer

[N, M] = size(A);

x = complex(zeros(N,1)); 

rowIndexCycle = 1:M;

% calculate the energy of each frequency component
% 
% Calculate the inner product of itself and sum by columnxian
energy = rowEnergy(A);

% may use a randomized Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
end

% estimate regularization parameter

for l = 1:iterations
    for m = 1:M
        k = rowIndexCycle(m);
        
        if energy(k) > 0
            tmp = A(:,k).'*x;
            beta = (b(k) - tmp ) / (energy(k)^2 );
            x = x + beta*conj(A(:,k));
        end
    end
    
    
    if enforceReal && ~isreal(x)
        x = complex(real(x),0);
    end
    
    if enforcePositive
        x(real(x) < 0) = 0;
    end
end


 function [ energy ] = rowEnergy(A)		
 % Calculate the norm of each row of the input		
 % This is an image of the energy stored by each frequency component		
 % but the scaling is unclear.		
 		
 energy = sqrt(sum(abs(A.*A),1));		
 		
 end		


end
