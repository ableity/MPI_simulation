function [ x ] = kaczmarzReg2(A,b,iterations,lambd,corr_lambd,shuff,enforceReal,enforcePositive)
% A has n rows and m columns（n * m）
% b must be 1 * m, m * 1 will cause a wrong answer
 [m,n] = size(A);
    energy = rowEnergy(A);
    lambdZero = sum(energy.^2)/m;
    lambdIter = lambd*lambdZero;


    A_plus = ((lambdIter)^0.5)*ones(1,n);
    A_plus = diag(A_plus);

    A = [A;A_plus];
    x = kaczmarz( A,b,iterations,lambd,corr_lambd,shuff,enforceReal,enforcePositive );
    x = x(1:m);
    
    function [ energy ] = rowEnergy(A)
        % Calculate the norm of each row of the input
        % This is an image of the energy stored by each frequency component
        % but the scaling is unclear.

        energy = sqrt(sum(abs(A.*A),1));
    end

end
