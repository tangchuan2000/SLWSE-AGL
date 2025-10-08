function [w] = updateW(X, P, A, Z)
% Input
% X   d * n
    
    numview =length(X); 
    numsample =size(X{1}, 2); 
    w = cell(1, numview);
    Y = cell(numview,1);
    YY = cell(numview,1);
    gama = 0.2;
    parfor v = 1:numview               
        Y{v} = (X{v}-P{v}*A *Z); 
        YY{v} = Y{v} .* Y{v};
    end   

    W_n = cell(numview, 1);
    for v = 1:numview 
       W_n{v} = zeros(1, numsample);
    end   

    for n = 1:numsample         
        sumOfVYcol = 0;
        for v = 1:numview   
            sumOfYcol = sum(YY{v}(:, n));
            sumOfYcol = 1 / sumOfYcol;
            sumOfVYcol = sumOfVYcol + sumOfYcol;
        end  
        for v = 1:numview  
            sumOfYcol2 = sum(YY{v}(:, n));
            W_n{v}(1, n) = (1 / (sumOfYcol2 * sumOfVYcol)) ;
        end   
    end 
    Wpow = zeros(numview, numsample);
    for v = 1:numview  
        Wpow(v,:) = W_n{v};
    end
    %power-normalized    
    Wpow = (Wpow ).^gama;  
    Wpn  = Wpow ./ max(sum(Wpow,1), 0.00000001);
    for v = 1:numview 
      w{v} = Wpn(v,:); 
    end
end