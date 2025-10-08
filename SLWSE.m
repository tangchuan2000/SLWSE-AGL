function [y]=SLWSE(X ,Y, d ,numanchor, lambda)

    m = numanchor ;
    numclass = length (unique (Y));
    numview = length (X);
    numsample = size (Y ,1);

    w = cell(1, numview);
    ww_vector = ones(1,numsample) / numview;
    for i = 1:numview
        w{i} = ww_vector;
    end

    P =cell(numview ,1);
    A = zeros(d ,m);
    Z = zeros(m ,numsample);
    Q = zeros(numclass, m);

    Q(1,:) =  1;% first row set to 1
    parfor i = 1 :numview
        di = size (X{i},1);
        P{i}=zeros (di ,d);
        X{i}=mapstd (X{i} ,0 ,1);
    end

    Z(m,:) = 1;
    flag = 1 ;
    iter = 0 ;
    rng(5489,'twister');
    maxIter = 30;
    while flag
        iter = iter +1 ;

        %update P{v}
        AZ = A * Z ;
        parfor v = 1 :numview
            di = size(X{v},1);
            C = X{v} .* repmat(w{v}.^2,di,1) * AZ' ;
            [U ,~,V ] = svd(C ,'econ');
            P{v} = U * V' ;
        end
        clear AZ;

        %update A
        part1 = 0 ;
        termpart1 = cell(numview ,1);
        parfor v = 1 :numview
            termpart1{v} = P{v}' *  X{v} .*  repmat(w{v}.^2,d,1) * Z' ;
        end
        for v = 1 : numview
            part1 = part1 + termpart1{v}  ;
        end
        [Unew ,~,Vnew ] = svd(part1 ,'econ');
        clear part1;

        A = Unew * Vnew' ;

        %update Z
        C_cell = cell(numview , 1);
        for v = 1 : numview
            C_cell{v} = P{v} * A ;
        end
        options = optimset('Algorithm' ,'interior-point-convex' ,'Display' ,'off');
        parfor ji = 1 : numsample
            ff = 0 ;
            G = 0;
            for v = 1 : numview
                ff = ff - w{v}(ji).^2 * X{v}(:,ji)' * C_cell{v} ;
                G = G + w{v}(ji).^2;
            end
            GG = G * eye(m);
            GG =(GG +GG')/2 ;

            QQZ = Q' * Q * Z(:,ji);
            di = QQZ * (Z(:,ji)' * QQZ)^(-0.5);
            ff = ff - lambda * di';
            Z(:,ji) = quadprog(GG,ff',[],[],ones(1 ,m),1 ,zeros(m ,1),ones(m ,1),[],options);
        end

        %update Q
        Q = updateQ(Z', Q');
        Q = Q';

        %update w
        w = updateW(X, P, A, Z);

        term1 = 0 ;
        termobj = cell(numview , 1);
        parfor v = 1 : numview
            di = size(X{v}, 1);
            termobj{v} = norm((X{v} - P{v} * A * Z)  .* repmat(w{v},di,1), 'fro');
        end

        for v = 1 : numview
            term1 = term1 +  termobj{v};
        end

        obj(iter) = term1 -  lambda * (trace(Q * Z * Z' * Q'))^(0.5);

        if((iter >1 && abs(obj(iter) - obj(iter -1)) / obj(iter -1) < 0.0001)|| iter == maxIter)
            flag = 0 ;
            F = Z' * Q';
            [~,y]=max(F,[],2);
        end
    end
end



