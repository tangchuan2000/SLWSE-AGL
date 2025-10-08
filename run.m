clear;
%clc;
warning off;
addpath(genpath('./'));

DBDIR = 'H:/data/';
%DBDIR = 'dataset/';
%% dataset

i= 1;

% DataName{i} = 'Wiki'; i = i + 1;
DataName{i} = 'COIL20'; i = i + 1;
% DataName{i} = 'Caltech101-20'; i = i + 1;
% DataName{i} = 'NUSWIDEOBJ'; i = i + 1;
% DataName{i} = 'YouTubeFace10_4Views'; i = i + 1;
% DataName{i} = 'AwA'; i = i + 1;
% DataName{i} = 'YouTubeFace50_4Views'; i = i + 1;
dbNum = length(DataName);
for dsi = 1:dbNum
    
    clear X gt Y;
    dataName = DataName{dsi};
    dbfilename = sprintf('%s%s.mat',DBDIR,dataName);
    load(dbfilename);
    
    Y = gt;
    k = length(unique(Y));
    num_view = length(X);
    
    %% para setting
    anchor = [1,2,4]*k ;
    d = [1,2,4]*k ;
    lamb = [0.0001,0.001, 0.01,0.1,1];

    maxAcc = 0;
    maxStr = '';
    for i_m = 1:length(anchor)
        for i_d = 1:length(d)
            for i_lamb = 1:length(lamb)                
                [index] = SLWSE(X,Y,d(i_d),anchor(i_m), lamb(i_lamb)); % X,Y,lambda,d,numanchor
                res = Clustering8Measure(Y, index); % ACC nmi AR Fscore Purity  Precision Recall
                str = sprintf('db:%s\t  ACC:%.4f nmi:%.4f AR:%.4f Fscore:%.4f Purity:%.4f  Precision:%.4f Recall:%.4f \n',...
                    dataName,  res(1), res(2), res(3), res(4), res(5), res(6), res(7));
                fprintf(str);
                clear index;
                if (maxAcc < res(1))
                    maxAcc = res(1);
                    maxStr = str;
                end
            end
        end
    end
    fprintf('Max:%s', maxStr);
    clear X Y k;
end


