clear;
%clc;
warning off;
addpath(genpath('./'));

%DBDIR = 'H:/data/';
DBDIR = 'dataset/';
%% dataset

i= 1;

% DataName{i} = 'Wiki'; i = i + 1;
% DataName{i} = 'COIL20'; i = i + 1;
DataName{i} = 'Caltech101-20'; i = i + 1;
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
    anchor_list = [1,2,4]*k ;
    d_list = [1,2,4]*k ;
    lamb_list = [0.0001,0.001, 0.01,0.1,1];
    %%%%% for Caltech101-20 %%%%%%%%%%%%%
    d = d_list(2);
    m = anchor_list(1);
    lambda = lamb_list(4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [index] = SLWSE(X,Y,d,m, lambda); % X,Y,lambda,d,numanchor
    res = Clustering8Measure(Y, index); % ACC nmi AR Fscore Purity  Precision Recall
    str = sprintf('db:%s\t  ACC:%.4f nmi:%.4f AR:%.4f Fscore:%.4f Purity:%.4f  Precision:%.4f Recall:%.4f \n',...
        dataName,  res(1), res(2), res(3), res(4), res(5), res(6), res(7));
    fprintf(str);
end


