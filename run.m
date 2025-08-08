clear;
%clc;
warning off;
addpath(genpath('./'));

%DBDIR = 'D:/data/';
DBDIR = 'dataset/';
%% dataset

i= 1;
DataName{i} = 'Caltech101-20'; i = i + 1;
dbNum = length(DataName);


for dsi = 1:dbNum

    clear X gt Y;
    dataName = DataName{dsi};
    dbfilename = sprintf('%s%s.mat',DBDIR,dataName);
    load(dbfilename);

    Y = gt;    
    k = length(unique(Y));   
    num_view = length(X);
    for v=1:num_view
        X{v} = X{v}';
    end
        
    %% para setting
    anchor = [1,2]*k ;
    d = [1,2,4]*k ;   
    lamb = [0.0001,0.001, 0.01,0.1,1];  
 
    if contains(dataName,'Caltech101-20') 
        anchor=[1]*k ;
   	    d = [4]*k ;    
        lamb = [1];         
    end 
    
   [index] = SLWSE(X,Y,d(1),anchor(1), lamb(1)); % X,Y,lambda,d,numanchor

    res = Clustering8Measure(Y, index); % ACC nmi AR Fscore Purity  Precision Recall  
    str = sprintf('db:%s\t  ACC:%.4f nmi:%.4f AR:%.4f Fscore:%.4f Purity:%.4f  Precision:%.4f Recall:%.4f \n',...
        dataName,  res(1), res(2), res(3), res(4), res(5), res(6), res(7));
    fprintf(str);
    clear index;   
end


