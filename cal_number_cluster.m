function [num_cluster]=cal_number_cluster(TT,C,D)

kk=zeros(1,3);
N_upp=size(TT,1);
kArr = 10:10:N_upp; 
for loop=1:3
for kIdx = 1: length(kArr)
   
candK = kArr(kIdx);
options = []; 
options.start = 1;
options.repeat = 10; options.blockLen = 1;
[~,Wpre,~] = fastkmeans(TT',candK,options);
if size(Wpre,2) < C*candK
    k=size(Wpre,2);
    k = size(Wpre,2)*D;
    kk(1,loop)=k;
    break;
end    
end
end
num_cluster=ceil(mean(kk));
end