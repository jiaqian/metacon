%% load data set
[Head,ContigData]=fastaread('data/sharon/contigs_sharon_labeled.fasta');
V= csvread('data/sharon/label.txt');
contigLabel=V;
%contigId_ref={};
contigId_head={};
for oo=1:size(Head,2)
    str=Head{1,oo};
    contigId_head=[contigId_head;cell2mat(regexp(str,'[0-9]+\.+[0-9]|[0-9]','match'))];
end
contigId_ref=contigId_head;

%%
k=4;
v_alpha=['A','G','T','C'];
k_mers=permn(v_alpha,k);
num_contig=size(ContigData,2);
mers=[];
while ~isempty(k_mers)
    f=k_mers(1,:);
    mers=[mers; f];
    f_=seqrcomplement(f);
    k_mers(1,:)=[];
    k_mers(strmatch(f_,k_mers),:)=[];
end
num_mers=size(mers,1);
feature_matrix=zeros(num_contig,num_mers);
size_arr=zeros(1,num_contig);
for j=1:num_contig
l= length(ContigData{1,j});
size_arr(1,j)=l;
f=ContigData{1,j};
for h=1:num_mers
    mer=mers(h,:);
    mer_=seqrcomplement(mer);
    frq=sum(strfind(f,mer)>0)+sum(strfind(f,mer_)>0);
    feature_matrix(j,h)=frq;
end
end
feature_matrix = feature_matrix + 1;
FF=feature_matrix;
for er=1:size(FF,1)
    le=size_arr(1,er);
    FF(er,:)=FF(er,:)./le;
end
%% coverage matrix count
Y=csvread("data/sharon/covMat_inputtableR.csv"); 
Y = Y + 1e-2;
Y = Y ./ repmat(sum(Y,1),size(Y, 1),1);
Y=  Y ./ repmat(sum(Y, 2),1,size(Y, 2));

%% estimate the number of clusters
TT=[FF Y];
num_cluster=cal_number_cluster(TT,0.8,0.4);
display(['num_cluster is ',num2str(num_cluster)])
%% fisrt-stage: self-standarized k-mers statistics
cc=0;
R_=[];
for u=1:num_contig
     if size_arr(1,u)<=2000
         R_=[R_ u];
     else
         cc=cc+1;
        len=length(ContigData{1,u});
        prob_A=sum(strfind(ContigData{1,u},'A')>0)/len;% prob(A)
        prob_G=sum(strfind(ContigData{1,u},'G')>0)/len;% prob(G)
        prob_T=sum(strfind(ContigData{1,u},'T')>0)/len;% prob(T)
        prob_C=sum(strfind(ContigData{1,u},'C')>0)/len;%prob(C)
        for jj=1:num_mers
         str=mers(jj,:);
         num_a=sum(strfind(str,'A')>0);num_g=sum(strfind(str,'G')>0);
         num_t=sum(strfind(str,'T')>0);num_c=sum(strfind(str,'C')>0);
         coarse_prob=(prob_A^num_a)* (prob_G^num_g)* (prob_T^num_t)* (prob_C^num_c);
         average_frq=len*coarse_prob;
         variance_frq=len*coarse_prob*(1-coarse_prob);
         feature_matrix(u,jj)=abs(feature_matrix(u,jj)-average_frq)/variance_frq;
        end  
     end
end
%% fisrt-stage: k-medoids
left_part=setdiff([1:1:num_contig],R_);
F_part=feature_matrix(left_part,:);
F_part=F_part./repmat(sum(F_part,1),size(F_part,1),1);
X=[F_part Y(left_part,:)];
[belong,cluster_centroids]=kmedoids(X,num_cluster,'distance','sqeuclidean','Replicates',10);
[pre_first,rec_first] = precision_recall_computing(cluster_centroids,contigLabel(left_part),belong);

%% second-stage: assign the left contigs to the closest cluster centroid
B(left_part)=belong;
M2=Y(R_,:);
left_contigs=[feature_matrix(R_,:) M2];
n_left=size(left_contigs,1);
dist=zeros(n_left,num_cluster);
for b=1:n_left
    v=left_contigs(b,:);
    for bb=1:num_cluster
        v_=cluster_centroids(bb,:);
        dist(b,bb)=norm(v_-v,1);
    end  
end
[n_min,index_min]=min(dist,[],2);
B(R_)=index_min;
%% performance evaluation
[pre_overall,rec_overall] = precision_recall_computing(cluster_centroids,contigLabel,B);
