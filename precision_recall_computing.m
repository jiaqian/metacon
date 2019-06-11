function [pre_overall,rec_overall] = precision_recall_computing(C,contigLabel,label)  
num_cluster=size(C,1);
contigs_in_clusters=cell(num_cluster,2);
for nu=1:num_cluster
    contigs_in_clusters{nu,1}=nu;
    contigs_in_clusters{nu,2}=find(label==nu);
end

recorder=cell(num_cluster,5);
numerator=0;
denominator_p=0;
for u=1:num_cluster
    recorder{u,1}=u;
    if size(contigs_in_clusters{u,2},2)==1
        recorder{u,2}=size(contigs_in_clusters{u,2},1);
        recorder{u,3}=contigs_in_clusters{u,2};
    else
    recorder{u,2}=size(contigs_in_clusters{u,2},2);
    recorder{u,3}=contigs_in_clusters{u,2}';
    end
    s=unique(contigLabel(recorder{u,3}));
    ar=[];
    for g=1:size(s,1)
        ar=[ar sum(contigLabel(recorder{u,3})==s(g))];
    end
    recorder{u,4}=ar;
    recorder{u,5}=max(recorder{u,4})/recorder{u,2};
    numerator=numerator+max(recorder{u,4});
    denominator_p=denominator_p+recorder{u,2};
end
pre_overall=numerator/denominator_p;


[species,~]=unique(contigLabel);
recor_recall=cell(size(species,1),2);
ma=0;
denominator_re=0;
recall=zeros(1,size(species,1));
for s=1:size(species,1)
    id=find(contigLabel==species(s,1));
    cont=[];
    for re=1:size(id,1)
        id_sig=id(re,1);
        for u=1:num_cluster
        t_c=find(recorder{u,3}==id_sig);
        if t_c>0
            cont=[cont u];
            break;
        end
        end
    end
    uu=unique(cont);
    num_re=[];
    for g=1:size(uu,2)
        num_re=[num_re sum(cont==uu(1,g))];
    end
    recor_recall{s,1}=uu;
    recor_recall{s,2}=num_re;
    max_c=max(num_re);
    ma=ma+max_c;
    re=max_c/sum(num_re);
    recall(1,s)=re;
    denominator_re=denominator_re+size(id,1);
end
rec_overall=ma/denominator_re;
end