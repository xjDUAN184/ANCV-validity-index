function [Sep,Com,sep_clu,result] = validity_index(data,label,K,SNN_thr,Eu_dist)
% function result = validity_index(data,label,K,SNN_thr)
    k_value = max(unique(label));
    [n,m] = size(data); 
    Eu_dist = pdist2(data,data);
    [~,dist_index] = sort(Eu_dist,2);
    knn = dist_index(:,2:K+1); %每个点的K近邻
    Dataset_MST = makeMST_r1(Eu_dist); % 整个数据集的MST
    %% Compute the Compactness
    com = Compactness_r9(data,label,Dataset_MST,SNN_thr,knn,Eu_dist);
    Com = mean(com);
%     Com = sum(com); 
    %% Compute the Separation
    sep_clu = Separation_r1(data,label,Dataset_MST,knn,Eu_dist);
    
    for i = 1:k_value
        sep_c = sep_clu(i,find(sep_clu(i,:)~=0));
        index(i) = (max(sep_c) - com(i)) / max(max(sep_c),com(i));
    end
    result = mean(index);
    
    sep = tril(sep_clu);
    sep = sep(sep>0);
    Sep = mean(sep);
%     Sep = sum(sep);
end