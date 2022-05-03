function com = Compactness_r9(data_new,output,Dataset_MST,SNNthr,knn,Eu_dist)
    k_value = length(unique(output));
    com = zeros(k_value,1);
    for i = 1:k_value
        pairwise = 0; % 满足有相连边且共享近邻个数小于阈值的点对数
        tiaojian = [];
        idx = find(output==i);
        distance = []; point_wise_in = [];
        n = length(idx);
        clu_mst = Dataset_MST;
        clu_mst_tril = tril(clu_mst);
        if n == 1
            com(i) = 0;
            continue;
        end
        for ii = 1:n
            mst_neibor_clu = find(clu_mst_tril(idx(ii),:)~=0);
            no_home1 = find(output(mst_neibor_clu)~=i);
            mst_neibor_clu(no_home1) = []; %去掉该点在mst上相连的不属于同一个簇的点
            i_knn = knn(idx(ii),:);
%             
%             no_homo2 = find(output(i_knn)~=i);
%             i_knn(no_homo2) = []; % 去掉该点的k近邻点中不属于i簇的点
            len = length(mst_neibor_clu);
            for jj = 1:len
                near = mst_neibor_clu(jj);                
                j_knn = knn(near,:);
%                 no_homo3 = find(output(j_knn)~=i);
%                 j_knn(no_homo3) = []; % 去掉该点的k近邻点中不属于i簇的点
                
                i_point = i_knn;
                j_point = j_knn;
                
                sn = length(intersect(i_point,j_point));
                
                if sn < SNNthr
                    pairwise = pairwise + 1;
                    tiaojian(pairwise,1) = idx(ii);tiaojian(pairwise,2) = near;
%                     line([data_new(idx(ii),1),data_new(near,1)],[data_new(idx(ii),2),data_new(near,2)],'color','b','LineWidth',1);
                    [sdist,num,Pair_point] = isharedist(data_new,idx(ii),near,i_knn,j_knn,output,Eu_dist,1);
                    point_wise_in = [point_wise_in;Pair_point]; %去重前的点对
%                     distance = [distance;sdist];
                end
            end
        end
        
        P_Wmatrix_in = unique(point_wise_in,'row','stable'); %去重后的点对
        for row_num = 1:size(P_Wmatrix_in,1)
            a = P_Wmatrix_in(row_num,1);b = P_Wmatrix_in(row_num,2);
            distance(row_num) = Eu_dist(a,b); % 点对距离矩阵
        end
        
        CLU_MST = clu_mst_tril(idx,idx);
        if length(find(CLU_MST>0))~=0
            special = mean(CLU_MST(find(CLU_MST~=0)));
        else
            special = mean(mean(Eu_dist(idx,idx)));
        end
        if pairwise>=1
%             com(i) = special;
            com(i) = mean(distance);
        else
            if isempty(find(clu_mst_tril~=0)) ~= 1
                com(i) = special; % 非孤立点，但是较为紧凑的簇
            else
                com(i) = 0; % 孤立点
            end
        end
%         if com(i) < special
%             com(i) = special;
%         end
    end
    com(find(com==0)) = max(com);
end