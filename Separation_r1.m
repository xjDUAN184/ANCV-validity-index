function sep = Separation_r1(data_new,output,Dataset_MST,knn,Eu_dist,SNN_thr)
    k_value = length(unique(output));
    sep = zeros(k_value,k_value);
    % 首先判断两个簇的MST之间是否有边相连
    for i = 1:k_value-1
        idx = find(output==i);
        num = 0;
        clu_mst = Dataset_MST(idx,:);
        n = length(idx);
        for j = i+1:k_value
            juli = []; point_wise_out = [];
            for L = 1:n
                v1 = find(clu_mst(L,:)~=0);
                if isempty(v1)==1
                    continue;
                elseif isempty(find(output(v1)==j))~=1
                    v2 = find(output(v1)==j);
                    for J = 1:length(v2)
                        p = idx(L);
                        q = v1(v2(J));
                        
                        sn = length(intersect(knn(p,:),knn(q,:)));
%                         if sn < SNN_thr
%                             plot(data_new(p, 1), data_new(p, 2), '.','color','r','MarkerSize',15);
%                             plot(data_new(q, 1), data_new(q, 2), '.','color','r','MarkerSize',15);
%                             line([data_new(p,1),data_new(q,1)],[data_new(p,2),data_new(q,2)],'color','r','LineWidth',1);
%                         end
                        
                        [pq_dist,~,Pair_point] = isharedist(data_new,p,q,knn(p,:),knn(q,:),output,Eu_dist,2);
                        point_wise_out = [point_wise_out;Pair_point]; %去重前的点对
                        juli = [juli;pq_dist];num = num + 1;
%                         juli(sn) = sum(sum(pq_dist))/num;
%                         sn = sn + 1;
                    end
                end
            end
            distance = [];
            P_Wmatrix_out = unique(point_wise_out,'row','stable'); %去重后的点对
            for row_num = 1:size(P_Wmatrix_out,1)
                a = P_Wmatrix_out(row_num,1);b = P_Wmatrix_out(row_num,2);
                distance(row_num) = Eu_dist(a,b); % 点对距离矩阵
            end
            if length(juli) ~= 0 
                sep(i,j) = mean(distance);
                sep(j,i) = sep(i,j);
            end
        end
    end
end