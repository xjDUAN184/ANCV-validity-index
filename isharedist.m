function [dist_result,number,pairwise] = isharedist(data,p,q,p_knn,q_knn,output,Eu_dist,pp)
    % 输入参数：数据集data  p点序号  q点序号  p点k-近邻序号  q点k-近邻序号  聚类结果
    % 输出结果：p点和q点的非共享近邻距离（序列）
    p_point = data(p,:);
    q_point = data(q,:);
    
%     % 去掉p点的k近邻中和p处于不同簇的点
%     p_nohome = find(output(p_knn)~=output(p));
%     p_knn(p_nohome) = [];
%     
%     % 去掉q点的k近邻中和q处于不同簇的点
%     q_nohome = find(output(q_knn)~=output(q));
%     q_knn(q_nohome) = [];
    
    snn = intersect(p_knn,q_knn);
    
    p_q_sharedneibor = setdiff([p_knn,p],[snn,q]);  % p相对于q的非共享近邻点
    q_p_sharedneibor = setdiff([q_knn,q],[snn,p]);  % q相对于p的非共享近邻点
    
    % 画图
%     if pp == 1
%         for i = 1:length(p_q_sharedneibor)
%             a = p_q_sharedneibor(i);
%             plot(data(a,1),data(a,2),'^','color','b','MarkerSize',3);
%         end
%         for j = 1:length(q_p_sharedneibor)
%             b = q_p_sharedneibor(j);
%             plot(data(b,1),data(b,2),'^','color','b','MarkerSize',3);
%         end
%     else
%         for i = 1:length(p_q_sharedneibor)
%             a = p_q_sharedneibor(i);
%             plot(data(a,1),data(a,2),'^','color','b','MarkerSize',3);
%         end
%         for j = 1:length(q_p_sharedneibor)
%             b = q_p_sharedneibor(j);
%             plot(data(b,1),data(b,2),'^','color','b','MarkerSize',3);
%         end
%     end
    pairwise = [];pw = 1;
    for i = 1:length(p_q_sharedneibor)
        for j = 1:length(q_p_sharedneibor)
            if p_q_sharedneibor(i) < q_p_sharedneibor(j)
                pairwise(pw,1) = p_q_sharedneibor(i);
                pairwise(pw,2) = q_p_sharedneibor(j);
            else
                pairwise(pw,1) = q_p_sharedneibor(j);
                pairwise(pw,2) = p_q_sharedneibor(i);  
            end
            pw = pw + 1;
        end
    end
    distance = Eu_dist(p_q_sharedneibor,q_p_sharedneibor);
    dist_result = distance(:);
    number = length(dist_result);
end