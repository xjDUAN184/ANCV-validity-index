function [ IDX ] = NTHC_clustering( X, K )
% NTHC clustering.
%   IDX = NTHC_clustering(X, K) partitions the points in the N-by-P data matrix X into K clusters and 
% returns an N-by-1 vector IDX containing the cluster indices of each point.
% 
%   input_args:
%       X : Data sets to be clustered (Rows of X correspond to points, columns correspond to variables)
%       K : The number of clusters
%
%   output_args:
%       IDX : The final clusters' label
%     
%   Demo for UCI data sets : 
%       dataSet = importdata('Sticks.mat');   
%       K = 3;
%       [ clusterLabel ] = NTHC_clustering( dataSet, K );

k = 10;
r = 0.2;
[N, d] = size(X);
cluster_label = zeros(N, 1);
antihub_label = zeros(N, 1);

[ nearest_neighbors_matrix, Nx_of_point, ~ ] = nei_Nx_sharenei( squareform( pdist(X, 'euclidean') ), k, N, 1, 0 );
[B, ~]=sort(Nx_of_point, 1, 'ASCEND');
n_antihub = floor(N*r);
thres_anti = B(n_antihub);


%% anti points
body_point = [];
anti_point = [];
temp_point = [];
for i=1:N  
    if Nx_of_point(i)>thres_anti
        body_point = [body_point i];
    else if Nx_of_point(i)==thres_anti
            temp_point = [temp_point i];
        else
            anti_point = [anti_point i];
        end
    end  
end

distanceMatrix = squareform( pdist(X, 'euclidean') );
if numel(temp_point)<=1
    anti_point = [anti_point temp_point];
else
    dist_tp = zeros(numel(temp_point), 1);
    for j=1:numel(temp_point)
        p = temp_point(j);
        
        nearest_nei = nearest_neighbors_matrix(p).nearest_nei;        
        dist_tp(j) = sum(distanceMatrix(p, nearest_nei))/numel(nearest_nei);
    end     
    while(numel(anti_point)<n_antihub)
        max_i = find(dist_tp==max(max(dist_tp)));
        
        anti_point = [anti_point temp_point(max_i)];
        temp_point(max_i)=[];
        dist_tp(max_i)=[];
    end
    body_point = [body_point temp_point];  
end



%% body point
dataSet_body = X(body_point, :);
N2 = size(dataSet_body, 1);
distance_matrix_body=squareform( pdist(dataSet_body, 'euclidean') );
for j=1:N2
    distance_matrix_body(j,j)=Inf;
end
[ nearest_neighbors_matrix_body, ~, share_nei_matrix_body ] = nei_Nx_sharenei( distance_matrix_body, k, N2, 0, 1 );


adjacency_matrix = zeros(N2, N2);
for i=1:N2
    nearest_nei = find(distance_matrix_body(i, :)==min(min(distance_matrix_body(i, :))) );
    for j=1:numel(nearest_nei)
        n = nearest_nei(j);
        adjacency_matrix(i, n) = 1;
    end   
end


digraph_1 = digraph(adjacency_matrix);
bins = conncomp(digraph_1, 'Type', 'strong');
for i=1:N2
    be_p = i;  
    en_ps = find(adjacency_matrix(be_p, :)==1);
    for j=1:numel(en_ps)
        en_p = en_ps(j);
        
        if bins(be_p)~=bins(en_p)
            if be_p<en_p
                jj = en_p;
                ii = be_p;
            else
                jj = be_p;
                ii = en_p;
            end
            loc = ((jj-2)*(jj-1))/2 + ii;
            share_nei = share_nei_matrix_body(loc).share_nei;
            n_nearest_nei1 = numel(nearest_neighbors_matrix_body(be_p).nearest_nei);
            n_nearest_nei2 = numel(nearest_neighbors_matrix_body(en_p).nearest_nei);
            if 2*numel(share_nei)/(n_nearest_nei1+n_nearest_nei2) <0.65
                adjacency_matrix(be_p, en_p)=0;     
            end   
        end
    end
end
            
 

digraph_2 = digraph(adjacency_matrix);
bins = conncomp(digraph_2, 'Type', 'weak');
isolate_po = [];
for i=1:max(bins)
    conn = find(bins==i);  
    sub_c(i).points = conn;
end
for i=1:size(sub_c, 2)
    subcluster = sub_c(i).points;
    cluster_label(subcluster) = i;   
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

%% Merge Subcluster
n_num = size(sub_c, 2);
dist_bet_subc = Inf(n_num, n_num);
for i=1:n_num
    for j=i+1:n_num
        dist_bet_subc(i,j) = compute_dist_bet_subclusters( sub_c(i).points, sub_c(j).points, distance_matrix_body, nearest_neighbors_matrix_body);
    end
end
while(n_num>K)
    min_adist = min(min(dist_bet_subc));  
    [ci, cj] = find(dist_bet_subc == min_adist);  
    
    if numel(ci)>=1
        min_dist2 = Inf;
        for ii = 1:numel(ci)
            cii = ci(ii);
            cjj = cj(ii);
            
            a_dist  = sum(sum( distance_matrix_body(sub_c(cii).points, sub_c(cjj).points) )) / ( numel(sub_c(cii).points)*numel(sub_c(cjj).points) );
            if a_dist<min_dist2
                min_dist2 = a_dist;
                subc1_tom = cii;
                subc2_tom = cjj;
            end
        end      
    else
        subc1_tom = ci(1);
        subc2_tom = cj(1);
    end
    
    sub_c(subc1_tom).points = [sub_c(subc1_tom).points, sub_c(subc2_tom).points];
    sub_c(subc2_tom) = [];
    n_num = size(sub_c, 2);

    dist_bet_subc(subc2_tom, :) = [];
    dist_bet_subc(:, subc2_tom) = [];
    for i = 1:n_num
        for j = i+1:n_num
            if i==subc1_tom || j==subc1_tom         
                dist_bet_subc(i,j) = compute_dist_bet_subclusters( sub_c(i).points, sub_c(j).points, distance_matrix_body, nearest_neighbors_matrix_body);               
            end
        end
    end    
end

cluster(K).points = [];
for i=1:K
    cluster(i).points = body_point(sub_c(i).points);
end


%% assign antihub 
n_antihub = numel(anti_point);
while n_antihub>0   
    md_antihub = Inf(n_antihub, 1);
    l_antihub = zeros(n_antihub, 1);    
    for i=1:n_antihub
        ap = anti_point(i);       
        for j=1:K
            cur_dist = min(min( distanceMatrix(ap, cluster(j).points) ));
            if cur_dist<md_antihub(i)
                md_antihub(i) = cur_dist;
                l_antihub(i) = j;
            end
            
        end
    end   
    res_i = find(md_antihub == min(min(md_antihub)) );
    if numel(res_i)>1
        min_adist = Inf;
        for i =1:numel(res_i)
            a_dist = sum(sum(distanceMatrix(anti_point(res_i(i)), cluster(l_antihub(res_i(i))).points) )) / numel( cluster(l_antihub(res_i(i))).points );
            if a_dist<min_adist
                min_adist = a_dist;
                gi = res_i(i);
            end
        end
        res_i = gi;
    end    
    gi = l_antihub(res_i);
    cluster( gi ).points = [ cluster( gi ).points, anti_point(res_i) ];
    anti_point(res_i) = [];
    n_antihub = numel(anti_point);
end
for i=1:K
    points = cluster(i).points;   
    for j=1:numel(points)
        p = points(j);
        cluster_label(p) = i;
    end
end
IDX = cluster_label;


function [ dist ] = compute_dist_bet_subclusters( subc1, subc2, distance_matrix, nearest_neighbors_matrix)
    
    p_in_subc1 = [];
    p_in_subc2 = [];
    for i_p1 = 1:numel(subc1)
        for i_p2 = 1:numel(subc2)           
            p1 = subc1(i_p1);
            p2 = subc2(i_p2);           
            if ismember(p1, nearest_neighbors_matrix(p2).nearest_nei) && ismember(p2, nearest_neighbors_matrix(p1).nearest_nei)
                p_in_subc1 = [p_in_subc1, p1];
                p_in_subc2 = [p_in_subc2, p2];
            end
            
        end
    end
    
    if numel(p_in_subc1)<=0
        dist = Inf;
        return;
    end   
    
    p_in_subc1 = unique(p_in_subc1);
    p_in_subc2 = unique(p_in_subc2);
    
    p_in_subc1_add = [];
    p_in_subc2_add = [];
    for i_p1 = 1:numel(subc1)
        p1 = subc1(i_p1);             
        if ~ismember(p1, p_in_subc1)          
            for iii = 1:numel(p_in_subc1)
                if ismember(p1, nearest_neighbors_matrix(p_in_subc1(iii)).nearest_nei)
                    p_in_subc1_add = [p_in_subc1_add, p1];
                    break;
                end
            end
        end
    end
    for i_p2 = 1:numel(subc2)
        p2 = subc2(i_p2);               
        if ~ismember(p2, p_in_subc2)           
            for iii = 1:numel(p_in_subc2)
                if ismember(p2, nearest_neighbors_matrix(p_in_subc2(iii)).nearest_nei)
                    p_in_subc2_add = [p_in_subc2_add, p2];
                    break;
                end
            end
        end
    end       
    p_in_subc1 = [p_in_subc1, p_in_subc1_add]; 
    p_in_subc2 = [p_in_subc2, p_in_subc2_add];
    p_in_subc1 = unique(p_in_subc1);
    p_in_subc2 = unique(p_in_subc2);    
    acc_dist = 0;
    for iii = 1:numel(p_in_subc1)
        for jjj = 1:numel(p_in_subc2)
            acc_dist = acc_dist + distance_matrix(p_in_subc1(iii), p_in_subc2(jjj));
        end
    end   
    dist = acc_dist / ( numel(p_in_subc1)*numel(p_in_subc2) );       
end

end

