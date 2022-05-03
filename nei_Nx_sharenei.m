function [ nearest_neighbors_matrix, n_rknn, share_nei_matrix ] = nei_Nx_sharenei( distance_matrix, num_of_neighbor, N, is_need2, is_need3 )
% nearest_neighbors


for j=1:N
    distance_matrix(j,j)=Inf;
end
%% nearest_neighbors
nearest_neighbors_matrix(N).nearest_nei = [];
n_rknn = zeros(N, 1);

[B, ~]=sort(distance_matrix, 2, 'ASCEND');
thres_nei = B(:, num_of_neighbor);
for i=1:N
    %find nearest neighbors
    neighbors = find( distance_matrix(i, :)<=thres_nei(i) );
    nearest_neighbors_matrix(i).nearest_nei = neighbors;
    
    %compute n_rknn(Nx)
    if is_need2
        for j=1:numel(neighbors)
            nei = neighbors(j);
            n_rknn(nei) = n_rknn(nei)+1;
        end
    end
end



%% share_nei_matrix 
if is_need3
    share_nei_matrix( (N*(N-1))/2 ).share_nei = [];           
    for j=2:N
        for i=1:(j-1)

            l_nei1 = nearest_neighbors_matrix(j).nearest_nei;
            l_nei2 = nearest_neighbors_matrix(i).nearest_nei;

            loc = ((j-2)*(j-1))/2 + i;
            share_nei_matrix(loc).share_nei = intersect(l_nei1, l_nei2); 
        end
    end
else
    share_nei_matrix = [];
end




end

