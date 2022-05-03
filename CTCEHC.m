function [label, link] = CTCEHC(data, k, minNum, knn, T)

% MST_IN_CLUSTERING
%    generate a MST in the complete graph that created by dataset.
%    

dist_tmp = get_euclid_dist(data');
fprintf("Got dist...\n");
graph = get_complete_graph(dist_tmp);
fprintf("Got graph...\n");
[link, dist] = get_mst_in_complete_graph(graph);
fprintf("Got MST...\n");

if (nargin == 1)
label = bfs_merged_with_edges(link, dist, graph);
end
if (nargin == 2)
label = bfs_merged_with_edges(link, dist, graph, k);
end
if (nargin == 3)
label = bfs_merged_with_edges(link, dist, graph, k, minNum);
end
if (nargin == 4)
label = bfs_merged_with_edges(link, dist, graph, k, minNum, knn);
end
if (nargin == 5)
label = bfs_merged_with_edges(link, dist, graph, k, minNum, knn, T);
end
fprintf("Finished...\n");
end
