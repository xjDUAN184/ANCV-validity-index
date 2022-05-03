//
// Author: Hongren Lin
//

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <list>
#include <map>
#include <queue>

// #define DEBUG

#ifndef DEBUG
#include "mex.h"
#endif

#define EPS 1e-6
typedef long long ULL;

// 节点的度 <标记位, 度>
typedef std::pair<int, int> Degree;

// 度的比较函数
bool Degree_greater(const Degree &d1, const Degree &d2) {
    if (d1.second > d2.second) {
        return true;
    } else if (d1.second == d2.second) {
        return d1.first < d2.first;
    } else {
        return false;
    }
}

bool Degree_lesser(const Degree &d1, const Degree &d2) {
    if (d1.second < d2.second) {
        return true;
    } else if (d1.second == d2.second) {
        return d1.first < d2.first;
    } else {
        return false;
    }
}

// 浮点数转整数
static inline int ftoi(double x) { return floor(x + EPS); }

// 交换函数
template <class T>
static inline void _swap(T &a, T &b) {
    T t = a;
    a = b;
    b = t;
}

// 边结构体
typedef struct Edge {
    int s;
    int e;
    double w;
} Edge;

struct Edge_greater {
    bool operator()(const Edge &edge1, const Edge &edge2) {
        if (edge1.w < edge2.w) {
            return false;
        } else if (edge1.w > edge2.w) {
            return true;
        } else {
            if (edge1.s < edge2.s) {
                return false;
            } else if (edge1.s > edge2.s) {
                return true;
            } else {
                return edge1.e > edge2.e;
            }
        }
    }
};

struct Edge_lesser {
    bool operator()(const Edge &edge1, const Edge &edge2) {
        if (edge1.w < edge2.w) {
            return true;
        } else if (edge1.w > edge2.w) {
            return false;
        } else {
            if (edge1.s < edge2.s) {
                return true;
            } else if (edge1.s > edge2.s) {
                return false;
            } else {
                return edge1.e < edge2.e;
            }
        }
    }
};

struct Bridge {
    Edge edge;
    int c1, c2;
    int c1_count, c2_count;
    Bridge(Edge e, int _c1, int _c2, int _c1_count, int _c2_count)
        : edge(e), c1(_c1), c2(_c2), c1_count(_c1_count), c2_count(_c2_count) {}
};

struct Bridge_greater {
    bool operator()(const Bridge &bridge1, const Bridge &bridge2) {
        if (std::abs(bridge1.edge.w - bridge2.edge.w) < EPS) {
            if ((bridge1.c1_count < bridge2.c1_count &&
                 bridge1.c1_count < bridge2.c2_count) ||
                (bridge1.c2_count < bridge2.c1_count &&
                 bridge1.c2_count < bridge2.c2_count))
                return false;
            else
                return true;
        } else {
            return bridge1.edge.w > bridge2.edge.w;
        }
    }
};

typedef std::list<Edge> Graph;

// 全局最小生成树
Graph *mst = NULL;

// 全局完全图(矩阵表示)
double *complete_graph_dist = NULL;

// 全局节点度数
Degree *degree = NULL, *degree_sorted = NULL;

// 类别数量，总边数，总节点数
int cluster = 1, edge_N, node_N;

// 标签
int *label = NULL;

// 每个类别的数量
int *cluster_count = NULL;

// BFS标记数组
bool *vis = NULL;

// 每一个类别中的任意一个点的编号
int *any_in_cluster = NULL;

// 超参数
int knn;
double minNum, minDist, T;

// 更新被使用的节点的度
bool earse_degree(int cur) {
    for (int i = 0; i < node_N; i++) {
        if (degree_sorted[i].first == cur) {
            degree_sorted[i].first = -1;
            return true;
        }
    }
    return false;
}

// 找到当前未使用的最大的度的节点编号
Degree find_max_degree() {
    for (int i = 0; i < node_N; i++) {
        if (degree_sorted[i].first != -1) return degree_sorted[i];
    }
    return Degree(-1, 0);
}

// 从输入的边和距离数组构建最小生成树（输入已经是最小生成树）
void create_mst(double *edge, double *dist) {
    mst = new Graph[node_N];
#ifndef DEBUG
    if (mst == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    degree = new Degree[node_N];
#ifndef DEBUG
    if (degree == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(degree, 0, sizeof(Degree) * (node_N));
    int x, y;
    Edge e1, e2;
    for (int i = 0; i < edge_N; i++) {
        x = ftoi(edge[i]) - 1;
        y = ftoi(edge[i + edge_N]) - 1;
        e1.s = x;
        e1.e = y;
        e1.w = dist[i];
        mst[x].push_back(e1);

        e2.s = y;
        e2.e = x;
        e2.w = dist[i];
        mst[y].push_back(e2);

        degree[x].first = x;
        degree[x].second++;
        degree[y].first = y;
        degree[y].second++;
    }
}

// 找到s点的k临近的点的标号
std::vector<int> get_k_nearest_index(int k, int s) {
    std::priority_queue<std::pair<double, int> > Q;
    k++;
    for (int i = 0; i < node_N; i++) {
        if (Q.size() < k ||
            complete_graph_dist[s * node_N + i] < Q.top().first) {
            Q.push(std::make_pair(complete_graph_dist[s * node_N + i], i));
        }
        if (Q.size() > k) {
            Q.pop();
        }
    }
    std::vector<int> ans;
    while (Q.size() > 1) {
        ans.push_back(Q.top().second);
        Q.pop();
    }
    Q.pop();
    return ans;
}

// 两个点的k临近的标号的相同的个数
int same_node_in_k_nearest_index(int k, int s, int e) {
    std::vector<int> s_k = get_k_nearest_index(k, s);
    std::sort(s_k.begin(), s_k.end());
    std::vector<int> e_k = get_k_nearest_index(k, e);
    std::sort(e_k.begin(), e_k.end());

    int i = 0, j = 0;
    int res = 0;
    while (i < k && j < k) {
        if (s_k[i] == e_k[j]) {
            res++, i++, j++;
        } else if (s_k[i] < e_k[j]) {
            i++;
        } else {
            j++;
        }
    }
    return res;
}

// 首次聚类，单个节点（按度扩展）
void preliminary_clustering_one_step(int s, int k, double T) {
    std::queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int cur = q.front();
        q.pop();
        label[cur] = cluster;
        vis[cur] = 1;
        earse_degree(cur);
        for (auto i : mst[cur]) {
            if (vis[i.e] == 0) {
                bool max = true;
                double weight = i.w;
                for (auto j : mst[i.e]) {
                    if (weight >= j.w && j.e != cur) {
                        max = false;
                        break;
                    }
                }
                double m = same_node_in_k_nearest_index(k, cur, i.e);
                if (max) {
                    q.push(i.e);
                }
            }
        }
    }
}

// 首次聚类，所有节点
void preliminary_clustering(int k, double T) {
    label = new int[node_N];
#ifndef DEBUG
    if (label == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(label, 0, sizeof(int) * (node_N));
    vis = new bool[node_N];
#ifndef DEBUG
    if (vis == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(vis, 0, sizeof(bool) * (node_N));
    degree_sorted = new Degree[node_N];
#ifndef DEBUG
    if (degree_sorted == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memcpy(degree_sorted, degree, sizeof(Degree) * node_N);
    std::sort(degree_sorted, degree_sorted + node_N, Degree_greater);

    cluster = 1;
    while (true) {
        Degree cur = find_max_degree();
        if (cur.first == -1) break;
        preliminary_clustering_one_step(cur.first, k, T);
        cluster++;
    }
    delete[] vis;
}

// dfs寻找重心，无根树转有根
void find_barycentre_dfs(int node, int parent, int c, int *dp, int &minNode,
                         int &minBalance) {
    dp[node] = 1;
    int maxSubTree = 0;
    for (auto i : mst[node]) {
        int son = i.e;
        if (son != parent && label[son] == c) {
            find_barycentre_dfs(son, node, c, dp, minNode, minBalance);
            dp[node] += dp[son];
            maxSubTree = std::max(maxSubTree, dp[son]);
        }
    }
    maxSubTree = std::max(maxSubTree, cluster_count[c] - dp[node]);
    if (maxSubTree < minBalance) {
        minBalance = maxSubTree;
        minNode = node;
    }
}

// 计算类别c生成树的重心
int find_barycentre(const int &c) {
    int *dp = new int[node_N];
#ifndef DEBUG
    if (dp == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    int s = any_in_cluster[c];
    int minNode = -1;
    int minBalance = node_N + 1;
    find_barycentre_dfs(s, -1, c, dp, minNode, minBalance);
    if (dp != NULL) delete[] dp;
    return minNode;
}

double get_max_edge(int c) {
    int s = any_in_cluster[c];
    double max = 0.0;
    vis = new bool[node_N];
#ifndef DEBUG
    if (vis == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(vis, 0, sizeof(bool) * node_N);
    std::queue<int> Q;
    Q.push(s);
    while (!Q.empty()) {
        int cur = Q.front();
        Q.pop();
        vis[cur] = 1;
        for (auto i : mst[cur]) {
            if (!vis[i.e]) {
                if (label[i.e] == c) {
                    max = std::max(i.w, max);
                    Q.push(i.e);
                }
            }
        }
    }
    if (vis != NULL) delete[] vis;
    return max;
}

// 经过c1，c2类的两点的最短路
void shortest_path(int s, int e, int fa, int c1, int c2, double dist,
                   bool &find, double &ans, double &max) {
    if (s == e) {
        find = true;
        ans = dist;
        return;
    } else {
        for (auto i : mst[s]) {
            int v = i.e;
            if (v != fa && (label[v] == c1 || label[v] == c2)) {
                if (label[v] != label[s]) {
                    max = std::max(max, i.w);
                }
                shortest_path(v, e, s, c1, c2, dist + i.w, find, ans, max);
            }
        }
    }
}

// 最近类别 <编号, 距离>
struct Closest_cluster {
    int c;
    double d;
    Closest_cluster(int _c = -1, double _d = 0.0) : c(_c), d(_d){};
};

struct Closest_cluster_greater {
    bool operator()(const Closest_cluster &cc1, const Closest_cluster &cc2) {
        if (cc1.d > cc2.d) {
            return true;
        } else if (cc1.d < cc2.d) {
            return false;
        } else {
            return cc1.c > cc2.c;
        }
    }
};

// 寻找最近的类别（基于重心的测定距离）
Closest_cluster find_closest_cluster_in_tree(const int &c) {
    Closest_cluster res(-1, 0);
    std::priority_queue<Closest_cluster, std::vector<Closest_cluster>,
                        Closest_cluster_greater>
        dists, cache_dist;
    int s = find_barycentre(c);
    for (int i = 1; i < cluster; ++i) {
        if (i == c || cluster_count[i] == 0) continue;
        int e = find_barycentre(i);
        bool find = false;
        double dist = 0.0;
        double bridge = 0.0;
        double avg_max = get_max_edge(label[s]) + get_max_edge(label[e]);
        shortest_path(s, e, -1, c, i, 0, find, dist, bridge);

        if (find && bridge <= avg_max) {
            Closest_cluster cc(i, dist);
            dists.push(cc);
        } else if (find) {
            Closest_cluster cc(i, dist);
            cache_dist.push(cc);
        }
    }
    if (!dists.empty()) {
        res = dists.top();
        dists.pop();
    } else {
        if (!cache_dist.empty()) {
            res = cache_dist.top();
            cache_dist.pop();
        }
    }
    while (!dists.empty()) {
        dists.pop();
    }
    while (!cache_dist.empty()) {
        cache_dist.pop();
    }
    return res;
}

// 寻找最近的类别，基于最近距离
Closest_cluster find_closest_cluster(const int &c) {
    Closest_cluster res = Closest_cluster(-1, 0.0);
    int s = any_in_cluster[c];
    std::queue<int> Q;
    vis = new bool[node_N];
#ifndef DEBUG
    if (vis == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(vis, 0, sizeof(bool) * (node_N));
    Q.push(s);
    double minDist = 1e9;
    int minIndex = -1;
    while (!Q.empty()) {
        int cur = Q.front();
        Q.pop();
        vis[cur] = 1;
        for (auto i : mst[cur]) {
            if (!vis[i.e]) {
                if (label[i.e] == label[cur]) {
                    Q.push(i.e);
                } else {
                    if (i.w < minDist) {
                        minDist = i.w;
                        minIndex = i.e;
                    }
                }
            }
        }
    }
    res.c = label[minIndex];
    res.d = minDist;
    delete[] vis;
    return res;
}

// 合并类别标签
void merge_cluster(int c, int c_) {
    for (int i = 0; i < node_N; i++) {
        if (label[i] == c_) label[i] = c;
    }
    cluster_count[c] += cluster_count[c_];
    cluster_count[c_] = 0;
}

// 大于k类时，合并不同类别之间的最小边，基于重心的测定距离以及桥的长度
void merge_with_edges_until_k_cluster(int k) {
    if (cluster_count == NULL) cluster_count = new int[cluster + 1];
#ifndef DEBUG
    if (cluster_count == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    if (any_in_cluster == NULL) any_in_cluster = new int[node_N + 1];
#ifndef DEBUG
    if (any_in_cluster == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    std::priority_queue<Edge, std::vector<Edge>, Edge_greater> Q;
    std::vector<int> barycentres;
    memset(cluster_count, 0, sizeof(int) * (cluster + 1));
    memset(any_in_cluster, -1, sizeof(int) * (cluster + 1));
    for (int i = 0; i < node_N; i++) {
        cluster_count[label[i]]++;
        any_in_cluster[label[i]] = i;
    }
    for (int i = 1; i < cluster; ++i) {
        int g = find_barycentre(i);
        barycentres.push_back(g);
    }
    for (int i = 0; i < barycentres.size(); ++i) {
        int count = 0;
        std::priority_queue<Edge, std::vector<Edge>, Edge_greater> cache;
        for (int j = i + 1; j < barycentres.size(); ++j) {
            Edge edge;
            bool find = false;
            double dist = 0.0;
            double bridge = 0.0;
            int s = barycentres[i];
            int e = barycentres[j];
            shortest_path(s, e, -1, label[s], label[e], 0, find, dist, bridge);
            double avg_max =
                std::min(get_max_edge(label[s]), get_max_edge(label[e]));
            if (find && bridge <= avg_max) {
                edge.s = s;
                edge.e = e;
                edge.w = dist;
                Q.push(edge);
                count++;
            }
        }
        if (count == 0) {
            if (!cache.empty()) {
                Q.push(cache.top());
            }
        }
    }
    for (int cc = cluster; cc > k + 1 && !Q.empty(); --cc) {
        Edge cur = Q.top();
        Q.pop();
        while ((label[cur.s] == label[cur.e] ||
                cluster_count[label[cur.s]] == 0 ||
                cluster_count[label[cur.e]] == 0) &&
               !Q.empty()) {
            cur = Q.top();
            Q.pop();
        }
        int c = label[cur.s];
        int g = find_barycentre(c);
        int c2 = label[cur.e];
        int g2 = find_barycentre(c2);
        int p;
        for (p = 0; p < barycentres.size(); p++) {
            if (barycentres[p] == g) break;
        }
        barycentres.erase(barycentres.begin() + p);
        for (p = 0; p < barycentres.size(); p++) {
            if (barycentres[p] == g2) break;
        }
        barycentres.erase(barycentres.begin() + p);
        merge_cluster(c, c2);
        while (!Q.empty()) {
            Q.pop();
        }
        g = find_barycentre(c);
        barycentres.push_back(g);
        for (int i = 0; i < barycentres.size(); ++i) {
            int count = 0;
            std::priority_queue<Edge, std::vector<Edge>, Edge_greater> cache;
            for (int j = i + 1; j < barycentres.size(); ++j) {
                Edge edge;
                bool find = false;
                double dist = 0.0;
                double bridge = 0.0;
                int s = barycentres[i];
                int e = barycentres[j];
                shortest_path(s, e, -1, label[s], label[e], 0, find, dist,
                              bridge);
                double avg_max =
                    std::min(get_max_edge(label[s]), get_max_edge(label[e]));
                if (find && bridge <= avg_max) {
                    edge.s = s;
                    edge.e = e;
                    edge.w = dist;
                    Q.push(edge);
                    count++;
                }
            }
            if (count == 0) {
                if (!cache.empty()) {
                    Q.push(cache.top());
                }
            }
        }
    }
}

// 大于k类合并，基于最短桥
void merge_with_edges_until_k_cluster_by_closest_d(int k) {
    std::priority_queue<Bridge, std::vector<Bridge>, Bridge_greater> Q;
    memset(cluster_count, 0, sizeof(int) * (cluster + 1));
    memset(any_in_cluster, -1, sizeof(int) * (cluster + 1));
    for (int i = 0; i < node_N; i++) {
        cluster_count[label[i]]++;
        any_in_cluster[label[i]] = i;
    }
    for (int i = 0; i < node_N; i++) {
        for (auto e : mst[i]) {
            int c1 = label[e.s];
            int c2 = label[e.e];
            if (c1 != c2) {
                Q.push(Bridge(e, c1, c2, cluster_count[c1], cluster_count[c2]));
            }
        }
    }
    for (int cc = cluster; cc > k + 1; cc--) {
        Bridge cur_b = Q.top();
        Edge cur = cur_b.edge;
        Q.pop();
        while (label[cur.s] == label[cur.e] && !Q.empty()) {
            cur = Q.top().edge;
            Q.pop();
        }
        // std::cout << label[cur.s] << " " << label[cur.e] << " " << cur.w
        //           << std::endl;
        merge_cluster(label[cur.s], label[cur.e]);
        while (!Q.empty()) Q.pop();
        for (int i = 0; i < node_N; i++) {
            for (auto e : mst[i]) {
                int c1 = label[e.s];
                int c2 = label[e.e];
                if (c1 != c2) {
                    Q.push(Bridge(e, c1, c2, cluster_count[c1],
                                  cluster_count[c2]));
                }
            }
        }
    }
    while (!Q.empty()) Q.pop();
}

// 类别计数
typedef struct Count {
    int index;
    int value;
} Count;

struct Count_greater {
    bool operator()(const Count &c1, const Count &c2) {
        if (c1.value < c2.value) {
            return false;
        } else if (c1.value == c2.value) {
            return c1.index > c2.index;
        } else {
            return true;
        }
    }
};

struct Count_lesser {
    bool operator()(const Count &c1, const Count &c2) {
        if (c1.value < c2.value) {
            return true;
        } else if (c1.value == c2.value) {
            return c1.index < c2.index;
        } else {
            return false;
        }
    }
};

// 按照最小的数量进行合并，从节点数最小的类别开始扩展，基于重心的测定距离和桥的长度
void merge_cluster_below_minNum_minDist(int minNum) {
    cluster_count = new int[cluster + 1];
#ifndef DEBUG
    if (cluster_count == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    any_in_cluster = new int[node_N + 1];
#ifndef DEBUG
    if (any_in_cluster == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }
#endif
    memset(cluster_count, 0, sizeof(int) * (cluster + 1));
    memset(any_in_cluster, -1, sizeof(int) * (cluster + 1));
    for (int i = 0; i < node_N; i++) {
        cluster_count[label[i]]++;
        any_in_cluster[label[i]] = i;
    }
    std::priority_queue<Count, std::vector<Count>, Count_greater> Q;
    for (int i = 1; i < cluster; i++) {
        if (cluster_count[i] <= minNum) {
            Count c;
            c.index = i;
            c.value = cluster_count[i];
            Q.push(c);
        }
    }
    Closest_cluster cc;
    int pppp = cluster;
    while (!Q.empty()) {
        int cur = Q.top().index;
        int cnt = Q.top().value;
        Q.pop();
        // if (cluster_count[cur] != cnt) continue;
        while ((cluster_count[cur] > minNum || cluster_count[cur] == 0) &&
               !Q.empty()) {
            cur = Q.top().index;
            cnt = Q.top().value;
            Q.pop();
        }
        if (cluster_count[cur] > minNum) break;
        
        cc = find_closest_cluster_in_tree(cur);
        if (cc.c != -1) {
            int c_ = cc.c;
            if (cluster_count[cur] > cluster_count[c_]) {
                _swap(cur, c_);
            }
            merge_cluster(c_, cur);
            Count c;
            c.index = c_;
            c.value = cluster_count[c_];
            Q.push(c);
        }
    }
    // for (int i = 0; i < node_N; ++i) {
    //     if (cluster_count[i]) 
    //     std::cout << cluster_count[i] << std::endl;
    // }
}

// 重新组织标签
void rebuild_label() {
    std::map<int, int> m;
    int c = 1;
    for (int i = 0; i < node_N; i++) {
        if (m.find(label[i]) != m.end()) {
            label[i] = m[label[i]];
        } else {
            m.insert(std::make_pair(label[i], c));
            label[i] = c;
            c++;
        }
    }
    cluster = c;
}

// 析构函数
static inline void destory_graph() {
    if (mst != NULL) delete[] mst;
    if (degree != NULL) delete[] degree;
    if (degree_sorted != NULL) delete[] degree_sorted;
    if (label != NULL) delete[] label;
    if (cluster_count != NULL) delete[] cluster_count;
    if (any_in_cluster != NULL) delete[] any_in_cluster;
}

// 主函数
void solve(double *edge, double *dist, int k, int minNum, int knn, double T) {
    create_mst(edge, dist);
    preliminary_clustering(knn, T);
    merge_cluster_below_minNum_minDist(minNum);
    rebuild_label();
    merge_with_edges_until_k_cluster(k);
    rebuild_label();
    merge_with_edges_until_k_cluster_by_closest_d(k);
    rebuild_label();
}

#ifndef DEBUG
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *output, *edge, *dist;
    double *p1, *p2, *p3, *p4;

    if (nrhs < 3 || nrhs > 7)
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: wrong number of input arguments.");
    if (nlhs > 1)
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: wrong number of output arguments.");

    edge_N = mxGetM(prhs[0]);
    node_N = edge_N + 1;
    int N_ = mxGetM(prhs[1]);
    if (edge_N != N_)
        mexErrMsgTxt(
            "BFS_AND_CLUSTER.MEX: edges'row and dists'row must agree.");

    edge = mxGetPr(prhs[0]);
    dist = mxGetPr(prhs[1]);
    complete_graph_dist = mxGetPr(prhs[2]);

    int k;
    if (nrhs > 3) {
        p1 = mxGetPr(prhs[3]);
        k = *p1;
    } else {
        mexWarnMsgTxt("Using default k!");
        k = 2;
    }

    if (nrhs > 4) {
        p2 = mxGetPr(prhs[4]);
        minNum = *p2;
    } else {
        mexWarnMsgTxt("Using default minNum!");
        minNum = 0.25 * edge_N / k;
    }

    if (nrhs > 5) {
        p3 = mxGetPr(prhs[5]);
        knn = *p3;
    } else {
        mexWarnMsgTxt("Using default knn!");
        knn = 3;
    }

    if (nrhs > 6) {
        p4 = mxGetPr(prhs[6]);
        T = *p4;
        if (T < 0)
            T = 0.0;
        else if (T > 1)
            T = 1.0;
    } else {
        mexWarnMsgTxt("Using default knn!");
        T = 0.0;
    }

    plhs[0] = mxCreateDoubleMatrix(node_N, 1, mxREAL);
    output = mxGetPr(plhs[0]);
    if (output == NULL) {
        mexErrMsgTxt("BFS_AND_CLUSTER.MEX: Do not have enough memory.");
    }

    solve(edge, dist, k, minNum, knn, T);

    for (int i = 0; i < node_N; i++) {
        output[i] = label[i];
    }

    destory_graph();
}
#endif

#ifdef DEBUG
int main() {
    edge_N = 5;
    node_N = edge_N + 1;
    alpha = 1.5;
    int k = 2;
    int mn = 1;  // 0.25 * node_N / k;
    double in[10] = {1, 2, 3, 4, 5, 2, 3, 4, 5, 6};
    double dist[5] = {1, 1, 1.414, 1, 1};
    create_mst((double *)in, dist);
    preliminary_clustering();
    merge_cluster_below_minNum_minDist(mn, 1e9);
    rebuild_label();
    merge_with_edges_until_k_cluster(k);
    rebuild_label();

    for (int i = 0; i < node_N; i++) {
        std::cout << i << ": " << label[i] << std::endl;
    }

    bool find = false;

    destory_graph();
    return 0;
}
#endif