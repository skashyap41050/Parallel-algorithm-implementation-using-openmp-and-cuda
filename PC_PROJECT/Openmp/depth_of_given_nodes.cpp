#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

class Graph
{
private:
    vector<vector<int>> G;
    vector<vector<int>> edges_number;
    int prefix_sum(vector<int> &a, int low, int hi)
    {
        if (low == hi)
            return a[low];
        int left_prefix_sum = 0;
        int right_prefix_sum = 0;
        int mid = (hi - low) / 2 + low;
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            left_prefix_sum = prefix_sum(a, low, mid);
#pragma omp section
            right_prefix_sum = prefix_sum(a, mid + 1, hi);
        }
#pragma omp parallel for
        for (int i = mid + 1; i <= hi; i++)
        {
            a[i] = a[i] + left_prefix_sum;
        }
        return left_prefix_sum + right_prefix_sum;
    }
    vector<int> listRanking(vector<int> next, int root)
    {
        vector<int> distance(next.size());
        int sink;
#pragma omp parallel for num_threads(next.size())
        for (int i = 0; i < next.size(); i++)
        {
            if (i == next[i])
            {
                sink = i;
                distance[i] = 0;
            }
            else
            {
                distance[i] = 1;
            }
        }
        vector<int> temp_distance(next.size());
        vector<int> temp_next(next.size());

        for (int j = 0; j <= log2(next.size()); j++)
        {
#pragma omp parallel for
            for (int i = 0; i < next.size(); i++)
            {
                if (next[i] != sink)
                {
                    temp_distance[i] = distance[i] + distance[next[i]];
                    temp_next[i] = next[next[i]];
                }
                else
                {
                    temp_distance[i] = distance[i];
                    temp_next[i] = next[i];
                }
            }
            temp_distance.swap(distance);
            temp_next.swap(next);
        }
        cout << "rank of edges" << endl;
        for (int i = 0; i < distance.size(); i++)
            cout << i << " " << distance[i] << endl;
        return distance;
    }

public:
    void addEdge(int u, int v)
    {
        G[u].push_back(v);
        G[v].push_back(u);
    }
    Graph(int n)
    {
        G.resize(n, vector<int>());
        edges_number.resize(n, vector<int>(n));
    }
    void cal_degree(vector<int> &deg)
    {
        for (int i = 0; i < G.size(); i++)
        {
            deg.push_back(G[i].size()); 
        }
    }
    void make_labels(vector<pair<int, int>> &edge_label, vector<int> &deg)
    {
#pragma omp parallel for num_threads(G.size())
        for (int i = 0; i < G.size(); i++)
        {
            if (i == 0)
            {
                int k = 0;
                for (int j = 0; j < deg[i]; j++)
                {
                    edge_label[j] = make_pair(i, G[i][k]);
                    edges_number[i][G[i][k]] = j;
                    k++;
                }
            }
            else
            {
                int k = 0;
                for (int j = deg[i - 1]; j < deg[i]; j++)
                {
                    edge_label[j] = make_pair(i, G[i][k]);
                    edges_number[i][G[i][k]] = j;
                    k++;
                }
            }
        }
    }
    vector<int> find_parent(vector<int> &rank, vector<pair<int, int>> &edge_label, int root)
    {
        vector<int> parent(G.size());
#pragma omp parallel for num_threads(edge_label.size())
        for (int i = 0; i < edge_label.size(); i++)
        {
            int u = edge_label[i].first;
            int v = edge_label[i].second;
            int u_to_v = edges_number[u][v];
            int v_to_u = edges_number[v][u];
            if (rank[u_to_v] < rank[v_to_u])
                parent[u] = v;
            else
                parent[v] = u;
        }
        parent[root] = root;
        cout << "parents are" << endl;
        for (int i = 0; i < parent.size(); i++)
            cout << i << "parent is" << parent[i] << endl;
        cout << endl;
        return parent;
    }
    void find_successor(vector<int> &successor, const vector<pair<int, int>> &edge_label, int root)
    {
    
        for (int i = 0; i < edge_label.size(); i++)
        {
            int u = edge_label[i].first;
            int v = edge_label[i].second;
        
            int t;
            for (int k = 0; k < G[v].size(); k++)
            {
                if (G[v][k] == u)
                {
                
                    t = G[v][(k + 1) % (G[v].size())];
                    break;
                }
            }
        
            successor[i] = edges_number[v][t];
        }
        cout << "Successor" << endl;
        
        int last_node = G[root][G[root].size() - 1];
        int last_node_edge_number = edges_number[last_node][root];
        successor[last_node_edge_number] = last_node_edge_number;
        for (int i = 0; i < successor.size(); i++)
        {
            cout << i << " " << successor[i] << endl;
        }
    }
    void unrootedToRooted(int root)
    {
        vector<int> deg;
        cal_degree(deg);
        prefix_sum(deg, 0, deg.size() - 1);

        vector<pair<int, int>> edge_label(deg[deg.size() - 1]);
        make_labels(edge_label, deg);

        vector<int> successor(edge_label.size());
        find_successor(successor, edge_label, root);

        vector<int> rank = listRanking(successor, root);
        vector<int> parent = find_parent(rank, edge_label, root);

        vector<int> depth = cal_depth(edge_label, rank, parent, root);
    }

    vector<int> cal_depth(vector<pair<int, int>> &edge_label, vector<int> &rank, vector<int> &parent, int root)
    {
        vector<int> edge_weight(edge_label.size());
#pragma omp parallel for
        for (int i = 0; i < edge_label.size(); i++)
        {
            int u, v;
            u = edge_label[i].first;
            v = edge_label[i].second;
            if (parent[v] == u)
            {
                edge_weight[i] = 1;
            }
            else
                edge_weight[i] = -1;
        }
        vector<int> new_edge_weight(edge_label.size());
#pragma omp parallel for
        for (int i = 0; i < edge_label.size(); i++)
        {
            int k = (edge_label.size() - 1) - rank[i];
            new_edge_weight[k] = edge_weight[i];
        }
        prefix_sum(new_edge_weight, 0, new_edge_weight.size() - 1);
        vector<int> depth(G.size());
#pragma omp parallel for
        for (int i = 0; i < edge_label.size(); i++)
        {
            int pos = (edge_label.size() - 1) - rank[i];
            edge_weight[i] = new_edge_weight[pos];
        }

#pragma omp parallel for
        for (int i = 0; i < G.size(); i++)
        {
            int no = edges_number[parent[i]][i];
            depth[i] = edge_weight[no];
        }
        depth[root] = 0;
        cout << "depth of nodes are" << endl;
        for (int i = 0; i < depth.size(); i++)
            cout << i << " depth is " << depth[i] << endl;
        return depth;
    }
};

int main(void)
{
    int N, E;
    cout << "Enter no of vertices" << endl;
    cin >> N;
    cout << "Enter no of edges" << endl;
    cin >> E;
    cout << "Enter edges pair_wise" << endl;
    Graph G(N);
    for (int i = 0; i < E; i++)
    {
        int u, v;
        cin >> u >> v;
        G.addEdge(u, v);
    }
    cout << "Enter root" << endl;
    int root;
    cin >> root;
    G.unrootedToRooted(root);
    return 0;
}
