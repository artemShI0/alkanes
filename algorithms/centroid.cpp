#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

void print(vector<int> v){
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << ' ';
    }
    cout << endl;
}

void print(vector<vector<int>> v){
    for(int i = 0; i < v.size(); ++i){
        for(int j = 0; j < v[i].size(); ++j){
            cout << v[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void dfs(vector<vector<int>>& graph, vector<int>& sz, int v){
	sz[v] = 1;
	for(int i = 0; i < graph[v].size(); ++i){
        if(sz[graph[v][i]] != 0){
            continue;
        }
        dfs(graph, sz, graph[v][i]);
        sz[v] += sz[graph[v][i]];
    }
}

int getCentroid(vector<vector<int>>& graph, vector<int>& sz, int v){
	dfs(graph, sz, v);
	while(true) {
		int w = -1;
		for(int i = 0; i < graph[v].size(); ++i){
			if(sz[graph[v][i]] > sz[v]){
                continue;
            }
            if(2 * sz[graph[v][i]] > graph.size()){
                w = graph[v][i];
                break;
            }
		}
		if (w == -1) {
            break;
        }
		v = w;
	}
	return v;
}

vector<int> getCentroids(vector<vector<int>>& graph, int v){ //v - any vertex of tree
	vector<int> sz(graph.size());
    v = getCentroid(graph, sz, v);
	vector<int> res;
	res.push_back(v);
    for(int i = 0; i < graph[v].size(); ++i){
        if(2 * sz[graph[v][i]] == graph.size()){
            res.push_back(graph[v][i]);
        }
    }
	return res;
}


void build_AHU(vector<vector<int>>& graph, vector<string>& AHU, int v){
	AHU[v] = "0";
    vector<string> children;
	for(int i = 0; i < graph[v].size(); ++i){
        if(AHU[graph[v][i]] != ""){
            continue;
        }     
        build_AHU(graph, AHU, graph[v][i]);
        children.push_back(AHU[graph[v][i]]);
    }
    sort(children.begin(), children.end());
    string sm = "";
    for(int i = 0; i < children.size(); ++i){
        sm +=children[i];
    }
    AHU[v] = "(0" + sm + ")";
}
// children переходит в AHU в узлах и листьях по-разному


bool CheckIsomorphicTrees(vector<vector<int>> graph1, vector<vector<int>> graph2){
    if(graph1.size() != graph2.size()){
        return false;
    }
    int n = graph1.size();
    vector<int> root1 = getCentroids(graph1, 0);  // 1 или 2 центроида
    vector<int> root2 = getCentroids(graph2, 0);
    print(root1);
    print(root2);
    if(root1.size() == 1 && root2.size() == 1){
        vector<string> AHU1(n);
        vector<string> AHU2(n);
        build_AHU(graph1, AHU1, root1[0]);
        build_AHU(graph1, AHU2, root2[0]);
        return AHU1[root1[0]] == AHU2[root2[0]];
    } else if(root1.size() == 2 && root2.size() == 2){
        vector<string> AHU11(n);             
        vector<string> AHU12(n); 
        vector<string> AHU21(n);
        vector<string> AHU22(n); 
        build_AHU(graph1, AHU11, 0);
        build_AHU(graph1, AHU12, 0);
        build_AHU(graph2, AHU21, 0);
        build_AHU(graph2, AHU22, 0);
        return AHU11[root1[0]] == AHU21[root2[0]] || AHU12[root1[1]] == AHU21[root2[0]] || AHU11[root1[0]] == AHU22[root2[1]] || AHU12[root1[1]] == AHU22[root2[1]];
    } else {
        return false;
    }
}





int main(){
    int n, m, k;
    cin >> n;
    vector<vector<int>> graph1(n);
    vector<vector<int>> graph2(n);
    for(int i = 0; i < n; ++i){
        cin >> m;
        for(int j = 0; j < m; ++j){
            cin >> k;
            graph1[i].push_back(k);
        }
    }
    cin >> n;
    for(int i = 0; i < n; ++i){
        cin >> m;
        for(int j = 0; j < m; ++j){
            cin >> k;
            graph2[i].push_back(k);
        }
    }
    cout << CheckIsomorphicTrees(graph1, graph2);

}





/*
3
1 1
2 0 2
1 1
3
1 0
2 1 2
1 0
*/








