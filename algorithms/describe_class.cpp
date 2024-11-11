#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <cmath>


using namespace std;


void print(vector<int> v){
    cout << "size = " << v.size() << endl;
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << ' ';
    }
    cout << endl;
}


void mistake(){
    throw length_error("carbon forms only 4 bonds");
}


class Atom{
public:
    bool diametre;
    int number = 0;
    int x = 0, y = 0;
    vector<int> connection;


    Atom(int _number = 0){
        number = _number;
    }


    void connect(Atom& other){                                    // добавляет атому новую связь
        this->connection.push_back(other.number);
        other.connection.push_back(this->number);
        if (this->size() > 4 || other.size() > 4){
            mistake();
        }
    }

    int size(){
        return connection.size();
    }
    

    int& operator[](int i){
        return connection[i];
    }

};


class Molecule{
public:
    vector<Atom> graph;
    vector<int> root;
    vector<string> root_AHU;
    vector<int> sz;              
    int L = 10;


    Molecule(){
        Atom C;
        graph.push_back(C);
        sz.resize(this->size());
        getCentroids();
        if(root.size() == 1){
            vector<string> AHU(this->size());
            build_AHU(AHU, root[0], -1);
            root_AHU.push_back(AHU[root[0]]);
        } else if(root.size() == 2){
            vector<string> AHU1(this->size());
            vector<string> AHU2(this->size());
            build_AHU(AHU1, root[0], -1);
            build_AHU(AHU2, root[1], -1);
            root_AHU.push_back(AHU1[root[0]]);
            root_AHU.push_back(AHU2[root[1]]);
        }
        give_coordinates();
    }

    
    Molecule(Molecule& other, int number){
        Atom C(other.size());
        this->graph = other.graph;
        this->graph.push_back(C);
        this->graph[number].connect(this->graph[this->size() - 1]);
        this->sz.resize(this->size());
        this->getCentroids();
        if(root.size() == 1){
            vector<string> AHU(this->size());
            build_AHU(AHU, root[0], -1);
            root_AHU.push_back(AHU[root[0]]);
        } else if(root.size() == 2){
            vector<string> AHU1(this->size());
            vector<string> AHU2(this->size());
            build_AHU(AHU1, root[0], -1);
            build_AHU(AHU2, root[1], -1);
            root_AHU.push_back(AHU1[root[0]]);
            root_AHU.push_back(AHU2[root[1]]);
        }
        give_coordinates();
    }

    void centroid_dfs(int v){                                                         
        sz[v] = 1;
        for (int i = 0; i < graph[v].size(); ++i){
            if (sz[graph[v][i]] != 0){
                continue;
            }
            centroid_dfs(graph[v][i]);
            sz[v] += sz[graph[v][i]];
        }
    }


    int getCentroid(){                                                                          // ищет первый центроид
        int v = 0;
        centroid_dfs(v);
        while(true){
            int w = -1;
            for(int i = 0; i < graph[v].size(); ++i){
                if(sz[graph[v][i]] > sz[v]){
                    continue;
                }
                if (2 * sz[graph[v][i]] > graph.size()){
                    w = graph[v][i];
                    break;
                }
            }
            if (w == -1){
                break;
            }
            v = w;
        }
        return v;
    }


    void getCentroids(){                                                                          // определяет оба центроида
        int v = getCentroid();
        root.push_back(v);
        for (int i = 0; i < graph[v].size(); ++i){
            if (2 * sz[graph[v][i]] == graph.size()){         
                root.push_back(graph[v][i]);
            }
        }
    }


    void build_AHU(vector<string>& AHU, int v, int p){                                          // определяет AHU для всего графа
        vector<string> children;
        for (int i = 0; i < graph[v].size(); ++i){
            if (graph[v][i] == p){
                continue;
            }
            build_AHU(AHU, graph[v][i], v);
            children.push_back(AHU[graph[v][i]]);
        }
        sort(children.begin(), children.end());
        string sm = "";
        for (int i = 0; i < children.size(); ++i){
            sm += children[i];
        }
        AHU[v] = "(0" + sm + ")";
    }


    bool operator==(Molecule other){
        if(this->size() != other.size()){
            return false;
        }
        if(this->root_AHU.size() == 1 && other.root_AHU.size() == 1){
            return this->root_AHU[0] == other.root_AHU[0];
        } else if(this->root_AHU.size() == 2 && other.root_AHU.size() == 2){
            return this->root_AHU[0] == other.root_AHU[0] || this->root_AHU[1] == other.root_AHU[1] || this->root_AHU[0] == other.root_AHU[1] || this->root_AHU[1] == other.root_AHU[0];
        } else {
            return false;
        }
    }   





    void diameter_dfs(int v, int d, int p, pair<int, int>& mx){                                                     
        if(mx.first < d){
            mx.first = d;
            mx.second = v;
        }
        for(int i = 0; i < graph[v].size(); ++i){
            if(graph[v][i] != p){
                diameter_dfs(graph[v][i], d + 1, v, mx);
            }
        }
    }



    vector<int> coordinates_bfs(int s,int f) {                                                          // возвращает массив с номерами вершин диаметра
        vector<int> parent(this->size());
        vector<int> dist(this->size(), 1000000000);
        vector<int> rez;
        queue<int> q;
        parent[s] = -1;
        dist[s] = 0;
        q.push(s);
        while(q.size()){
            int v = q.front();
            q.pop();
            for(int i = 0; i < graph[v].size(); ++i){
                if(dist[graph[v][i]] > dist[v] + 1){
                    dist[graph[v][i]] = dist[v] + 1;
                    q.push(graph[v][i]);
                    parent[graph[v][i]] = v;
                }
            }
        }
        
        int cur = f;
        while(cur != -1){
            rez.push_back(cur);
            cur = parent[cur];
        }
        return rez;

    }


    void coordinate_dfs(int v, int p, int st, vector<int>& diameter){            // тут доделай
        if(p == -1){
            graph[v].x = 0;
            graph[v].y = 0;
        } else {

        }
        for(int i = 0; i < graph[v].size(); ++i){
            if(graph[v][i] != p){
                coordinate_dfs(graph[v][i], v, !st, diameter);
            }
        }


    }


    void give_coordinates(){                                                            // задает координаты атомам сразу с созданием молекулы
        vector<int> diameter;                                                           
        pair<int, int> edge1 = {0, 0};
        pair<int, int> edge2 = {0, 0};                                                  // доделай
        diameter_dfs(this->root[0], 0, -1, edge1);
        diameter_dfs(edge1.second, 0, -1, edge2);
        diameter = coordinates_bfs(edge1.second, edge2.second);
        coordinate_dfs(diameter[0], -1, 0, diameter);
    }


    int size(){         
        return graph.size();
    }


    Atom& operator[](int i){
        return graph[i];
    }



    void draw(int dx = 0, int dy = 0){                                                  // рисует все связи
        for(int i = 0; i < this->size(); ++i){                                   //доделай            
                                                                                    
        }
    }


    void print(){
        for(int i = 0; i < size(); ++i){
            cout << graph[i].number << ": ";
            for(int j  = 0; j < graph[i].size(); ++j){
                cout << graph[i][j] << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }

    void print(vector<int> v){
        cout << "size = " << v.size() << endl;
        for(int i = 0; i < v.size(); ++i){
            cout << v[i] << ' ';
        }
        cout << endl;
    }
};







int main(){
    Molecule C1;
    Molecule C2(C1, 0);
    Molecule C3(C2, 0);
    Molecule C40(C3, 0);
    Molecule C41(C3, 1);
    Molecule C42(C3, 2);
    C40.print();
    C41.print();
    C42.print();




}