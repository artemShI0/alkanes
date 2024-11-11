#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "raylib.h"

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
    int number = 0;
    int x = 0, y = 0;
    vector<int> connection;


    Atom(int _number = 0){
        number = _number;
    }


    void connect(Atom& other){
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
    vector<int> sz;               // sz.size() = graph.size();
    
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


    int getCentroid(){
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


    void getCentroids(){
        int v = getCentroid();
        root.push_back(v);
        for (int i = 0; i < graph[v].size(); ++i){
            if (2 * sz[graph[v][i]] == graph.size()){         
                root.push_back(graph[v][i]);
            }
        }
    }


    void build_AHU(vector<string>& AHU, int v, int p){
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





    void coordinates_dfs(int v, int d, int p, pair<int, int>& mx){
        if(mx.first < d){
            mx.first = d;
            mx.second = v;
        }
        for(int i = 0; i < graph[v].size(); ++i){
            if(graph[v][i] != p){
                coordinates_dfs(graph[v][i], d + 1, v, mx);
            }
        }
    }

    void give_coordinates(){
        vector<int> diametre;
        pair<int, int> edge1 = {0, 0};
        pair<int, int> edge2 = {0, 0};
        coordinates_dfs(this->root[0], 0, -1, edge1);
        coordinates_dfs(edge1.second, 0, -1, edge2);
        cout << edge1.first << ' ' << edge1.second << endl;
        cout << edge2.first << ' ' << edge2.second << endl;
    }




    int size(){
        return graph.size();
    }


    Atom& operator[](int i){
        return graph[i];
    }



    void draw(int dx = 0; int dy = 0){
        for(int i = 0; i < n; ++i){
            
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
};







int main(){
    Molecule C1;
    Molecule C2(C1, 0);
    Molecule C3(C2, 0);
    cout << "djkaflsjkdfkas" << endl << endl << endl;
    Molecule C40(C3, 0);
    Molecule C41(C3, 1);
    Molecule C42(C3, 2);
    C40.print();
    C41.print();
    C42.print();

    cout << (C40 == C41) << endl;
    cout << (C41 == C42) << endl;
    cout << (C42 == C40) << endl;


    InitWindow(800, 600, "Hello, World!");
    
    while (!WindowShouldClose())
    {
        BeginDrawing();
        
        ClearBackground(RAYWHITE);
        C4.draw();
        
        EndDrawing();
    }

    CloseWindow();


}
