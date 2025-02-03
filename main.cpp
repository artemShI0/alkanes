#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <cmath>
#include <fstream>
#include <windows.h>
#include <time.h>
using namespace std;


const double pi = 3.1415926535;


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
    int L = 100;
    double alpha = pi / 6;



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



    bool include(vector<int>& diameter, int v){                   // говорит, есть ли элемент в массиве
        for(int i = 0; i < diameter.size(); ++i){
            if(v == diameter[i]){
                return true;
            }
        }
        return false;
    }

    int next_diameter(vector<int>& diameter, int v){                   // индекс следующего атома, принадлежащего диаметру
        for(int i = 0; i < diameter.size(); ++i){
            if(v == diameter[i]){
                return i;
            }
        }
        return -1;
    }

    double get_engle(Atom& f, Atom& t){
        double dx = t.x - f.x;
        double dy = t.y - f.y;
        if(dx == 0 && t.y >= f.y){
            return pi / 2;
        }
        if(dx == 0 && t.y < f.y){
            return - pi / 2;
        }
        if(dy == 0 && t.x >= f.y){
            return 0;
        }
        if(dy == 0 && t.x < f.y){
            return pi;
        }
        double engle = atan(dy / dx);
        if(dx < 0){
            engle += pi;
        }
        return engle;
    }



    void coordinate_dfs(int v, int p, int st, vector<int>& diameter){            // тут доделай
        this->L -= 1;
        if(p == -1 && graph[v].size()){
            graph[v].x = 0;
            graph[v].y = 0;
            graph[graph[v][0]].x = graph[v].x + L * cos(alpha + 2 * alpha * (st == 0 ? -1 : 1));
            graph[graph[v][0]].y = graph[v].y + L * sin(alpha + 2 * alpha * (st == 0 ? -1 : 1));
            coordinate_dfs(graph[v][0], v, !st, diameter);
            return;
        }
        int used = 0;
        for(int i = 0; i < graph[v].size(); ++i){
            if(graph[v][i] != p){
                double fi = get_engle(graph[p], graph[v]);
                double _fi = 2 * pi - fi;
                if(used == 0){
                    if(st == 0){
                        graph[graph[v][i]].x = graph[v].x + L * cos(- 2 * alpha + fi);
                        graph[graph[v][i]].y = graph[v].y + L * sin(- 2 * alpha + fi);
                        coordinate_dfs(graph[v][i], v, 1, diameter);
                    } else {
                        graph[graph[v][i]].x = graph[v].x + L * cos(2 * alpha - _fi);
                        graph[graph[v][i]].y = graph[v].y + L * sin(2 * alpha - _fi);
                        coordinate_dfs(graph[v][i], v, 0, diameter);
                    }
                } else if(used == 1){
                    if(st == 0){
                        graph[graph[v][i]].x = graph[v].x + L * cos(pi / 2 + fi - 2 * alpha);
                        graph[graph[v][i]].y = graph[v].y + L * sin(pi / 2 + fi - 2 * alpha);
                        coordinate_dfs(graph[v][i], v, 1, diameter);
                    } else {
                        graph[graph[v][i]].x = graph[v].x + L * cos(- pi / 2 - _fi + 2 * alpha);
                        graph[graph[v][i]].y = graph[v].y + L * sin(- pi / 2 - _fi + 2 * alpha);
                        coordinate_dfs(graph[v][i], v, 0, diameter);
                    }
                } else if(used == 2){
                    if(st == 0){
                        graph[graph[v][i]].x = graph[v].x + L * cos(- pi / 2 + fi);
                        graph[graph[v][i]].y = graph[v].y + L * sin(- pi / 2 + fi);
                        coordinate_dfs(graph[v][i], v, 0, diameter);
                    } else {
                        graph[graph[v][i]].x = graph[v].x + L * cos(pi / 2 - _fi);
                        graph[graph[v][i]].y = graph[v].y + L * sin(pi / 2 - _fi);
                        coordinate_dfs(graph[v][i], v, 1, diameter);
                        // graph[graph[v][i]].x = graph[v].x + L * cos(fi + pi / 2);
                        // graph[graph[v][i]].y = graph[v].y + L * sin(fi + pi / 2);
                    }
                }
                used++;
            }
        }
    }


    void give_coordinates(){                                                            // задает координаты атомам сразу с созданием молекулы
        vector<int> diameter;                                                           
        pair<int, int> edge1 = {0, 0};
        pair<int, int> edge2 = {0, 0};                                                  // доделай
        diameter_dfs(this->root[0], 0, -1, edge1);
        diameter_dfs(edge1.second, 0, -1, edge2);
        diameter = coordinates_bfs(edge1.second, edge2.second);                               // в diameter хранятся индексы вершин, принадлежащих одному из диаметров
        coordinate_dfs(diameter[0], -1, 0, diameter);
        int x0 = 0, y0 = 0;
        for(int i = 0; i < this->size(); i++){
            if(graph[i].x < x0){
                x0 = graph[i].x;
            }
            if(graph[i].y < y0){
                y0 = graph[i].y;
            }
        }
        for(int i = 0; i < this->size(); i++){
            graph[i].x -= x0;
            graph[i].y -= y0;
        }
    }


    int size(){         
        return graph.size();
    }


    Atom& operator[](int i){
        return graph[i];
    }




    void draw(int dx = 0, int dy = 0){                                                  // рисует все связи
        for(int i = 0; i < this->size(); ++i){                                   //доделай       
            cout << i << ": (" << graph[i].x << ":" << graph[i].y << ")" << endl;                                                               
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



    void print(ofstream& file_out){
        for(int i = 0; i < size(); ++i){
            file_out << graph[i].number << ": ";
            for(int j  = 0; j < graph[i].size(); ++j){
                file_out << graph[i][j] << ' ';
            }
            file_out << '\n';
        }
        file_out << '\n';
    }

};


class SVG_picture{
public:
    string s = "";


    void draw_dfs(int v, int p, vector<int>& used, Molecule& mol){
        for (int i = 0; i < mol[v].size(); ++i){
            if(mol[v][i] == p){
                continue;
            }
            s += "<line x1=\"";
            s += to_string(mol[v].x);
            s += "\" y1=\"";
            s += to_string(mol[v].y);
            s += "\" x2=\"";
            s += to_string(mol[mol[v][i]].x);
            s += "\" y2=\"";
            s += to_string(mol[mol[v][i]].y);
            s += "\" style=\"stroke:rgb(0,0,0); stroke-width:1\" />";
            s += '\n';
            draw_dfs(mol[v][i], v, used, mol);
        }
    }


    void read_molecule(Molecule& mol){
        vector<int> used(mol.size());
        s = "";
        s += "<svg xmlns=\"http://www.w3.org/2000/svg\">";
        s += '\n';
        draw_dfs(0, -1, used, mol);
        s += "</svg>";
        s += '\n';
    }
};



int main(){
    clock_t tStart = clock();
    cout << "start time = " << tStart << endl;
    int n = 9;



    for(int i = 1; i < n + 1; ++i){
        string s = ".\\skelets\\C";
        s += char(i / 10 + '0');
        s += char(i % 10 + '0');
        CreateDirectoryA(s.c_str(), NULL);
    }



    vector<vector<Molecule>> molecules(n);
    Molecule C1;
    molecules[0].push_back(C1);
    for(int i = 1; i < molecules.size(); ++i){
        for(int j = 0; j < molecules[i - 1].size(); ++j){
            for(int k = 0; k < molecules[i - 1][j].size(); ++k){
                if(molecules[i - 1][j][k].size() > 3){
                    continue;
                }
                Molecule C(molecules[i - 1][j], k);
                bool already = false;
                for(int l = 0; l < molecules[i].size(); ++l){
                    if(molecules[i][l] == C){
                        already = true;
                    }
                }
                if(!already){
                    molecules[i].push_back(C);
                }
            }
        }
        cout << "C" << i + 1 << ": " << molecules[i].size() << " isomers" << endl;
    }


    vector<vector<SVG_picture>> pictures(n);
    for(int i = 0; i < molecules.size(); ++i){
        for(int j = 0; j < molecules[i].size(); ++j){
            SVG_picture pict;
            pict.read_molecule(molecules[i][j]);
            pictures[i].push_back(pict);
        }
    }
    for(int i = 0; i < pictures.size(); ++i){
        for(int j = 0; j < pictures[i].size(); ++j){
            string s = ".\\skelets\\C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += "\\";
            s += to_string(j);
            s += ".svg";
            ofstream file;
            file.open(s.c_str()); // <- here
            file << pictures[i][j].s;
            file.close();
        }
    }


    // SVG_picture picture;
    // picture.read_molecule(molecules[3][0]);
    // molecules[3][0].print();
    // molecules[3][0].draw();
    // ofstream file;
    // string name = "demo2";
    // string file_name = ".\\" + name + ".svg";
    // file.open(file_name.c_str()); // <- here
    // file << picture.s;
    // file.close();



    cout << "end time = " << clock() << endl;
    cout << "all time: " << 1.0 * (clock() - tStart)/CLOCKS_PER_SEC;

    return 0;
}