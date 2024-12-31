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
        return connection.size();                                   // возвращает количество связей атома
    }
    

    int& operator[](int i){                                        // возвращает номер i-того связаного атома
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
    int x1 = 1000000000, y1 = 1000000000;
    int x2 = -1000000000, y2 = -1000000000;
    double alpha = pi / 6;



    Molecule(){                                                 // конструктор метана
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

    
    Molecule(Molecule& other, int number){                                  // конструктор присоединяющий к n-ному атому еще один
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

    void centroid_dfs(int v){                                                          // dfs для поиска центроида
        sz[v] = 1;
        for (int i = 0; i < graph[v].size(); ++i){
            if (sz[graph[v][i]] != 0){
                continue;
            }
            centroid_dfs(graph[v][i]);
            sz[v] += sz[graph[v][i]];
        }
    }


    int getCentroid(){                                                                  // ищет первый центроид
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


    void getCentroids(){                                                                   // определяет оба центроида
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


    bool operator==(Molecule other){                                            // сравнивает AHU молекул
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



    void diameter_dfs(int v, int d, int p, pair<int, int>& mx){                            // запускает bfs для определения атомов, принадлежащих к диаметру                      
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



    vector<int> coordinates_bfs(int s,int f) {                                               // возвращает массив с номерами вершин диаметра
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

    double get_engle(Atom& f, Atom& t){                              // определяет угол между горизонталью и линией, соединяющей координаты двух атомов
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



    void coordinate_dfs(int v, int p, int st, vector<int>& diameter){            // dfs определяющий координаты каждого атома
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
        for(int i = 0; i < this->size(); i++){
            if(graph[i].x < this->x1){
                this->x1 = graph[i].x;
            }
            if(graph[i].y < this->y1){
                this->y1 = graph[i].y;
            }
            if(graph[i].x > this->x2){
                this->x2 = graph[i].x;
            }
            if(graph[i].y > this->y2){
                this->y2 = graph[i].y;
            }
        }
    }


    int size(){         
        return graph.size();
    }


    Atom& operator[](int i){
        return graph[i];
    }




    void draw(int dx = 0, int dy = 0){                                                  // выводит координаты всех атомов
        for(int i = 0; i < this->size(); ++i){                                                 
            cout << i << ": (" << graph[i].x << ":" << graph[i].y << ")" << endl;                                                               
        }
    }

    void print(){                                                                   // выводит список смежности
        for(int i = 0; i < size(); ++i){
            cout << graph[i].number << ": ";
            for(int j  = 0; j < graph[i].size(); ++j){
                cout << graph[i][j] << ' ';
            }
            cout << endl;
        }
        cout << endl;
    }

    void print(vector<int> v){                                                     // выводит вектор
        cout << "size = " << v.size() << endl;
        for(int i = 0; i < v.size(); ++i){
            cout << v[i] << ' ';
        }
        cout << endl;
    }



    void print(ofstream& file_out){                                               // вводит список смежности в файл       
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


    void draw_dfs(int v, int p, vector<int>& used, Molecule& mol){                     // dfs для записи в s координат линий в формате для svg
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


    void read_molecule(Molecule& mol){                                          // заполняет файл координатами линий в формате svg
        vector<int> used(mol.size());
        s = "";
        s += "<svg xmlns=\"http://www.w3.org/2000/svg\">";
        s += '\n';
        draw_dfs(0, -1, used, mol);
        s += "</svg>";
        s += '\n';
    }
};



class Page{
public:
    string s = "";
    int y_max = 0;

    void draw_dfs(int v, int p, Molecule& mol, int rx, int ry){                     // dfs для записи в s координат линий в формате для svg
        for (int i = 0; i < mol[v].size(); ++i){
            if(mol[v][i] == p){
                continue;
            }
            s += "        <line x1=\"";
            s += to_string(mol[v].x + rx);
            s += "\" y1=\"";
            s += to_string(mol[v].y + ry);
            s += "\" x2=\"";
            s += to_string(mol[mol[v][i]].x + rx);
            s += "\" y2=\"";
            s += to_string(mol[mol[v][i]].y + ry);
            s += "\" style=\"stroke:rgb(0,0,0); stroke-width:1\" />";
            s += '\n';
            draw_dfs(mol[v][i], v, mol, rx, ry);
        }
    }


    void read_molecules(vector<Molecule>& molecules){                                          // заполняет файл координатами линий в формате svg
        s = "";
        s += "<!DOCTYPE html>\n";
        s += "<html lang=\"en\">\n";
        s += "<head>\n";
        s += "    <meta charset=\"UTF-8\">\n";
        s += "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        s += "    <title>Document</title>\n";
        s += "</head>\n";
        s += "<body>\n";
        int rx = 0, ry = 0;
        int dy = 0;
        s += "    <svg height=\"\" width=\"1500\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        s += '\n';
        for(int i = 0; i < molecules.size(); ++i){
            if(molecules[i].x2 + rx > 1500){
                rx = 0;
                ry += dy + 25;
                dy = 0;
            }
            draw_dfs(0, -1, molecules[i], rx, ry);
            dy = max(dy, molecules[i].y2);
            rx += molecules[i].x2;
            rx += 25;
            s += '\n';
            y_max = max(y_max, ry + dy);
        }
        s += "    </svg>\n";
        s += "</body>\n";
        s += "</html>\n";
        s += '\n';
    }
};


class Unite_page{
public:
    string s = "";

    Unite_page(int n){                                          
        s = "";
        s += "<!DOCTYPE html>\n";
        s += "<html lang=\"en\">\n";
        s += "<head>\n";
        s += "    <meta charset=\"UTF-8\">\n";
        s += "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        s += "    <title>Document</title>\n";
        s += "    <style>\n";
        for(int i = 0; i < n; ++i){
            s += "      #C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += " {\n";
            s += "          position: fixed;\n";
            s += "          top: ";
            s += to_string(335 - (n / 2 + n % 2) * 30 + i * 30);
            s += "px;\n";
            s += "          left: 600px;\n";
            s += "          height: 30px;\n";
            s += "          width: 300px;\n";
            s += "          font-size: 1em;\n";
            s += "      }\n";
        }
        s += "    </style>\n";
        s += "</head>\n";
        s += "<body>\n";
        
        for(int i = 0; i < n;  ++i){
            s += "    <form id=\"C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += "\" action=\"./pages/C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += ".html\">\n";
            s += "        <input id=\"C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += "\" type=\"submit\" value=\"C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += "\" />\n";
            s += "    </form>\n";
        }

        s += "</body>\n";
        s += "</html>\n";
        s += '\n';
    }
};


int main(){
    clock_t tStart = clock();
    cout << "start time = " << tStart << endl;
    int n = 16;

    Unite_page unite_page(n);
    if(true){
        string s = ".\\index.html";
        ofstream file;
        file.open(s.c_str()); // <- here
        file << unite_page.s;
        file.close();
    }

    string folder_name = "pages";
    CreateDirectoryA(folder_name.c_str(), NULL);



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


    vector<Page> pages;
    for(int i = 0; i < molecules.size(); ++i){
        Page pg;
        pg.read_molecules(molecules[i]);
        pg.s.insert(202, to_string(pg.y_max));
        pages.push_back(pg);
    }

    for(int i = 0; i < pages.size(); ++i){
        string s = ".\\pages\\C";
        s += char((i + 1) / 10 + '0');
        s += char((i + 1) % 10 + '0');
        s += ".html";
        ofstream file;
        file.open(s.c_str()); // <- here
        file << pages[i].s;
        file.close();
    }

    cout << "molecules.size() = " << molecules.size() << endl;
    cout << "pages.size() = " << pages.size() << endl;

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