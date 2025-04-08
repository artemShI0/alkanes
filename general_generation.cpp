#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <cmath>
#include <fstream>
#include <time.h>
#include <thread>
#include <mutex>
#include <set>
#include <map>
using namespace std;

const double pi = 3.1415926535;
const long long p = 13;
const long long m = 18014398241046527;
map<string, int> idx{
    {"C", 0},
    {"F", 1},
    {"Br", 2},
    {"Cl", 3},
    {"I", 4},
};
        

void print_v(vector<int> v){
    for(int i = 0; i < v.size(); ++i){
        cout << v[i] << ' ';
    }
    cout << endl;
}


void mistake()
{
    throw length_error("carbon forms only 4 bonds");
}



class Atom
{
public:
    bool diametre;
    int number = 0;
    int x = 0, y = 0;
    vector<int> connection;

    Atom(int _number = 0)
    {
        number = _number;
    }

// добавляет атому новую связь
    void connect(Atom &other)
    { 
        this->connection.push_back(other.number);
        other.connection.push_back(this->number);
        if (this->size() > 4 || other.size() > 4)
        {
            mistake();
        }
    }

// возвращает количество связей атома
    int size()
    {
        return connection.size(); 
    }

// возвращает номер i-того связанного атома
    int &operator[](int i)
    { 
        return connection[i];
    }
};



class Molecule
{
public:
    vector<int> chiral_atom;
    vector<int> part_chiral_atom;
    vector<Atom> graph;
    int root1 = -1, root2 = -1;
    string AHU1 = "", AHU2 = "";
    vector<int> sz;
    int L = 100;
    int x1 = 1000000000, y1 = 1000000000;
    int x2 = -1000000000, y2 = -1000000000;
    double alpha = pi / 6;

    // конструктор метана
    Molecule()
    {
        Atom C;
        graph.push_back(C);
        root1 = 0;
        AHU1 = "()";
    }

    // конструктор присоединяющий к n-ному атому еще один
    Molecule(Molecule &other, int number)
    {
        Atom C(other.size());
        graph = other.graph;
        graph.push_back(C);
        graph[number].connect(graph[size() - 1]);
        sz.resize(size());
        getCentroids();
        if (root2 == -1)
        {
            vector<string> ahu(size());
            build_AHU(ahu, root1, -1);
            AHU1 = ahu[root1];
        }
        else
        {
            vector<string> ahu1(size());
            vector<string> ahu2(size());
            build_AHU(ahu1, root1, -1);
            build_AHU(ahu2, root2, -1);
            AHU1 = ahu1[root1];
            AHU2 = ahu2[root2];
        }
    }

    // запускает все методы после того, как определили, что берем молекулу
    void build()
    {
        give_coordinates();
        chiral_atom.resize(size());
        part_chiral_atom.resize(size());
        for(int i = 0; i < size(); ++i){
            string AHU = chiral_AHU_dfs(i, -1);
            vector<string> AHUS;
            set<string> s;
            bool st = true;
            int open = 0; 
            for(int j = 1; j < AHU.size() - 1; ++j){
                open += (AHU[j] == '(' ? 1 : -1);
                if(open >= 1 && !st){
                    AHUS.back().push_back(AHU[j]);
                } else if(open == 1 && st == true){
                    AHUS.push_back("(");
                    st = !st;
                } else if(open == 0) {
                    AHUS.back().push_back(AHU[j]);
                    st = !st;
                }
            }
            for(int j = 0; j < AHUS.size(); ++j){
                s.insert(AHUS[j]);
            }
            if(s.size() > 3 || s.size() == 3 && graph[i].size() == 3){
                chiral_atom[i] = 1;
            }
        }
    }
    string chiral_AHU_dfs(int v, int p)
    {
        vector<string> children;
        for (int i = 0; i < graph[v].size(); ++i)
        {
            if (graph[v][i] == p)
            {
                continue;
            }
            children.push_back(chiral_AHU_dfs(graph[v][i], v));
        }
        sort(children.begin(), children.end());
        string sm = "";
        for (int i = 0; i < children.size(); ++i)
        {
            sm += children[i];
        }
        return "(" + sm + ")";
    }
    // dfs для поиска центроида
    void centroid_dfs(int v)
    {
        sz[v] = 1;
        for (int i = 0; i < graph[v].size(); ++i)
        {
            if (sz[graph[v][i]] != 0)
            {
                continue;
            }
            centroid_dfs(graph[v][i]);
            sz[v] += sz[graph[v][i]];
        }
    }

    // ищет первый центроид
    int getCentroid()
    {
        int v = 0;
        centroid_dfs(v);
        while (true)
        {
            int w = -1;
            for (int i = 0; i < graph[v].size(); ++i)
            {
                if (sz[graph[v][i]] > sz[v])
                {
                    continue;
                }
                if (2 * sz[graph[v][i]] > graph.size())
                {
                    w = graph[v][i];
                    break;
                }
            }
            if (w == -1)
            {
                break;
            }
            v = w;
        }
        return v;
    }

    // определяет оба центроида
    void getCentroids()
    {
        root1 = getCentroid();
        for (int i = 0; i < graph[root1].size(); ++i)
        {
            if (2 * sz[graph[root1][i]] == graph.size())
            {
                root2 = graph[root1][i];
            }
        }
    }

    // определяет AHU для всего графа
    void build_AHU(vector<string> &AHU, int v, int p)
    {
        vector<string> children;
        for (int i = 0; i < graph[v].size(); ++i)
        {
            if (graph[v][i] == p)
            {
                continue;
            }
            build_AHU(AHU, graph[v][i], v);
            children.push_back(AHU[graph[v][i]]);
        }
        sort(children.begin(), children.end());
        string sm = "";
        for (int i = 0; i < children.size(); ++i)
        {
            sm += children[i];
        }
        AHU[v] = "(" + sm + ")";
    }

    // сравнивает AHU молекул
    bool operator==(Molecule other)
    {
        if (AHU1 == other.AHU1)
        {
            return true;
        }
        if (AHU1 == other.AHU2)
        {
            return true;
        }
        return false;
    }

    // запускает dfs для определения атомов одного из края диаметра
    void diameter_dfs(int v, int d, int p, pair<int, int> &mx)
    {
        if (mx.first < d)
        {
            mx.first = d;
            mx.second = v;
        }
        for (int i = 0; i < graph[v].size(); ++i)
        {
            if (graph[v][i] != p)
            {
                diameter_dfs(graph[v][i], d + 1, v, mx);
            }
        }
    }

    // возвращает массив с номерами вершин диаметра
    vector<int> coordinates_bfs(int s, int f)
    {
        vector<int> parent(this->size());
        vector<int> dist(this->size(), 1000000000);
        vector<int> rez;
        queue<int> q;
        parent[s] = -1;
        dist[s] = 0;
        q.push(s);
        while (q.size())
        {
            int v = q.front();
            q.pop();
            for (int i = 0; i < graph[v].size(); ++i)
            {
                if (dist[graph[v][i]] > dist[v] + 1)
                {
                    dist[graph[v][i]] = dist[v] + 1;
                    q.push(graph[v][i]);
                    parent[graph[v][i]] = v;
                }
            }
        }

        int cur = f;
        while (cur != -1)
        {
            rez.push_back(cur);
            cur = parent[cur];
        }
        return rez;
    }

    // определяет угол между горизонталью и линией, соединяющей координаты двух атомов
    double get_engle(Atom &f, Atom &t)
    {
        double dx = t.x - f.x;
        double dy = t.y - f.y;
        if (dx == 0 && t.y >= f.y)
        {
            return pi / 2;
        }
        if (dx == 0 && t.y < f.y)
        {
            return -pi / 2;
        }
        if (dy == 0 && t.x >= f.y)
        {
            return 0;
        }
        if (dy == 0 && t.x < f.y)
        {
            return pi;
        }
        double engle = atan(dy / dx);
        if (dx < 0)
        {
            engle += pi;
        }
        return engle;
    }

    // возвращает индекс вхождения элемента
    int find(Atom &C, int x)
    {
        for (int i = 0; i < C.size(); ++i)
        {
            if (C[i] == x)
            {
                return i;
            }
        }
        return -1;
    }

    // dfs определяющий координаты каждого атома
    void coordinate_dfs(int v, int p, int st, vector<int> &diameter, int dist)
    { 
        L--;
        if (p == -1 && graph[v].size())
        {
            graph[v].x = 0;
            graph[v].y = 0;
            graph[graph[v][0]].x = graph[v].x + L * cos(alpha + 2 * alpha * (st == 0 ? -1 : 1));
            graph[graph[v][0]].y = graph[v].y + L * sin(alpha + 2 * alpha * (st == 0 ? -1 : 1));
            coordinate_dfs(graph[v][0], v, !st, diameter, 1);
            return;
        }
        int used = 0;
        int next_diam = find(graph[v], diameter[dist]);
        if (next_diam != -1)
        {
            double fi = get_engle(graph[p], graph[v]);
            double _fi = 2 * pi - fi;
            graph[graph[v][next_diam]].x = graph[v].x + L * cos((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
            graph[graph[v][next_diam]].y = graph[v].y + L * sin((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
            coordinate_dfs(graph[v][next_diam], v, !st, diameter, dist + 1);
            used++;
        }
        for (int i = 0; i < graph[v].size(); ++i)
        {
            if (graph[v][i] == p || graph[v][i] == next_diam)
            {
                continue;
            }
            double fi = get_engle(graph[p], graph[v]);
            double _fi = 2 * pi - fi;
            if (used == 0)
            {
                graph[graph[v][i]].x = graph[v].x + L * cos((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
                graph[graph[v][i]].y = graph[v].y + L * sin((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
                coordinate_dfs(graph[v][i], v, !st, diameter, dist + 1);
            }
            else if (used == 1)
            {
                graph[graph[v][i]].x = graph[v].x + L * cos((pi / 2 + fi - 2 * alpha) * (!st) + (-pi / 2 - _fi + 2 * alpha) * (st));
                graph[graph[v][i]].y = graph[v].y + L * sin((pi / 2 + fi - 2 * alpha) * (!st) + (-pi / 2 - _fi + 2 * alpha) * (st));
                coordinate_dfs(graph[v][i], v, !st, diameter, dist + 1);
            }
            else if (used == 2)
            {

                graph[graph[v][i]].x = graph[v].x + L * cos((-pi / 2 + fi) * (!st) + (pi / 2 - _fi) * (st));
                graph[graph[v][i]].y = graph[v].y + L * sin((-pi / 2 + fi) * (!st) + (pi / 2 - _fi) * (st));
                coordinate_dfs(graph[v][i], v, st, diameter, dist + 1);
            }
        used++;
    }
}


// задает координаты атомам 
void give_coordinates()
{ 
    vector<int> diameter;
    pair<int, int> edge1 = {0, 0};
    pair<int, int> edge2 = {0, 0};
    diameter_dfs(root1, 0, -1, edge1);
    diameter_dfs(edge1.second, 0, -1, edge2);
    diameter = coordinates_bfs(edge1.second, edge2.second); // в diameter хранятся индексы вершин, принадлежащих одному из диаметров
    coordinate_dfs(edge2.second, -1, 0, diameter, 0);
    int x0 = 0, y0 = 0;
    for (int i = 0; i < this->size(); i++)
    {
        x0 = min(graph[i].x, x0);
        y0 = min(graph[i].y, y0);
    }
    for (int i = 0; i < this->size(); i++)
    {
        graph[i].x -= x0;
        graph[i].y -= y0;
    }
    for (int i = 0; i < this->size(); i++)
    {
        this->x1 = min(this->x1, graph[i].x);
        this->y1 = min(this->y1, graph[i].y);
        this->x2 = max(this->x2, graph[i].x);
        this->y2 = max(this->y2, graph[i].y);
    }
}

// возвращает количество атомов в молекуле
int size()
{
    return graph.size();
}

// возвращает ссылку i-тый атом
Atom &operator[](int i)
{
    return graph[i];
}


// возвращает строку с координатами атомов
string coordinates()
{ 
    string s = "";
    s += "[";
    for (int i = 0; i < size(); ++i)
    {
        s += "[" + to_string(graph[i].x) + ", " + to_string(graph[i].y) + "], ";
    }
    s += "]";
    return s;
}

// возвращает строку со списком смежности
string connectivity()
{ 
    string s = "";
    s += "[";
    for (int i = 0; i < size(); ++i)
    {
        s += "[";
        for (int j = 0; j < graph[i].size(); ++j)
        {
            s += to_string(graph[i][j]) + ", ";
        }
        s += "], ";
    }
    s += "]";
    return s;
}
};



class SVG_picture
{
public:
    string s = "";

    void draw_dfs(int v, int p, vector<int> &used, Molecule &mol)
    { // dfs для записи в s координат линий в формате для svg
        for (int i = 0; i < mol[v].size(); ++i)
        {
            if (mol[v][i] == p)
            {
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

    void read_molecule(Molecule &mol)
    { // заполняет файл координатами линий в формате svg
        vector<int> used(mol.size());
        s = "";
        s += "<svg xmlns=\"http://www.w3.org/2000/svg\">";
        s += '\n';
        draw_dfs(0, -1, used, mol);
        s += "</svg>";
        s += '\n';
    }
};



class Page
{
public:
    string s = "";
    int y_max = 0;

// dfs для записи в s координат линий в формате для svg
    void draw_dfs(int v, int p, Molecule &mol, int rx, int ry, int number)
    { 
        for (int i = 0; i < mol[v].size(); ++i)
        {
            if (mol[v][i] == p)
            {
                continue;
            }
            s += "<line x1=\"";
            s += to_string(mol[v].x + rx);
            s += "\" y1=\"";
            s += to_string(mol[v].y + ry);
            s += "\" x2=\"";
            s += to_string(mol[mol[v][i]].x + rx);
            s += "\" y2=\"";
            s += to_string(mol[mol[v][i]].y + ry);
            s += "\" ";
            s += "id=\"" + to_string(number) + "_" + to_string(v) + "_" + to_string(mol[v][i]) + "\" ";
            s += "/>";
            s += '\n';
            if(mol.chiral_atom[v]){
                s += "<circle cx=\"";
                s += to_string(mol[v].x + rx);
                s += "\" cy=\"";
                s += to_string(mol[v].y + ry);
                s += "\" r=\"5\" ";
                s += "id=\"" + to_string(number) + "_" + to_string(v) + "_" + to_string(mol[v][i]) + "_h" + "\" ";
                s += "/>";
                s += '\n';
            }
            draw_dfs(mol[v][i], v, mol, rx, ry, number);
        }
    }

// заполняет файл координатами линий в формате svg
    void read_molecules(vector<Molecule> &molecules)
    { 
        s = "";
        s += "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"UTF-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\"><title>Document</title></head>\n";
        s += "<style>#line {stroke: #000000;stroke-width: 1;}#back{background-color:#ffffff;}#theme{position:fixed;top:5%;left:90%;height:30px;width:100px;font-size:1em;background-color:#000000;}#small{position:fixed;top:10%;left:90%;height:30px;width:100px;font-size:1em;background-color:#482048;}#big{position: fixed;top: 15%;left: 90%;height: 30px;width: 100px;font-size: 1em;background-color:#00db56;}</style>\n";
        s += "<body id=\"back\">\n";
        s += "    <div id=\"line\">\n";
        int rx = 0, ry = 0;
        int dy = 0;
        s += "<svg height=\"\" width=\"1500\" xmlns=\"http://www.w3.org/2000/svg\"id=\"svg\">\n";
        s += '\n';
        for (int i = 0; i < molecules.size(); ++i)
        {
            if (molecules[i].x2 + rx > 1500)
            {
                rx = 0;
                ry += dy + 25;
                dy = 0;
            }
            draw_dfs(0, -1, molecules[i], rx, ry, i);
            dy = max(dy, molecules[i].y2);
            rx += molecules[i].x2;
            rx += 25;
            s += '\n';
            y_max = max(y_max, ry + dy);
        }
        s += "    </svg>\n";

        s += "</div><button id=\"theme\"></button>";
        s += "<button id=\"small\"></button>";
        s += "<button id=\"big\"></button>\n";
        s += "</body>\n";
        s += "<script>";
        s += "var k = 1;";
        s += "let molecules = [";
        for (int i = 0; i < molecules.size(); ++i)
        {
            s += "{";
            s += "chiral: [";
            for(int j = 0; j < molecules[i].chiral_atom.size(); ++j){
                s += to_string(molecules[i].chiral_atom[j]) + ", ";
            }
            s += "], \n";
            s += "x2: " + to_string(molecules[i].x2) + ",";
            s += "y2: " + to_string(molecules[i].y2) + ",";
            s += "coordinates: " + molecules[i].coordinates() + ",";
            s += "connectivity: " + molecules[i].connectivity();
            s += "},";
        }
        s += "];";
        s += "    function draw_dfs(v,p,mol,rx,ry,number){for(let i=0;i<mol.connectivity[v].length;++i){if(mol.connectivity[v][i]==p){continue;}let id=number.toString()+\"_\"+v.toString()+\"_\"+mol.connectivity[v][i].toString();let line=document.getElementById(id);line.setAttribute(\"x1\",(Math.floor(mol.coordinates[v][0]*k)+rx).toString());line.setAttribute(\"y1\",(Math.floor(mol.coordinates[v][1]*k)+ry).toString());line.setAttribute(\"x2\",(Math.floor(mol.coordinates[mol.connectivity[v][i]][0]*k)+rx).toString());line.setAttribute(\"y2\",(Math.floor(mol.coordinates[mol.connectivity[v][i]][1]*k)+ry).toString());if(mol.chiral[v]){let circle = document.getElementById(id + \"_h\");circle.setAttribute( \"cx\", (Math.floor(mol.coordinates[v][0] * k) + rx).toString());circle.setAttribute( \"cy\", (Math.floor(mol.coordinates[v][1] * k) + ry).toString());circle.setAttribute(\"r\",(Math.ceil(5 * k)).toString());} draw_dfs(mol.connectivity[v][i],v,mol,rx,ry,number);}}function change(){let y_max=0;let rx=0,ry=0;let dy=0;for(let i=0;i<molecules.length;++i){if(Math.floor(molecules[i].x2*k)+rx>1500){rx=0;ry+=dy+Math.floor(25*k);dy=0;}draw_dfs(0,-1,molecules[i],rx,ry,i);dy=Math.max(dy,Math.floor(molecules[i].y2*k));rx+=Math.floor(molecules[i].x2*k);rx+=Math.floor(25*k);y_max=Math.max(y_max,ry+dy);}document.getElementById(\"svg\").setAttribute(\"height\",y_max.toString());}document.getElementById(\"small\").onclick=function(){k-=0.1;change();};document.getElementById(\"big\").onclick=function(){k+=0.1;change();};";
        s += "document.getElementById(\"theme\").onclick = function () {if(window.getComputedStyle(document.getElementById(\"line\"), null).getPropertyValue(\"stroke\") ==  \"rgb(255, 255, 255)\"){document.getElementById(\"line\").style.stroke = \"#000000\";document.getElementById(\"back\").style.backgroundColor = \"white\";document.getElementById(\"theme\").style.backgroundColor = \"black\";}else{document.getElementById(\"line\").style.stroke = \"#ffffff\";document.getElementById(\"back\").style.backgroundColor = \"black\";document.getElementById(\"theme\").style.backgroundColor = \"white\";}};";
        s += "</script>";
        s += "</html>\n";
    }
};



class Unite_page
{
public:
    string s = "";

    Unite_page(int n)
    {
        s = "";
        s += "<!DOCTYPE html>\n";
        s += "<html lang=\"en\">\n";
        s += "<head>\n";
        s += "    <meta charset=\"UTF-8\">\n";
        s += "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
        s += "    <title>Document</title>\n";
        s += "    <style>\n";
        for (int i = 0; i < n; ++i)
        {
            s += "      #C";
            s += char((i + 1) / 10 + '0');
            s += char((i + 1) % 10 + '0');
            s += " {\n";
            s += "          position: fixed;\n";
            s += "          top: ";
            s += to_string(10 + i * 5);
            s += "%;\n";
            s += "          left: 40%;\n";
            s += "          height: 5%;\n";
            s += "          width: 20%;\n";
            s += "          font-size: 1em;\n";
            s += "      }\n";
        }
        s += "    </style>\n";
        s += "</head>\n";
        s += "<body>\n";

        for (int i = 0; i < n; ++i)
        {
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

mutex mtx;
void generate_molecule(vector<vector<Molecule>>& molecules, vector<Molecule>& prepear_molecules, int i, int j){
    for (int k = 0; k < molecules[i - 1][j].size(); ++k){
        if (molecules[i - 1][j][k].size() > 3){
            continue;
        }
        Molecule C(molecules[i - 1][j], k);
        mtx.lock();
        prepear_molecules.push_back(C);
        mtx.unlock();
    }
}



void clean_molecules(vector<Molecule>& prepear_molecules, vector<Molecule>& clear_molecules){
    bool already = false;
    clear_molecules.push_back(prepear_molecules[0]);
    for(int i = 1; i < prepear_molecules.size(); ++i){
        already = false;
        for(int j = 0; j < clear_molecules.size(); ++j){
            if(prepear_molecules[i] == clear_molecules[j]){
                already = true;
                break;
            }
        }
        if(!already){
            clear_molecules.push_back(prepear_molecules[i]);
        }
    }
}


class Molecule_container{
public:
    vector<vector<Molecule>> molecules;
    int n;
    int m;
    clock_t tStart;
    

    Molecule_container(int _n, int _m, int generation_type){
        molecules.resize(1);
        tStart = clock();
        n = _n;
        m = _m;
        Molecule C1;
        molecules[0].push_back(C1);

        if(generation_type == 0){
            step_generation();
        } else if(generation_type == 1){
            general_generation();
        }
        cout << "all time: " << 1.0 * (clock() - tStart) / CLOCKS_PER_SEC;
        cout << endl << endl;
    }

    void step_generation(){
        for (int i = 1; i < n; ++i){
            generate_layer_skelets(i);
            build_layer(i);
            draw_layer(i);
        }
    }


    void general_generation(){
        for (int i = 1; i < n; ++i){
            generate_layer_skelets(i);
        }
        for (int i = 1; i < n; ++i){
            build_layer(i);
        }
        for (int i = 1; i < n; ++i){
            draw_layer(i);
        }
    }


    void generate_layer_skelets(int i){
        molecules.resize(molecules.size() + 1);
        vector<Molecule> prepear_molecules;
        for (int j = 0; j < molecules[i - 1].size(); j += m){
            vector<thread> th;
            for(int l = 0; l < m && j + l < molecules[i - 1].size(); ++l){
                thread t(generate_molecule, ref(molecules), ref(prepear_molecules), i, j + l);
                th.push_back(move(t)); 
            }
            for(int l = 0; l < th.size(); ++l){
                th[l].join();
            }  
        }
        clean_molecules(prepear_molecules, molecules[i]);
    }


    void build_layer(int i){
        for (int j = 0; j < molecules[i].size(); j += m){
            vector<thread> th;
            for(int l = 0; l < m && j + l < molecules[i].size(); ++l){
                thread t(&Molecule::build, &molecules[i][j + l]);
                th.push_back(move(t)); 
            }
            for(int l = 0; l < th.size(); ++l){
                th[l].join();
            }  
        }
    }


    void draw_layer(int i){
        Page page;
        page.read_molecules(molecules[i]);
        page.s.insert(610, to_string(page.y_max));

        string s;
        ofstream file;

        s = ".\\pages\\C";
        s += char((i + 1) / 10 + '0');
        s += char((i + 1) % 10 + '0');
        s += ".html";
        file.open(s.c_str());
        file << page.s;
        file.close();

        Unite_page unite_page(i + 1);
        s = ".\\index.html";
        file;
        file.open(s.c_str());
        file << unite_page.s;
        file.close();
    
        molecules[i - 1].clear();
        
        cout << "C" << i + 1 << ": " << molecules[i].size() << " isomers" << "   ";
        cout << "time: " << 1.0 * (clock() - tStart)/CLOCKS_PER_SEC << endl;
    }

};









int main(){
    int n;
    int m;
    cout << "\n\nCreate a directory and name it \"pages\"\n\n";
    cout << "if you have created a directory, with name \"pages\" press Enter";
    cin.get();

    cout << "\nHow many threads do you want to use?\n";
    cout << "If you don't know what it is, you take 1\n";
    cin >> m;
    
    cout << "\nEnter the amount of carbon in the largest alkane\n";
    cin >> n;
    Molecule_container(n, m, 0);
    int x;
    cin >> x; 
}