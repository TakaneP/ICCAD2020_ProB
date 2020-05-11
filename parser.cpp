#include<bits/stdc++.h>
#define F first
#define S second
using namespace std;
 
struct master {
    set<pair<int,int>>pin;
    map<int,int>blkg;
    map<pair<int,int>,int>same;
    map<pair<int,int>,int>dif;
};
struct cell {
    int type;
    bool moveable;
    int x, y;
};
struct seg
{
    int sx,ex,sy,ey,sz,ez;
};
struct segtree
{
    int mi;
}node[7][55][205];
struct net {
    int min_lay;
    vector<pair<int,int>>net;
    vector<seg>seg;
};
int max_move;
int default_supply;
string str;
int row, col, lay;
int lay_lim[32];
vector<vector<vector<int>>>fix_supply;
vector<master>master_cell;
vector<cell>cell;
vector<net>net;
vector<vector<map<int,int>>>cell_cont;
int to_int(string &s, int k) {
    return stoi(s.substr(k)) - 1;
}
void pushup(int x, int lay, int pos) {
    node[lay][pos][x].mi = min(node[lay][pos][x<<1].mi, node[lay][pos][x<<1|1].mi);
}
void build(int x, int l, int r,int lay, int pos) {
    if (l == r) {
        if (lay % 2 == 0)
            node[lay][pos][x].mi = fix_supply[pos][l][lay];
        else
            node[lay][pos][x].mi = fix_supply[l][pos][lay];
        return;
    }
    int mid = (l + r) >> 1;
    build(x<<1, l, mid, lay, pos);
    build(x<<1|1, mid + 1, r, lay, pos);
    pushup(x, lay, pos);
}
int query(int x, int l, int r, int ql, int qr, int lay, int pos) {
    if (ql <= l && qr >= r) {
        return node[lay][pos][x].mi;
    }
    int mid = (l + r) >> 1;
    if (qr <= mid) {
        return query(x<<1, l, mid, ql, qr, lay, pos);
    }
    else if (ql > mid) {
        return query(x<<1|1, mid + 1, r, ql, qr, lay, pos);
    }
    else {
        return min(query(x<<1, l, mid, ql, qr, lay, pos), query(x<<1|1, mid + 1, r, ql, qr, lay, pos));
    }
}
int getsupply(int x1, int x2, int y1, int y2, int z) {
    if (z % 2 == 0) {
        return query(1, 1, col, y1, y2, z, x1);
    }
    else {
        return query(1, 1, row, x1, x2, z, y1);
    }
}
void build_ini_seg() {
    for (int i = 0 ; i < lay ; i++) {
        if (i % 2 == 0) {
            for (int j = 0 ; j < row ; j++) {
                build(1, 1, col, i, j);
            }
        }
        else {
            for (int j = 0 ; j < col ; j++) {
                build(1, 1, row, i, j);
            }
        }
    }
}
void add_mc(int x, int y, int tp) {
    for (auto i:master_cell[tp].same) {
        if (cell_cont[x][y].count(i.F.F)) {
            fix_supply[x][y][i.F.S] -= cell_cont[x][y][i.F.F] * i.S;
        }
    }
    for (auto i:master_cell[tp].dif) {
        if (i.F.S % 2 == 0) {
            if (y >= 1 && cell_cont[x][y-1].count(i.F.F)) {
                fix_supply[x][y][i.F.S] -= cell_cont[x][y-1][i.F.F] * i.S;
                fix_supply[x][y-1][i.F.S] -= cell_cont[x][y-1][i.F.F] * i.S;
            }
            if (y < col - 1 && cell_cont[x][y+1].count(i.F.F)) {
                fix_supply[x][y][i.F.S] -= cell_cont[x][y+1][i.F.F] * i.S;
                fix_supply[x][y+1][i.F.S] -= cell_cont[x][y+1][i.F.F] * i.S;
            }
        }
        else {
            if (x >= 1 && cell_cont[x-1][y].count(i.F.F)) {
                fix_supply[x][y][i.F.S] -= cell_cont[x-1][y][i.F.F] * i.S;
                fix_supply[x-1][y][i.F.S] -= cell_cont[x-1][y][i.F.F] * i.S;
            }
            if (x < row - 1 && cell_cont[x+1][y].count(i.F.F)) {
                fix_supply[x][y][i.F.S] -= cell_cont[x+1][y][i.F.F] * i.S;
                fix_supply[x+1][y][i.F.S] -= cell_cont[x+1][y][i.F.F] * i.S;
            }
        }
        
    }
    cell_cont[x][y][tp]++;
}
void del_mc(int x, int y, int tp) {
    cell_cont[x][y][tp]--;
    if (cell_cont[x][y][tp] == 0) {
        cell_cont[x][y].erase(tp);
    }
    for (auto i:master_cell[tp].blkg) {
        fix_supply[x][y][i.F] += i.S;
    }
    for (auto i:master_cell[tp].same) {
        if (cell_cont[x][y].count(i.F.F)) {
            fix_supply[x][y][i.F.S] += cell_cont[x][y][i.F.F] * i.S;
        }
    }
    for (auto i:master_cell[tp].dif) {
        if (i.F.S % 2 == 0) {
            if (y >= 1 && cell_cont[x][y-1].count(i.F.F)) {
                fix_supply[x][y][i.F.S] += cell_cont[x][y-1][i.F.F] * i.S;
                fix_supply[x][y-1][i.F.S] += cell_cont[x][y-1][i.F.F] * i.S;
            }
            if (y < col - 1 && cell_cont[x][y+1].count(i.F.F)) {
                fix_supply[x][y][i.F.S] += cell_cont[x][y+1][i.F.F] * i.S;
                fix_supply[x][y+1][i.F.S] += cell_cont[x][y+1][i.F.F] * i.S;
            }
        }
        else {
            if (x >= 1 && cell_cont[x-1][y].count(i.F.F)) {
                fix_supply[x][y][i.F.S] += cell_cont[x-1][y][i.F.F] * i.S;
                fix_supply[x-1][y][i.F.S] += cell_cont[x-1][y][i.F.F] * i.S;
            }
            if (x < row - 1 && cell_cont[x+1][y].count(i.F.F)) {
                fix_supply[x][y][i.F.S] += cell_cont[x+1][y][i.F.F] * i.S;
                fix_supply[x+1][y][i.F.S] += cell_cont[x+1][y][i.F.F] * i.S;
            }
        }  
    }
}
int main()
{
    freopen("case3.txt","r",stdin);
    int k;
    //Cell Max Move
    cin >> str >> max_move;
    //Grid Boundary
    cin >> str >> k >> k >> k;
    row = k;
    cin >> k;
    col = k;
    //Layer
    cin >> str >> lay;
    for (int i = 0 ; i < lay ; i++) {
        cin >> str >> str >> str >> str >> lay_lim[i];
    }
    //Initial each layer's routing supply
    for (int i = 0 ; i < row ; i++) {
        vector<vector<int>>add;
        vector<map<int,int>>add2;
        fix_supply.push_back(add);
        cell_cont.push_back(add2);
        for (int j = 0 ; j < col ; j++) {
            vector<int>add;
            map<int,int>add2;
            fix_supply[i].push_back(add);
            cell_cont[i].push_back(add2);
            for (int k = 0 ; k < lay ; k++) {
                fix_supply[i][j].push_back(lay_lim[k]);
            }
        }
    }
    //Non Default Supply
    cin >> str >> k;
    for (int i = 0 ; i < k ; i++) {
        int x, y, z;
        cin >> x >> y >> z;
        x--, y--, z--;
        int del;
        cin >> del;
        fix_supply[x][y][z] += del;
    }
    //Master Cell
    cin >> str >> k;
    master_cell.resize(k);
    for (int i = 0 ; i < k ; i++) {
        int p, b;
        cin >> str >> str >> p >> b;
        for (int j = 0 ; j < p ; j++) {
            string s1, s2;
            cin >> str >> s1 >> s2;
            int x = to_int(s1, 1); //to_int will sub 1 and return 
            int y = to_int(s2, 1);
            master_cell[i].pin.insert({x,y});
        }
        for (int j = 0 ; j < b ; j++) {
            string s1, s2;
            int c;
            cin >> str >> s1 >> s2 >> c;
            //int x = to_int(s1, 1); Blockage name can be discarded
            int y = to_int(s2, 1);
            master_cell[i].blkg[y] += c; //Calculate the demand caused by the blockage on this master cell, same layer will be summed up
        }
    }
    //Extra Demand
    cin >> str >> k;
    for (int i = 0 ; i < k ; i++) {
        bool op;
        int x, y, pen, pos; //MC1 MC2 penalty layer
        cin >> str;
        op = (str[0] == 'a'); //same: 0 , adj : 1
        cin >> str;
        x = to_int(str, 2);
        cin >> str;
        y = to_int(str, 2);
        cin >> str;
        pos = to_int(str, 1);
        cin >> pen;
        if (op == 0){
            master_cell[x].same[{y,pos}] += pen;
            //master_cell[y].same[{x,pos}] += pen;
        }
        else {
            master_cell[x].dif[{y,pos}] += pen;
            //master_cell[y].dif[{x,pos}] += pen; Reverse when calculating
        }
    }
    //Cell Instance
    cin >> str >> k;
    for (int i = 0 ; i < k ; i++) {
        int x, y, tp;
        cin >> str >> str >> str;
        tp = to_int(str,2); //MC
        cin >> x >> y >> str;
        x--, y--;
        bool mov = (str[0] == 'M');
        cell.push_back({tp,mov,x,y});
        add_mc(x,y,tp); //Update graph
    }
    //Net
    cin >> str >> k;
    for (int i = 0 ; i < k ; i++) {
        int x;
        cin >> str >> str >> x >> str;
        int Min = 0;
        if (str[0] != 'N') {
            Min = to_int(str,1); //MinRoutingConstraint
        }
        vector<pair<int,int>>add;
        for (int j = 0 ; j < x ; j++) {
            cin >> str;
            cin >> str;
            size_t pos = str.find('/');
            string s1, s2;
            s1 = str.substr(0,pos);
            s2 = str.substr(pos+1);
            int x = to_int(s1,1);
            int y = to_int(s2,1);
            add.push_back({x,y});
        } 
        net.push_back({Min,add,{}});
    }
    //Build segment tree
    build_ini_seg();
    //Initial routing segments
    cin >> str >> k;
    for (int i = 0 ; i < k ; i++) {
        int x1, x2 , y1, y2, z1, z2;
        cin >> x1 >> y1 >> z1 >> x2 >> y1 >> z2 >> str;
        int x = to_int(str,1);
        net[x].seg.push_back({min(x1,x2)-1,max(x1,x2)-1,min(y1,y1)-1,max(y1,y2)-1,min(z1,z2)-1,max(z1,z2)-1});
    }

}
