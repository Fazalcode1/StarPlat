#include <bits/stdc++.h>
#include <omp.h>
#include"../graph.hpp"

using namespace std;
    
vector<vector<vector<double>>> prims_mst(int V,graph g,propNode<pair<double,double>> points,propNode<double> Demand,propEdge<int> weight){
        int v[V];
        memset(v,-1,sizeof(v));
        set<pair<double,pair<double,double>>> s;
        s.insert({0,{0,-1}});
        vector<vector<vector<double>>> z(V,vector<vector<double>>());
        while(!s.empty()){
            auto it=s.begin();
            s.erase(it);
            auto k = *it;
            auto y = k.second;
            if(v[(int)y.first]==1)continue;
                v[(int)y.first]=1;
               if(y.second!=-1){
                z[y.first].push_back({y.second,k.first});
                z[y.second].push_back({y.first,k.first});
                }
                // for(auto x:adj[y.first]){
                //    if(v[(int)x[0]]==-1)s.insert({x[1],{x[0],y.first}});
                // }
                forall(x in g.neighbors(y.first)){
                    if(v[x]==-1){
                        edge e = g.get_edge(x,y.first);
                        s.insert({e.weight,{x,y.first}});
                    }
                }
        }
        return z;
}

vector<vector<vector<double>>> randomise(vector<vector<vector<double>>> z){
     for(int i=0;i<z.size();i++){
       auto rng = default_random_engine {};
       shuffle(begin(z[i]), end(z[i]), rng);
     }
     return z;
}

vector<int> dfs_visit(int node, int visited[], const vector<vector<vector<double>>>& z) {
    visited[node] = 1;
    vector<int> permutation;
    permutation.push_back(node);
    int index = 0;
    while (index < permutation.size()) {
        int current = permutation[index++];
        for (auto x : z[current]) {
            if (visited[(int)x[0]] == -1) {
                visited[(int)x[0]] = 1;
                permutation.push_back((int)x[0]);
            }
        }
    }
    return permutation;
}


vector<vector<int>> convert_to_roots(int Depot,vector<int> p, double Demand[],double Capacity){
    vector<vector<int>> Routes;
    vector<int> OneRoute;
    int ResidueCap = Capacity;
    for(auto v:p){
        if(v==Depot)continue;
        if(ResidueCap-Demand[v]>=0){
            OneRoute.push_back(v);
            ResidueCap -= Demand[v];
        }
        else {
            Routes.push_back(OneRoute);
            OneRoute.clear();
            OneRoute.push_back(v);
            ResidueCap = Capacity-Demand[v];
        }
    }
    Routes.push_back(OneRoute);
    return Routes;
}

double distance(pair<int,int>p1,pair<int,int>p2){
    double d;
    d = (double)sqrt((p1.first-p2.first)*(p1.first-p2.first)+(p1.second-p2.second)*(p1.second-p2.second));
    return d;
}

double cost_of_oneroute(vector<int> OneRoute,int node,vector<pair<double,double>> c){
    double cost=0;
    cost+=distance(c[node],c[OneRoute[0]]);
    int previous_node=OneRoute[0];
    for(auto current_node:OneRoute){
            cost+=distance(c[previous_node],c[current_node]);
            previous_node=current_node;
    }
    cost+=distance(c[previous_node],c[node]);
    return cost;
}

double cost(vector<vector<int>> Routes,int node,vector<pair<double,double>> c){
    double cost=0;
    for(auto OneRoute:Routes){
        cost+=cost_of_oneroute(OneRoute,0,c);
    }
    return cost;
}

vector<vector<int>> Nearest_Neighbour(vector<vector<int>> Routes,int V,vector<pair<double,double>> c){
    int visited[V];
    memset(visited,-1,sizeof(visited));
    vector<vector<int>> Modified_Routes;
    for(auto OneRoute:Routes){
        int n=OneRoute.size();
        vector<int> tour;
        int current = 0;
        visited[current]=1;
        while(tour.size()<n){
            int nearest=-1;
            double minDist = numeric_limits<double>::max();
            for(auto i:OneRoute){
                if(visited[i]==-1){
                    double dist = distance(c[i],c[current]);
                    if(dist<minDist){
                        minDist = dist;
                        nearest = i;
                    }
                }
            }
            current = nearest;
            visited[current] = 1;
            tour.push_back(current);
        }
        Modified_Routes.push_back(tour);
    }
    return Modified_Routes;     
}

vector<vector<int>> two_OPT(vector<vector<int>> Routes,int V,vector<pair<double,double>> p){
    vector<vector<int>> Modified_Routes;
    for(auto OneRoute:Routes){
        vector<int> tour,t;
        tour.push_back(0);
        for(auto x:OneRoute)tour.push_back(x);
        int n = tour.size();
        bool improvement = true;
        while (improvement) {
            improvement = false;
            for (int i = 0; i < n - 1; ++i) {
                for (int j = i + 1; j < n; ++j) {
                int a = tour[i];
                int b = tour[(i + 1) % n];
                int c = tour[j];
                int d = tour[(j + 1) % n];
                double distBefore = distance(p[a],p[b]) + distance(p[c], p[d]);
                double distAfter = distance(p[a], p[c]) + distance(p[b], p[d]);
                if (distAfter < distBefore) {
                    reverse(tour.begin()+(i + 1),tour.begin()+j+1);
                    improvement = true;
                }
                }
            }
        }  
        int p;
        for(int i=0;i<tour.size();i++){
            if(tour[i]==0){
                p=i;
                break;
            }
        }
        for(int i=p+1;i<tour.size();i++){
            t.push_back(tour[i]);
        }
        for(int i=0;i<p;i++){
            t.push_back(tour[i]);
        }
        Modified_Routes.push_back(t);
    }
    return Modified_Routes;
}

vector<vector<int>> Refine_Routes(vector<vector<int>> Routes,int V,vector<pair<double,double>> c){
    vector<vector<int>> Modified_Routes,RoutesOne,RoutesTwo,RoutesThree;
    RoutesOne = Nearest_Neighbour(Routes,V,c);
    RoutesTwo = two_OPT(RoutesOne,V,c);
    RoutesThree = two_OPT(Routes,V,c);
    for(int i=0;i<Routes.size();i++){
        if(cost_of_oneroute(RoutesTwo[i],0,c)>=cost_of_oneroute(RoutesThree[i],0,c)){
            Modified_Routes.push_back(RoutesThree[i]);
        }
        else Modified_Routes.push_back(RoutesTwo[i]);
    }
    return Modified_Routes;
}

int main(){
    graph g;
    int V,Capacity;
    cin>>V>>Capacity;
    for(int i=0;i<V;i++){
        g.add_node(i);
    }
    for(int i=0;i<V;i++){
        for(int j=0;j<V;j++){
            if(j==i)continue;
            g.add_edge([i,j]);
        }
    }
    string s1,s2,s3,s4;
    cin>>s1;
    double Dem[V];
    vector<pair<double,double>> po(V);
    propNode<pair<double,double>> points;
    propNode<double> Demand;
    propEdge<int> weight; /*declaring the edge property*/
    g.attachEdgeProperty(weight);
    g.attachNodeProperty(points);
    g.attachNodeProperty(Demand);
    // vector<pair<double,double>> points(V);
    // double Demand[V];
    for(int i=0;i<V;i++){
        double n,x,y;
        cin>>n>>x>>y;
        node k=n-1;
        k.points = {x,y};
        po[n-1]={x,y};
        // points[node-1]={x,y};
    }
    cin>>s2;
    for(int i=0;i<V;i++){
        int n,demand;
        cin>>n>>demand;
        node w=n-1;
        // Demand[node-1]=demand;
        w.Demand = demand;
        Dem[node-1]=demand;
    }
    cin>>s3;
    int depot;
    cin>>depot;
    depot--;
    forall(v in g.nodes()){
        forall(nbr in g.nodes_from(v)){
            edge e = get_edge(v,nbr);
            e.weight = distance(v.points,nbr.points);
        }
    }
    vector<vector<vector<double>>> z(V,vector<vector<double>>());
    int i=0;
    // for(int i=0;i<V;i++){
    //     for(int j=0;j<V;j++){
    //         if(j==i)continue;
    //         double w = distance(points[i],points[j]);
    //         z[i].push_back({(double)j,w});
    //         z[j].push_back({(double)i,w});
    //     }
    // }

    z = prims_mst(V,graph g,points,Demand,weight);
    double minCost = INT_MAX;
    vector<vector<int>> minRoutes;
    for(int i=1;i<=pow(10,5);i++){
        vector<int> permutation;
        int visited[V];
        memset(visited,-1,sizeof(visited));
        z = randomise(z);
        permutation = dfs_visit(depot,visited,z);
        vector<vector<int>> Routes = convert_to_roots(depot,permutation,Dem,Capacity);
        double Cost = cost(Routes,depot,po);
        if(Cost<minCost){
            minCost = Cost;
            minRoutes = Routes;
        }   
    }
    minRoutes = Refine_Routes(minRoutes,V,po);
    minCost = cost(minRoutes,depot,po);
    cout<<"minimum cost = "<<minCost<<endl;
    // cout<<"minimum cost routes from depot:"<<endl;
    // int t=1;
    // for(auto x:minRoutes){
    //     cout<<"Route "<<t<<": ";
    //     for(auto y:x){
    //         cout<<y<<" ";
    //     }
    //     cout<<endl;
    //     t++;
    // }
}
