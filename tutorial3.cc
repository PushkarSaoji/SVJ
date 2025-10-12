// sift_clustering_opt_with_rules.cpp
// g++ -O3 -march=native -std=c++17 -funroll-loops -pipe sift_clustering_opt_with_rules.cpp -o sift_opt
// Optional: -DUSE_SIMD -fopenmp

#include <bits/stdc++.h>
#include "Pythia8/Pythia.h"
#include "chrono"
#include "iostream"

using namespace std;
using int32 = int32_t;
using uint8  = uint8_t;


// ---------------- Physics SoA ----------------

struct FourVecSoA {
    vector<double> E, px, py, pz;
    FourVecSoA() = default;
    FourVecSoA(size_t n) { resize(n); }
    void resize(size_t n) {
        E.assign(n,0.0); px.assign(n,0.0); py.assign(n,0.0); pz.assign(n,0.0);
    }
    inline void set(size_t i,double Ee,double pxi,double pyi,double pzi){
        E[i]=Ee; px[i]=pxi; py[i]=pyi; pz[i]=pzi;
    }
    inline void add_into(size_t dst,size_t src){
        E[dst]+=E[src]; px[dst]+=px[src]; py[dst]+=py[src]; pz[dst]+=pz[src];
    }
};

inline double et2_comp(double E,double pz){ return E*E - pz*pz; }
inline double ET_comp(double E,double pz){ double v=et2_comp(E,pz); return v>0.0? sqrt(v) : 0.0; }
inline double dot_comp(double E1,double px1,double py1,double pz1,double E2,double px2,double py2,double pz2){
    return E1*E2 - px1*px2 - py1*py2 - pz1*pz2;
}
inline double siftDist_comp(double E1,double px1,double py1,double pz1,
                            double E2,double px2,double py2,double pz2){
    double num = 2.0 * dot_comp(E1,px1,py1,pz1, E2,px2,py2,pz2);
    double denom = et2_comp(E1,pz1) + et2_comp(E2,pz2);
    return num / denom;
}

// ---------------- History/event types ----------------

enum class Action { MERGE, DROP, FILTER };
struct MergeStep {
    int32 left_id, right_id;
    int32 new_id;
    double dist;
    int32 step;
    Action action;
};

// ---------------- Heap entry ----------------

struct HeapEntry {
    double dist;
    int32 i,j,ver;
    bool operator<(const HeapEntry& o) const { return dist > o.dist; } // reverse for min-heap
};

// ---------------- SIFT clustering optimized with rules ----------------

struct SIFTClusteringOptRules {
    int32 n0;
    int32 n_alloc;
    FourVecSoA vecs;
    vector<uint8> alive;    // 1 alive, 0 removed (dropped), but final_obj flagged separately
    vector<uint8> final_obj; // 1 means declared final by FILTER
    vector<int32> id_map;
    vector<int32> version;
    vector<double> dist; // flat row-major n_alloc * n_alloc
    vector<int32> nn;
    vector<double> nnDist;
    vector<HeapEntry> heap;
    vector<MergeStep> history;
    int32 next_new_id;
    int32 alive_count;

    SIFTClusteringOptRules(const vector<array<double,4>>& parts) {
        n0 = int32(parts.size());
        n_alloc = n0;
        vecs.resize(n_alloc);
        alive.assign(n_alloc,1);
        final_obj.assign(n_alloc,0);
        id_map.resize(n_alloc);
        version.assign(n_alloc,0);
        dist.assign(size_t(n_alloc)*size_t(n_alloc), 0.0);
        nn.assign(n_alloc,-1);
        nnDist.assign(n_alloc,1e300);
        next_new_id = n0;
        alive_count = n0;
        for (int32 i=0;i<n_alloc;i++){
            vecs.set(i, parts[i][0], parts[i][1], parts[i][2], parts[i][3]);
            id_map[i]=i;
        }
        compute_all_distances_and_nns();
        build_heap_from_nns();
    }

    inline size_t idx(int32 i,int32 j) const { return size_t(i)*size_t(n_alloc) + size_t(j); }

    void compute_all_distances_and_nns(){
        const int32 N = n_alloc;
        for (int32 i=0;i<N;i++){
            double Ei = vecs.E[i], pxi = vecs.px[i], pyi = vecs.py[i], pzi = vecs.pz[i];
            for (int32 j=i+1;j<N;j++){
                double d = siftDist_comp(Ei,pxi,pyi,pzi, vecs.E[j],vecs.px[j],vecs.py[j],vecs.pz[j]);
                dist[idx(i,j)] = d;
                dist[idx(j,i)] = d;
            }
            dist[idx(i,i)] = 0.0;
        }
        for (int32 i=0;i<N;i++){
            if (!alive[i]) { nn[i]=-1; nnDist[i]=1e300; continue;}
            if (final_obj[i]) { nn[i]=-1; nnDist[i]=1e300; continue; }
            double best=1e300; int32 arg=-1;
            const double *row = &dist[idx(i,0)];
            #ifdef USE_SIMD
            #pragma omp simd
            #endif
            for (int32 j=0;j<N;j++){
                if (i==j) continue;
                if (!alive[j] || final_obj[j]) continue;
                double d = row[j];
                if (d < best) { best=d; arg=j; }
            }
            nn[i]=arg; nnDist[i]=best;
        }
    }

    void build_heap_from_nns(){
        heap.clear(); heap.reserve(n_alloc);
        for (int32 i=0;i<n_alloc;i++){
            if (!alive[i] || final_obj[i]) continue;
            if (nn[i] != -1) heap.push_back(HeapEntry{nnDist[i], i, nn[i], version[i]});
        }
        make_heap(heap.begin(), heap.end());
    }

    // decide action according to SIFT rules for pair (i,j)
    pair<Action,int32> decide_action_pair(int32 i,int32 j){
        double Ei = vecs.E[i], pzi = vecs.pz[i];
        double Ej = vecs.E[j], pzj = vecs.pz[j];
        double ETi = ET_comp(Ei,pzi), ETj = ET_comp(Ej,pzj);
        double denom = ETi*ETi + ETj*ETj;
        double eps = 0.0;
        if (denom > 0.0) eps = (ETi*ETj) / denom;
        double c2 = (2.0*eps)*(2.0*eps);
        double delta = dist[idx(i,j)];
        // Use consistent interpretation:
        if (delta < c2 && c2 <= 1.0) return {Action::MERGE, -1};
        if (c2 <= 1.0 && delta >= c2 && delta < 1.0) {
            // drop lower pT (pT^2 comparison)
            double pTi2 = vecs.px[i]*vecs.px[i] + vecs.py[i]*vecs.py[i];
            double pTj2 = vecs.px[j]*vecs.px[j] + vecs.py[j]*vecs.py[j];
            return {Action::DROP, (pTi2 <= pTj2 ? i : j)};
        }
        if (delta >= 1.0) return {Action::FILTER, -1};
        return {Action::MERGE, -1};
    }

    void recompute_nn_for(int32 i){
        if (!alive[i] || final_obj[i]) { nn[i] = -1; nnDist[i]=1e300; return; }
        double best=1e300; int32 arg=-1;
        const double *row = &dist[idx(i,0)];
        const int32 N = n_alloc;
        #ifdef USE_SIMD
        #pragma omp simd
        #endif
        for (int32 j=0;j<N;j++){
            if (i==j) continue;
            if (!alive[j] || final_obj[j]) continue;
            double d = row[j];
            if (d < best) { best = d; arg = j; }
        }
        nn[i]=arg; nnDist[i]=best;
    }

    void run_all(){
        int32 step = 0;
        while (true) {
            // stop condition: <=1 clusters eligible to cluster remain
            int32 eligible = 0;
            for (int32 i=0;i<n_alloc;i++) if (alive[i] && !final_obj[i]) ++eligible;
            if (eligible <= 1) break;

            // pop valid heap entry
            HeapEntry top;
            bool found=false;
            while (!heap.empty()) {
                pop_heap(heap.begin(), heap.end()); top = heap.back(); heap.pop_back();
                if (!alive[top.i] || !alive[top.j]) continue;
                if (final_obj[top.i] || final_obj[top.j]) continue;
                if (version[top.i] != top.ver) continue;
                if (nn[top.i] != top.j) continue;
                found = true; break;
            }
            if (!found) {
                // recompute all nns and rebuild
                compute_all_distances_and_nns();
                build_heap_from_nns();
                continue;
            }

            int32 a = top.i, b = top.j;
            if (!alive[a] || !alive[b]) continue;

            auto [action, drop_idx] = decide_action_pair(a,b);
            double dval = dist[idx(a,b)];

            if (action == Action::MERGE) {
                int32 newID = next_new_id++;
                history.push_back(MergeStep{ id_map[a], id_map[b], newID, dval, step, Action::MERGE });
                // merge b into a
                vecs.add_into(a, b);
                id_map[a] = newID;
                alive[b] = 0;
                version[a]++; --alive_count;

                // update row/col for a
                const int32 N = n_alloc;
                double Ea = vecs.E[a], pxa = vecs.px[a], pya = vecs.py[a], pza = vecs.pz[a];
                for (int32 k=0;k<N;k++){
                    if (k==a) { dist[idx(a,k)] = 0.0; continue; }
                    if (!alive[k] || final_obj[k]) { dist[idx(a,k)] = dist[idx(k,a)] = 1e300; continue; }
                    double nd = siftDist_comp(Ea,pxa,pya,pza, vecs.E[k],vecs.px[k],vecs.py[k],vecs.pz[k]);
                    dist[idx(a,k)] = nd; dist[idx(k,a)] = nd;
                }
                recompute_nn_for(a);
                if (nn[a] != -1) { heap.push_back(HeapEntry{nnDist[a], a, nn[a], version[a]}); push_heap(heap.begin(), heap.end()); }
                for (int32 k=0;k<N;k++){
                    if (!alive[k] || final_obj[k] || k==a) continue;
                    if (nn[k]==a || nn[k]==b) {
                        recompute_nn_for(k);
                        if (nn[k] != -1) { heap.push_back(HeapEntry{nnDist[k], k, nn[k], version[k]}); push_heap(heap.begin(), heap.end()); }
                    } else {
                        double dak = dist[idx(k,a)];
                        if (dak < nnDist[k]) {
                            nn[k] = a; nnDist[k] = dak;
                            heap.push_back(HeapEntry{nnDist[k], k, nn[k], version[k]}); push_heap(heap.begin(), heap.end());
                        }
                    }
                }
            }
            else if (action == Action::DROP) {
                int32 drop = drop_idx;
                history.push_back(MergeStep{ id_map[drop], -1, -1, dval, step, Action::DROP });
                alive[drop] = 0; --alive_count;
                // update neighbors which pointed to drop
                for (int32 k=0;k<n_alloc;k++){
                    if (!alive[k] || final_obj[k]) continue;
                    if (nn[k] == drop) {
                        recompute_nn_for(k);
                        if (nn[k] != -1) { heap.push_back(HeapEntry{nnDist[k], k, nn[k], version[k]}); push_heap(heap.begin(), heap.end()); }
                    }
                }
            }
            else { // FILTER
                history.push_back(MergeStep{ id_map[a], id_map[b], -1, dval, step, Action::FILTER });
                final_obj[a] = 1; final_obj[b] = 1;
                // recompute nn for others that pointed to a or b
                for (int32 k=0;k<n_alloc;k++){
                    if (!alive[k] || final_obj[k]) continue;
                    if (nn[k] == a || nn[k] == b) {
                        recompute_nn_for(k);
                        if (nn[k] != -1) { heap.push_back(HeapEntry{nnDist[k], k, nn[k], version[k]}); push_heap(heap.begin(), heap.end()); }
                    }
                }
            }
            ++step;
        }
    }
};

// ---------------- Demo ----------------

int main(){
    ofstream myfile;
    myfile.open("heap_time.csv");

    Pythia8::Pythia pythia;
    pythia.readFile("QCD.dat");
    
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowInfo = 0");
    
    pythia.init();
    int nEvent=1000;
    int N;
    
    double cluster_time;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;
        
        vector<array<double,4>> parts;
        
        for (int i=1; i< pythia.event.size(); i++){
            if (pythia.event[i].status() > 0){
                if (abs(pythia.event[i].eta())>5) continue;
                if (!pythia.event[i].isVisible()) continue;
                parts.push_back({pythia.event[i].e(),pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz()});
            }
        }
        N=parts.size();
        auto t1 = chrono::high_resolution_clock::now();
        SIFTClusteringOptRules algo(parts);
        algo.run_all();
        auto t2 = chrono::high_resolution_clock::now();
        cluster_time = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
        
        //    cerr << "Performed " << algo.history.size() << " events\n";
        //    for (size_t i=0;i<algo.history.size(); ++i){
            //       auto &h = algo.history[i];
            //       string s = (h.action==Action::MERGE?"MERGE":(h.action==Action::DROP?"DROP":"FILTER"));
            //       cerr << "step " << h.step << ": " << s << " (" << h.left_id << "," << h.right_id << ") -> new=" << h.new_id
     //        << " dist=" << h.dist << "\n";
     //   }
    myfile << N << ", " << cluster_time << endl;
    }
    return 0;
}
