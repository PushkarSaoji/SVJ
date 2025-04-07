#include "Pythia8/Pythia.h"
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
#include <float.h>
#include <chrono>

using namespace Pythia8;

// Custom dot product function using Pythia8’s Particle four‐vector methods.
double dot(const Particle &p1, const Particle &p2) {
    // Compute the modified dot product:
    double prod = -p1.px() * p2.px() - p1.py() * p2.py() - p1.pz() * p2.pz() + p1.e() * p2.e();
    return 2 * prod;
}

// Updated custom distance function (del_ab) using the new metric: (del_ab)^0.03
double del_ab(const Particle &p1, const Particle &p2) {
    // Assume pT() is defined in Particle (or computed as sqrt(px^2+py^2))
    double metric = dot(p1, p2) / ( std::pow(p1.pT(), 2) + std::pow(p2.pT(), 2) );
    return std::pow(metric, 0.03);
}

// VP-tree node storing a pointer to a Particle.
struct VPNode {
    const Particle* point; 
    double threshold;
    std::unique_ptr<VPNode> inner;
    std::unique_ptr<VPNode> outer;

    VPNode(const Particle* pt) : point(pt), threshold(0.0), inner(nullptr), outer(nullptr) {}
};

// VP-tree class for Particle objects using our custom metric.
class VPTree {
public:
    // Build the tree from a vector of Particle pointers.
    VPTree(std::vector<const Particle*>& items) {
        root = build(items);
    }

    // Returns a pair: (pointer to nearest neighbor, distance to it)
    std::pair<const Particle*, double> nearest_neighbor(const Particle* query) {
        return search(root.get(), query, nullptr, std::numeric_limits<double>::max());
    }

private:
    std::unique_ptr<VPNode> root;
    std::mt19937 rng{std::random_device{}()};

    // Recursively build the VP-tree.
    std::unique_ptr<VPNode> build(std::vector<const Particle*>& items) {
        if (items.empty())
            return nullptr;

        // Choose a random vantage point.
        std::uniform_int_distribution<size_t> dist(0, items.size() - 1);
        size_t idx = dist(rng);
        const Particle* vp = items[idx];
        items.erase(items.begin() + idx);

        auto node = std::make_unique<VPNode>(vp);

        if (items.empty())
            return node;

        // Compute distances from vp to all other items.
        std::vector<double> distances;
        distances.reserve(items.size());
        for (const Particle* item : items) {
            distances.push_back(del_ab(*vp, *item));
        }

        // Find the median distance.
        std::vector<double> sorted = distances; // copy for median calculation
        size_t mid = sorted.size() / 2;
        std::nth_element(sorted.begin(), sorted.begin() + mid, sorted.end());
        double median = sorted[mid];
        node->threshold = median;

        // Partition items into inner (closer than median) and outer (farther than or equal to median).
        std::vector<const Particle*> inner_items;
        std::vector<const Particle*> outer_items;
        for (size_t i = 0; i < items.size(); i++) {
            if (distances[i] < median)
                inner_items.push_back(items[i]);
            else
                outer_items.push_back(items[i]);
        }

        // Recursively build subtrees.
        node->inner = build(inner_items);
        node->outer = build(outer_items);

        return node;
    }

    // Recursive nearest neighbor search.
    std::pair<const Particle*, double> search(VPNode* node, const Particle* query,
                                                const Particle* best, double best_dist) {
        if (!node)
            return { best, best_dist };

        double d = del_ab(*query, *(node->point));
        // Do not compare the query with itself.
        if (query != node->point && d < best_dist) {
            best = node->point;
            best_dist = d;
        }

        if (!node->inner && !node->outer)
            return { best, best_dist };

        // Decide which subtree to search first.
        if (d < node->threshold) {
            auto inner_result = search(node->inner.get(), query, best, best_dist);
            best = inner_result.first;
            best_dist = inner_result.second;
            if (d + best_dist >= node->threshold) {
                auto outer_result = search(node->outer.get(), query, best, best_dist);
                best = outer_result.first;
                best_dist = outer_result.second;
            }
        } else {
            auto outer_result = search(node->outer.get(), query, best, best_dist);
            best = outer_result.first;
            best_dist = outer_result.second;
            if (d - best_dist <= node->threshold) {
                auto inner_result = search(node->inner.get(), query, best, best_dist);
                best = inner_result.first;
                best_dist = inner_result.second;
            }
        }
        return { best, best_dist };
    }
};

void search(VPTree& tree, std::vector<const Particle*>& particle_ptrs, const Particle*& best_pair_first, const Particle*& best_pair_second, double& global_best){
   // For each particle, search for its nearest neighbor.
    for (const Particle* p : particle_ptrs) {
        auto res = tree.nearest_neighbor(p);
        if (res.first != nullptr && res.second < global_best) {
            global_best = res.second;
            best_pair_first = p;
            best_pair_second = res.first;
        }
    }
}

int main() {
    // For demonstration, we create a few dummy Particle objects.
    // In actual use with Pythia8, these Particle objects would be produced by the generator.
    bool pout =0;
    int nEvent  = 1000;
    Pythia pythia;
    ofstream myfile1, myfile2;
    
    myfile1.open("tree_time.csv");
    myfile2.open("search_time.csv");
    
    // Process selection.
  pythia.readFile("QCD.dat");
  
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowInfo = 0");
  
  pythia.init();
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    
    std::vector<Particle> particles;

    double t_tot_1 =0;
    double t_tot_2 =0;
    

    for (int i=1; i<pythia.event.size(); i++){
      if (pythia.event[i].status() > 0){
        if (abs(pythia.event[i].eta())>5) continue;
        if (!pythia.event[i].isVisible()) continue;
//  
	particles.push_back(pythia.event[i]);
	}
    }

//    vector <double> dist_t;
      int N_t= particles.size();
//    double min_dist=DBL_MAX;
//    double dist;
//    for (int i_t=0; i_t < N_t; i_t++){
//    	for (int j_t=i_t+1; j_t< N_t; j_t++){
//    	dist=del_ab(particles[i_t],particles[j_t]);
//    		dist_t.push_back(dist);
//    		if (dist< min_dist) min_dist=dist;
//    		
//    }}
//    cout << " sizes : " << N_t << "  " << dist_t.size() << "min_dist = " << min_dist << endl;
//    
//    
    
    auto t1 = std::chrono::high_resolution_clock::now();

    // Build a vector of pointers to the Particle objects.
    std::vector<const Particle*> particle_ptrs;
    for (const auto &p : particles) {
        particle_ptrs.push_back(&p);
    }

    // Build the VP-tree using our custom metric.
    VPTree tree(particle_ptrs);

    auto t2 = std::chrono::high_resolution_clock::now();

    t_tot_1= t_tot_1 + std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    myfile1 << N_t << " , " << t_tot_1 << " , ";

    const Particle* best_pair_first = nullptr;
    const Particle* best_pair_second = nullptr;
    double global_best = std::numeric_limits<double>::max();


    t1 = std::chrono::high_resolution_clock::now();
    search(tree, particle_ptrs, best_pair_first, best_pair_second, global_best);
    t2 = std::chrono::high_resolution_clock::now();

    t_tot_2= t_tot_2+ std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    myfile1 << t_tot_2 << endl;
    // Report the best pair.
    if (pout){
        if (best_pair_first && best_pair_second) {
            std::cout << "Pair with the smallest distance:\n";
            std::cout << "Particle 1: (px = " << best_pair_first->px() << ", py = " << best_pair_first->py()
                    << ", pz = " << best_pair_first->pz() << ", e = " << best_pair_first->e()
                    << ", pT = " << best_pair_first->pT() << ")\n";
            std::cout << "Particle 2: (px = " << best_pair_second->px() << ", py = " << best_pair_second->py()
                    << ", pz = " << best_pair_second->pz() << ", e = " << best_pair_second->e()
                    << ", pT = " << best_pair_second->pT() << ")\n";
            std::cout << "Minimum distance (del_ab): " << global_best << "\n";
        } else {
            std::cout << "Not enough particles to find a pair.\n";
        }
    }
  }
    return 0;
}

