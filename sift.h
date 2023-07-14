#ifndef SIFT_H
#define SIFT_H

#include "Pythia8/Pythia.h"
#include <float.h>
using namespace Pythia8;

bool sortFunction(const pair<Particle, vector<int>> &a, const pair<Particle, vector<int>> &b){
	return a.first.pT() > b.first.pT();
}
// this function takes the index from the vector delta, and gives the indices of the corresponding pair
// in vector part_temp.
array<int,2> ij_calc (const int& k, const int& N){
	int i = N - 2 - floor(sqrt(-8 * k + 4 * N * (N - 1) - 7) / 2 - 0.5);
	int j = k + i + 1 - N*i + (i*(i+1)) / 2;
	return {i,j};
}
// inverse of ij_calc
int k_calc (const int& i, const int& j, const int& N){
	return ((N-1) * i - (i * (i + 1) / 2) + j-1);
}
//gives the index of the minimum value of a vector.
int min_index(const vector<double>& vec, const int& N_delta){
  	int index=0;
  	for (int i=0; i< N_delta; i++){
    		if (vec[i] < vec[index]) index = i;
  	}
  	return index;
}

// this function replaces the ith particle with combined particle
void join (vector<pair<Particle,vector<int>>>& v, int i, int j){
  	Vec4 p_tot=v[i].first.p() + v[j].first.p();
  	Particle p_new;
  	p_new.p(p_tot);
  	p_new.m(p_tot.mCalc());
  	vector<int> a= v[i].second;
	vector<int> b= v[j].second;
	a.insert(a.end(), b.begin(), b.end());
	v[i]=make_pair(p_new,a);
}  
double dm_ab2 (const Particle& a1, const Particle& b1){
	return 2.0*(a1.p()*b1.p());
}
double eps_ab (const Particle& a2, const Particle& b2){
	double eT_a2= pow2(a2.pT())+pow2(a2.m());
  	double eT_b2= pow2(b2.pT())+pow2(b2.m());
	return sqrt(eT_a2 * eT_b2)/(eT_a2 + eT_b2);
}  
double dR_t_ab2 (const Particle& a2, const Particle& b2){
	double eT_a2= pow2(a2.pT())+pow2(a2.m());
	double eT_b2= pow2(b2.pT())+pow2(b2.m());
	return dm_ab2(a2,b2)/sqrt(eT_a2 * eT_b2);
}  

double del_ab (const Particle& a, const Particle& b){
	double eT_a2= pow2(a.pT())+pow2(a.m());
	double eT_b2= pow2(b.pT())+pow2(b.m());
	return dm_ab2(a,b)/(eT_a2 + eT_b2);
} 

vector<pair<Particle,vector<int>>> sift (const Event& event){
	vector<pair<Particle,vector<int>>> part_temp;
	vector<pair<Particle,vector<int>>> jet_list;
	for (int i=1; i<event.size(); i++){
      	if (event[i].status() > 0){
      		if (abs(event[i].eta())>5) continue;
      		if (!event[i].isVisible()) continue;
        		vector<int> v1={i};
      		part_temp.push_back(make_pair(event[i],v1));
      	}
	}
	int N=part_temp.size();
	vector<double> delta(N*(N-1)/2);
	int k=0;
	for (int i1=0; i1< N-1; i1++){
      	for (int j1=i1+1; j1< N; j1++){
      		delta[k]=del_ab(part_temp[i1].first , part_temp[j1].first);
		      k++;
      	}
    	}
    	
    	while (any_of(delta.begin(),delta.end(), [](double element){
    		return element < DBL_MAX;
    		}))
   	{
   		int k_min=min_index(delta,N*(N-1)/2);
   		int i_min=ij_calc(k_min,N)[0];
   		int j_min=ij_calc(k_min,N)[1];
   		double dR_t2= dR_t_ab2(part_temp[i_min].first,part_temp[j_min].first); 
   		double eps = eps_ab(part_temp[i_min].first,part_temp[j_min].first);
   		int k_test;
      	if (eps > (dR_t2/4.0)){
        		join (part_temp, i_min, j_min);
		      k_test=0;
			// this loop is set just to take into account the case when only one particle will be remaining,
			// in that case all values of delta will be DBL_MAX, but we don't want to miss the last jet.
			for (double test: delta){
        			if (test != DBL_MAX) k_test++;
        			if (k_test >1) break;
			}
			if (k_test==1){
				jet_list.push_back(part_temp[i_min]);
				fill(delta.begin(), delta.end(), DBL_MAX);
				continue;
			}
			for (int i2=0; i2< i_min; i2++){
				if (delta[k_calc(i2,i_min,N)] == DBL_MAX) continue;
				delta[k_calc(i2,i_min,N)]= del_ab(part_temp[i2].first , part_temp[i_min].first);
			}
			for (int j2=i_min+1; j2<N; j2++){
				if (delta[k_calc(i_min,j2,N)] == DBL_MAX) continue;
				delta[k_calc(i_min,j2,N)]= del_ab(part_temp[i_min].first, part_temp[j2].first);
			}
			for (int i3=0; i3< j_min; i3++){
				delta[k_calc(i3,j_min,N)]= DBL_MAX;
			}
			for (int j3=j_min+1; j3<N; j3++){
				delta[k_calc(j_min,j3,N)]= DBL_MAX;
			}
		}
		else if ( eps > (1.0/dR_t2)){
			k_test=0;
			for (double test: delta){
				if (test != DBL_MAX) k_test++;
				if (k_test >3) break;
			}
			if (k_test ==3){
				set<int> temp_set;
				for (int i8=0; i8< N*(N-1)/2; i8++){
					if (delta[i8] != DBL_MAX){
						temp_set.insert(ij_calc(i8,N)[0]);
						temp_set.insert(ij_calc(i8,N)[1]);
					}
				}
				for (int i9: temp_set){
					jet_list.push_back(part_temp[i9]);
				}
				fill(delta.begin(), delta.end(), DBL_MAX);
				continue;
      		} 
			jet_list.push_back(part_temp[i_min]);
			jet_list.push_back(part_temp[j_min]);        
			for (int i4=0; i4< i_min; i4++){
				delta[k_calc(i4,i_min,N)]= DBL_MAX;
			}
     			for (int j4=i_min+1; j4<N; j4++){
     				delta[k_calc(i_min,j4,N)]= DBL_MAX;
     			}
			for (int i5=0; i5< j_min; i5++){
				delta[k_calc(i5,j_min,N)]= DBL_MAX;
			}
     			for (int j5=j_min+1; j5<N; j5++){
				delta[k_calc(j_min,j5,N)]= DBL_MAX;
			}
		}
     		else {
			k_test=0;
        		for (double test: delta){
        			if (test != DBL_MAX) k_test++;
				if (k_test >1) break;
			}
       
			if (part_temp[i_min].first.pT() < part_temp[j_min].first.pT()){
				if (k_test==1){
					jet_list.push_back(part_temp[j_min]);
     					fill(delta.begin(), delta.end(), DBL_MAX);
					continue;
				}
				for (int i6=0; i6< i_min; i6++){
					delta[k_calc(i6,i_min,N)]= DBL_MAX;
				}
				for (int j6=i_min+1; j6<N; j6++){
					delta[k_calc(i_min,j6,N)]= DBL_MAX;
				}
         
       		} 
       		else {
         			if (k_test==1){
					jet_list.push_back(part_temp[i_min]);
					fill(delta.begin(), delta.end(), DBL_MAX);
					continue;
				}
				for (int i7=0; i7< j_min; i7++){
					delta[k_calc(i7,j_min,N)]= DBL_MAX;
				}
				for (int j7=j_min+1; j7<N; j7++){
					delta[k_calc(j_min,j7,N)]= DBL_MAX;
				}
			}
		}
	}
	return jet_list;
}

// Defining a separate funtion to get substructure.
// this is similar to sift function, but without isolation condition, and to get fixed number (currently 2) of output jets.
vector<pair<Particle,vector<int>>> sub_sift (const Event& event){
	vector<pair<Particle,vector<int>>> part_temp;
	vector<pair<Particle,vector<int>>> subjet_list;
	for (int i=0; i<event.size(); i++){
      	if (event[i].status() > 0){
      		if (abs(event[i].eta())>5) continue;
      		if (!event[i].isVisible()) continue;
        		vector<int> v1={i};
      		part_temp.push_back(make_pair(event[i],v1));
      	}
	}
	int N=part_temp.size();
	vector<double> delta(N*(N-1)/2);
	int k=0;
	for (int i1=0; i1< N-1; i1++){
      	for (int j1=i1+1; j1< N; j1++){
      		delta[k]=del_ab(part_temp[i1].first , part_temp[j1].first);
		      k++;
      	}
    	}
    	
    	while (any_of(delta.begin(),delta.end(), [](double element){
    		return element < DBL_MAX;
    		}))
   	{
   		int k_min=min_index(delta,N*(N-1)/2);
   		int i_min=ij_calc(k_min,N)[0];
   		int j_min=ij_calc(k_min,N)[1];
   		double dR_t2= dR_t_ab2(part_temp[i_min].first,part_temp[j_min].first); 
   		double eps = eps_ab(part_temp[i_min].first,part_temp[j_min].first);
   		int k_test;
   		k_test=0;
		for (double test: delta){
			if (test != DBL_MAX) k_test++;
			if (k_test >1) break;
		}
		if (k_test==1){
			set<int> temp_set;
			for (int i8=0; i8< N*(N-1)/2; i8++){
				if (delta[i8] != DBL_MAX){
					temp_set.insert(ij_calc(i8,N)[0]);
					temp_set.insert(ij_calc(i8,N)[1]);
				}
			}
			for (int i9: temp_set){
				subjet_list.push_back(part_temp[i9]);
			}
			break;
		}
			
      	if (eps > (dR_t2/4.0) or eps > (1.0/dR_t2)){
        		join (part_temp, i_min, j_min);
		      for (int i2=0; i2< i_min; i2++){
				if (delta[k_calc(i2,i_min,N)] == DBL_MAX) continue;
				delta[k_calc(i2,i_min,N)]= del_ab(part_temp[i2].first , part_temp[i_min].first);
			}
			for (int j2=i_min+1; j2<N; j2++){
				if (delta[k_calc(i_min,j2,N)] == DBL_MAX) continue;
				delta[k_calc(i_min,j2,N)]= del_ab(part_temp[i_min].first, part_temp[j2].first);
			}
			for (int i3=0; i3< j_min; i3++){
				delta[k_calc(i3,j_min,N)]= DBL_MAX;
			}
			for (int j3=j_min+1; j3<N; j3++){
				delta[k_calc(j_min,j3,N)]= DBL_MAX;
			}
		}
     		else {
			k_test=0;
        		
			if (part_temp[i_min].first.pT() < part_temp[j_min].first.pT()){
				for (int i6=0; i6< i_min; i6++){
					delta[k_calc(i6,i_min,N)]= DBL_MAX;
				}
				for (int j6=i_min+1; j6<N; j6++){
					delta[k_calc(i_min,j6,N)]= DBL_MAX;
				}
         
       		} 
       		else {
         			for (int i7=0; i7< j_min; i7++){
					delta[k_calc(i7,j_min,N)]= DBL_MAX;
				}
				for (int j7=j_min+1; j7<N; j7++){
					delta[k_calc(j_min,j7,N)]= DBL_MAX;
				}
			}
		}
	}
	return subjet_list;
}

#endif    
