/*
 *LinearPartition.h*
 header file for LinearPartition.cpp.

 author: He Zhang
 created by: 03/2019
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include <set>

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755

#define NEG_INF -2e20 
// #define testtime

using namespace std;

#ifdef lpv
  typedef float pf_type;
#else
  typedef double pf_type;
#endif

#include "logspace.h" // N.B.: uses pf_type

#ifdef lpv
  typedef int value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#endif

// A hash function used to hash a pair of any kind 
struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};



struct comp
{
    template<typename T>
    bool operator()(const T& l, const T& r) const
    {
        if (l.first == r.first)
            return l.second < r.second;
 
        return l.first < r.first;
    }
};

struct State {

    pf_type alpha;
    pf_type beta;

    State(): alpha(VALUE_MIN), beta(VALUE_MIN) {};

    inline void logplus(pf_type outside) {
        Fast_LogPlusEquals(beta, outside);
    }
  
};

// unified hyperedge
struct HEdge {
  State * left=NULL, * right=NULL; // right=null <=> Edge
  pf_type out;
  HEdge(State * l=NULL, State * r=NULL, pf_type o=VALUE_MIN) : left(l), right(r), out(o) {}; // default constructor needed
};

enum Type { // reverse topological order
  TYPE_C = 0,
  TYPE_M = 1 ,
  TYPE_M2 = 2,
  TYPE_P = 3,
  TYPE_MULTI = 4,
  TYPE_MAX
};

string type2str[5] = {"C", "M", "M2", "P", "MULTI"};

typedef unordered_map<int, State> *mypointer;
class BeamCKYParser {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    string bpp_file;
    string bpp_file_index;
    bool pf_only;
    float bpp_cutoff;
    string forest_file;
    bool mea_;
    float gamma;
    string mea_file_index;
    bool bpseq;
    bool threshknot_;
    float threshknot_threshold;
    string threshknot_file_index;
    bool is_fasta;
    bool is_lazy = false;

    // SHAPE
    bool use_shape = false;
    double m = 1.8;
    double b = -0.6;


    BeamCKYParser(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
                  string bppfile="",
                  string bppfileindex="",
                  bool pf_only=false,
                  float bpp_cutoff=0.0,
      	         string forestfile="",
                  bool mea_=false,
                  float gamma=3.0,
                  string mea_file_index="",
                  bool bpseq=false,
                  bool threshknot_=false,
                  float threshknot_threshold=0.3,
                  string threshknot_file_index="",
                  string shape_file_path="",
                  bool is_fasta=false,
                  bool is_lazy=false,
                  float threshold=9.91152);

    // DecoderResult parse(string& seq);
    void parse(string& seq);

private:
    void get_parentheses(char* result, string& seq);

    unsigned seq_length;
    string seq;

    unordered_map<int, State> *bestH, *bestP, *bestM2, *bestMulti, *bestM;
    unordered_map<int, State> bestC; // lhuang: caution: -1
    unordered_map<int, State> **best_states; // pointing to (bestC), bestM, bestM2, bestP, bestMulti

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;


    int *nucs;
    vector<int> * next_pair;
  vector<int> * prev_pair; // N.B.: lhuang: from linearsampling

  vector<int> *sortedP; // sorted by -i for P(i,j); used for M2=M+P backwards
  float inside_time;
  float deviation_threshold = 9.91152; // inside + outside > global - ..

    void prepare();
    void postprocess();

    void cal_PairProb(State& viterbi); 

    void PairProb_MEA(string & seq);

    void ThreshKnot(string & seq);

    string back_trace(const int i, const int j, const vector<vector<int> >& back_pointer);
    map<int, int> get_pairs(string & structure);

    void outside();

    void dump_forest(bool inside_only);
    void print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold);

    pf_type beam_prune(unordered_map<int, State>& beamstep);

    vector<pair<pf_type, int>> scores;

    unordered_map<pair<int,int>, pf_type, hash_pair> Pij;

    void output_to_file(string file_name, const char * type);
    void output_to_file_MEA_threshknot_bpseq(string file_name, const char * type, map<int,int> & pairs, string & seq);

    void lazyoutside();
    void backward_update(int i, int j, State & state, Type type);

    // SHAPE
    std::vector<double> SHAPE_data;

    std::vector<int> pseudo_energy_stack;

  // for lazyoutside
  float global_threshold, edge_threshold;
  int pruned = 0, local_pruned = 0, saved = 0;
  float best_inside, saved_inside;

  vector <HEdge> saved_hedges;
  HEdge best_edge; // needs default constructor
  inline void check_edge(State * left, float out, float edge_extra, State * right);
  int edge_counts[6] = {0, 0, 0, 0, 0, 0};
};

#endif //FASTCKY_BEAMCKYPAR_H
