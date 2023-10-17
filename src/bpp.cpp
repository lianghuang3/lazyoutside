/*
 *bpp.cpp*
 The main code for base pair probability calculation.

 author: He Zhang
 created by: 04/2019
*/

#include <stdio.h> 
#include <set>
#include <algorithm>
#include <stdexcept>
#include "LinearPartition.h"

using namespace std;

void BeamCKYParser::output_to_file(string file_name, const char * type) {
    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.5f\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::output_to_file_MEA_threshknot_bpseq(string file_name, const char * type, map<int,int>& pairs, string & seq) {

    int i,j;
    char nuc;
    if(!file_name.empty()) {
        printf("Outputing base pairs in bpseq format to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            fprintf(fptr, "%d %c %d\n", i, nuc, j);
        }

        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }
    else{
        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            printf("%d %c %d\n", i, nuc, j);
        }
        printf("\n");
    }

}

void BeamCKYParser::cal_PairProb(State& viterbi) {

    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            
            pf_type temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > pf_type(-9.91152)) {
                pf_type prob = Fast_Exp(temp_prob_inside);
                if(prob > pf_type(1.0)) prob = pf_type(1.0);
                if(prob < pf_type(bpp_cutoff)) continue;
                Pij[make_pair(i+1, j+1)] = prob;
            }
        }
    }

    // -o mode: output to a single file with user specified name;
    // bpp matrices for different sequences are separated with empty lines
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "a");
    } 

    // -prefix mode: output to multiple files with user specified prefix;
    else if (!bpp_file_index.empty()) {
        output_to_file(bpp_file_index, "w");
    }
    return;
}


string BeamCKYParser::back_trace(const int i, const int j, const vector<vector<int> >& back_pointer){

    if (i>j) return "";
    if (back_pointer[i][j] == -1){
        if (i == j) return ".";
        else return "." + back_trace(i+1,j, back_pointer);
    }else if (back_pointer[i][j] != 0){
        int k = back_pointer[i][j];
        assert(k + 1 > 0 && k + 1 <= seq_length);
        string temp;
        if (k == j) temp = "";
        else temp = back_trace(k+1,j, back_pointer);
        return "(" + back_trace(i+1,k-1, back_pointer) + ")" + temp;
    }
    assert(false);
    return "";
}

map<int, int> BeamCKYParser::get_pairs(string & structure){
    map<int, int> pairs;
    stack<int> s;
    int index = 1;
    int pre_index = 0;
    for (auto & elem : structure){
        if (elem == '(') s.push(index);
        else if(elem == ')'){
            pre_index = s.top();
            pairs[pre_index] = index;
            pairs[index] = pre_index;
            s.pop();
        }
        index++;
    }
    return pairs;
}

void BeamCKYParser::ThreshKnot(string & seq){
    
    map<int, pf_type> rowprob;
    vector<tuple<int, int, pf_type> > prob_list;

    map<int, int> pairs;
    set<int> visited;

    for(auto& pij : Pij){
        auto i = pij.first.first; //index starts from 1
        auto j = pij.first.second; 
        auto score = pij.second;

        if (score < threshknot_threshold) continue;

        prob_list.push_back(make_tuple(i,j,score));

        rowprob[i] = max(rowprob[i], score);
        rowprob[j] = max(rowprob[j], score);

    }

    for(auto& elem : prob_list){

        auto i = std::get<0>(elem);
        auto j = std::get<1>(elem);
        auto score =  std::get<2>(elem);

        if (score == rowprob[i] && score == rowprob[j]){

            if ((visited.find(i) != visited.end()) || (visited.find(j) != visited.end())) continue;
            visited.insert(i);
            visited.insert(j);

            pairs[i] = j;
            pairs[j] = i;
        }
    }

    fprintf(stderr, "%s\n", seq.c_str());
    output_to_file_MEA_threshknot_bpseq(threshknot_file_index, "w", pairs, seq);
}

void BeamCKYParser::PairProb_MEA(string & seq) {
    
    vector<vector<pf_type> > OPT;
    OPT.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) OPT[i].resize(seq_length);

    vector<vector<pf_type>> P;
    P.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) P[i].resize(seq_length);

    vector<vector<int> > back_pointer;
    back_pointer.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) back_pointer[i].resize(seq_length);

    vector<vector<int>> paired;
    paired.resize(seq_length);

    vector<pf_type> Q;
    for (int i = 0; i < seq_length; ++i) Q.push_back(pf_type(1.0));

    for(auto& pij : Pij){
        auto i = pij.first.first-1;
        auto j = pij.first.second-1;
        auto score = pij.second;

        P[i][j] = score;

        paired[i].push_back(j);
        Q[i] -= score;
        Q[j] -= score;
    }

    for (int i = 0; i < seq_length; ++i) std::sort (paired[i].begin(), paired[i].end());
    for (int l = 0; l< seq_length; l++){
        for (int i = 0; i<seq_length - l; i++){
            int j = i + l;
            if (i == j){
                OPT[i][j] = Q[i];
                back_pointer[i][j] = -1;
                continue;
            }
            OPT[i][j] = OPT[i][i] + OPT[i+1][j];
            back_pointer[i][j] = -1;
            for (int k : paired[i]){
                if (k>j) break;
                pf_type temp_OPT_k1_j;
                if (k<j) temp_OPT_k1_j = OPT[k+1][j];
                else temp_OPT_k1_j = pf_type(0.);
                auto temp_score = 2 * gamma * P[i][k] + OPT[i+1][k-1] + temp_OPT_k1_j;
                if (OPT[i][j] < temp_score){
                    OPT[i][j] = temp_score;
                    back_pointer[i][j] = k;
                }
            }
        }
    }

    auto structure = back_trace(0,seq_length-1, back_pointer);

    if (!bpseq){
        if(!mea_file_index.empty()) {
            FILE *fptr = fopen(mea_file_index.c_str(), "w"); 
            if (fptr == NULL) { 
                printf("Could not open file!\n"); 
                return; 
            }
            fprintf(fptr, "%s\n", seq.c_str());
            fprintf(fptr, "%s\n\n", structure.c_str());
        }

        else{
            printf("%s\n", seq.c_str());
            printf("%s\n\n", structure.c_str());
        }
    }

    else{
        auto pairs = get_pairs(structure);
        output_to_file_MEA_threshknot_bpseq(mea_file_index, "w", pairs, seq);
    }
}


void BeamCKYParser::outside() {
      
    struct timeval bpp_starttime, bpp_endtime;
    gettimeofday(&bpp_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
            Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
                    
#else
            newscore = score_external_unpaired(j+1, j+1);
            Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + newscore);
#endif
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + newscore);
#endif
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lpv
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
#else
                            newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + newscore);
#endif
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
#ifndef lpv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lpv
                                newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);
                                // SHAPE for Vienna only
                                if (use_shape)
                                {
                                    newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                }


                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            } else {
                                // single branch
#ifdef lpv
                                newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed + 
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lpv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore/kT);
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = newscore/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = newscore;
#endif
                    pf_type m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore/kT;

#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore;
#endif
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
#ifdef lpv
                        newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + newscore/kT));
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepC.beta + newscore);
#endif
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lpv
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore/kT);
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore);
#endif
                }
            }
        }
    }  // end of for-loo j

    gettimeofday(&bpp_endtime, NULL);
    double bpp_elapsed_time = bpp_endtime.tv_sec - bpp_starttime.tv_sec + (bpp_endtime.tv_usec-bpp_starttime.tv_usec)/1000000.0;

    if(is_verbose) fprintf(stderr,"Base Pairing Probabilities Calculation Time: %.2f seconds.\n", bpp_elapsed_time);

    fflush(stdout);

    return;
}

void BeamCKYParser::lazyoutside() {
    struct timeval parse_starttime, parse_endtime;

    printf("lazy\n");
    gettimeofday(&parse_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;
    bestC[-1].alpha = 0.0; // lhuang: N.B.: was 50001... BAD! hash not array!
    global_threshold = bestC[seq_length-1].alpha - deviation_threshold; //- 9.91152;
    int tot = 0, vis = 0, weird = 0;
    pruned = 0;
    saved = 0;
    memset(edge_counts, 0, sizeof(edge_counts));
    
    // from right to left
    for(int j = seq_length-1; j > 0; --j) {
      if (bestC[j].beta > -9.91152) {
	backward_update(0, j, bestC[j], TYPE_C);
	vis ++;
      }
      tot++;
      
      for (int type = TYPE_M; type < TYPE_MAX; type++) { // reverse topol order: C->M->M2->P->Multi        
        for (auto & item : best_states[type][j]) {
          int i = item.first;
          State & state = item.second;
          if (state.beta > -9.91152) {
            //if (state.alpha + state.beta > global_threshold) {
              backward_update(i, j, state, (Type)type);
              vis ++;
            }
          tot ++;
        }
      }
    }

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

// #ifdef testtime
    if(is_verbose) {
      float tot_edges = pruned + saved;
      printf("Lazy Outside Time: %.2f seconds (%.1f%%) (visited edges: pruned %d, saved %d) (weird %d).\n",
	     parse_elapsed_time, 100. * parse_elapsed_time / inside_time,
	     pruned, saved, weird);
      printf("nodes: visited %d; total %d (%.1f%%)\n", vis, tot, vis * 100. / tot);    
      printf("edge counts: ");
      for (int t = 0; t < 5; t++)
	printf("%s: %d (%.1f%%)  ", type2str[t].c_str(), edge_counts[t], (edge_counts[t] * 100. / tot_edges));
      printf("\n");
    }
// #endif

    fflush(stdout);
    }

inline void BeamCKYParser::check_edge(State * left, float out, float edge_extra=0, State * right=NULL) {

  float edge_inside = left->alpha + edge_extra + (right? right->alpha : 0);  
  float edge_outside = out + edge_extra;
  
  if (edge_inside > edge_threshold) {
    Fast_LogPlusEquals(saved_inside, edge_inside); 
    saved_hedges.push_back(HEdge(left, right, edge_outside)); 
  }
  else {
    local_pruned ++;
    if (saved_hedges.empty() && edge_inside > best_inside) {
      best_inside = edge_inside;
      best_edge = HEdge(left, right, edge_outside);
    }
  }
}

void BeamCKYParser::backward_update(int i, int j, State & state, Type type) {
  //printf("back %d %d %d\n", i, j, type);
  float out = state.beta;

  float global_total = bestC[seq_length-1].alpha; // (log)Z

  if (state.alpha + state.beta > global_total + 0.01) {
    printf("WARNING %s[%d, %d]: %.2f %.2f dev %.2f\n",
	   type2str[type].c_str(), i, j, state.alpha, state.beta,
	   state.alpha + state.beta - global_total);
  }

  // the following are (global); not ideal
  edge_threshold = global_threshold - out; // edge_inside + out > global - threshold
  saved_inside = VALUE_MIN;
  best_inside = VALUE_MIN;
  local_pruned = 0;
  saved_hedges.clear();
  // best_edge
  
  switch (type) {
  case TYPE_C: {
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
  
    // C = C + U
    if (j == 0) { //samplestate.append(alphalist, 0.0, MANNER_C_eq_C_plus_U); // hzhang: N.B. j == 0
    }

    else{
      // C = C + U
      check_edge(&bestC[j-1], out);      
      //edge_counts[type] ++;
	
      // C = C + P
      for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
        int i = item.first;
        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

        int k = i - 1;
        State& Pstate = item.second, & Cstate = bestC[k];

        int nuck = nuci_1;
        int nuck1 = nuci;
        float score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                              nucj, nucj1, seq_length) /kT;
	check_edge(&Cstate, out, score_external_paired, &Pstate);
	//edge_counts[type] ++;
      }
    }
  }
  break;
  case TYPE_P: {
    saved_inside = VALUE_MIN;

    check_edge(&bestH[j][i], out); // hidden H=>P hyperedge; special NULLary edge, but still H
    //edge_counts[type] ++;
              
    int newscore;
    int nuci = nucs[i];
    int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
    int nucj = nucs[j];
    int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

    int p, q, nucq, nucq1, nucp, nucp_1;
    // helix or single_branch
    for (q = j - 1; q >= std::max(j - SINGLE_MAX_LEN, i+5); --q) { // no sharp turn
      nucq = nucs[q];
      nucq1 = nucs[q + 1];
      p = next_pair[nucq][i]; 
      while (p != -1 && p <= q - 4 && ((p - i) + (j - q) - 2 <= SINGLE_MAX_LEN)) {
        auto iterator = bestP[q].find(p);
        //if(bestP[q].find (p) != bestP[q].end()) {
	if (iterator != bestP[q].end()) {
          nucp = nucs[p];
          nucp_1 = nucs[p - 1];
          
          float score_single = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
					       nucp_1, nucp, nucq, nucq1) / kT; // same for vienna

	  //check_edge(&bestP[q][p], out, score_single);
	  check_edge(&(iterator->second), out, score_single);	  
	  //edge_counts[type] ++;
        }
        p = next_pair[nucq][p];
      }
    }
    // hairpin
    // Multiloop
    auto iterator = bestMulti[j].find (i);
    //if(bestMulti[j].find(i) != bestMulti[j].end()) {
    if (iterator != bestMulti[j].end()) {
        float score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length)/kT;
	//check_edge(&bestMulti[j][i], out, score_multi);
	check_edge(&(iterator->second), out, score_multi);
	//edge_counts[type] ++;
    }
  }
  break;
  case TYPE_M: {
    //printf("M %d %d %d %.2f %.2f\n", i, j, type, state.alpha, state.beta);
    int nuci = nucs[i];
    int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
    float edge_inside;
    
    // M = M + U
    //if(j > i+1 && bestM[j-1].find(i) != bestM[j-1].end()) {
    if (j > i+1) {
      auto iterator = bestM[j-1].find(i);
      if (iterator != bestM[j-1].end()) {
	//check_edge(&bestM[j-1][i], out);
	check_edge(&(iterator->second), out);
	//edge_counts[type] ++;
      }
    }

    auto iterator = bestP[j].find(i);
    //if(bestP[j].find(i) != bestP[j].end()) {
    if (iterator != bestP[j].end()) {
      // M = P
        float M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)/kT;
	//check_edge(&bestP[j][i], out, M1_score);
	check_edge(&(iterator->second), out, M1_score);
	//edge_counts[type] ++;
    }

    iterator = bestM2[j].find(i);
    //if(bestM2[j].find(i) != bestM2[j].end()) {
    if (iterator != bestM2[j].end()) {
      // M = M2
      //check_edge(&bestM2[j][i], out);
      check_edge(&(iterator->second), out);
      //edge_counts[type] ++;
    }
  }
  break;
  case TYPE_M2: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

    if(sortedP[j].size() == 0){
      // lhuang: M2 might be short, and b might be big
      // so only explore P that fits M2 = M + P
      for(auto const& item : bestP[j])
        sortedP[j].push_back(-item.first);
      sort(sortedP[j].begin(), sortedP[j].end());
    }
    // M2 = M + P
    for (auto & item : sortedP[j]) { // map not unorderd_map
      int k = -item;
      if (k > i + 4) { // lhuang: +4
          int m = k - 1;
          auto iterator = bestM[m].find(i);
          //if(bestM[m].find(i) != bestM[m].end()) {
	  if (iterator != bestM[m].end()) {
              int nuck = nucs[k];
              int nuck_1 = (k-1>-1) ? nucs[k-1] : -1;
              // M2(i,j) = M(i,k-1) + P(k,j)
              float M1_score = - v_score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length)/kT;
              State & Mstate = iterator->second; // bestM[m][i], & Pstate = bestP[j][k];
	      State & Pstate = bestP[j][k];
	      check_edge(&Mstate, out, M1_score, &Pstate);
	      //edge_counts[type] ++;		    
          }
      }
      else break;
    }
  }
  break;
  case TYPE_MULTI: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int jprev = prev_pair[nuci][j];
    //if (jprev > i+10 && bestMulti[jprev].find(i) != bestMulti[jprev].end()) { // no sharp turn
    if (jprev > i+10) {
      auto iterator = bestMulti[jprev].find (i);
      if (iterator != bestMulti[jprev].end()) { // no sharp turn
	// Multi = Multi + jump
	//check_edge(&bestMulti[jprev][i], out);
	check_edge(&(iterator->second), out);
	//edge_counts[type] ++;
      }
    }

    for (int q = j - 1; q >= jprev; q--) { 
      for (int p = i+1; p <= q - 9 && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN; p++){
        auto iterator = bestM2[q].find(p);
        //if(bestM2[q].find(p) != bestM2[q].end()){
	if (iterator != bestM2[q].end()){
	  //check_edge(&bestM2[q][p], out);
	  check_edge(&(iterator->second), out);
	  //edge_counts[type] ++;
       }  
     }
    } 
  }
    break;
  case TYPE_MAX:
    throw std::invalid_argument("invalid state type");
  }

  pruned += local_pruned; //global
  float delta; // scaling factor to compensate for edge pruning
  if (!saved_hedges.empty()) 
    delta = state.alpha - saved_inside;
  else { // all pruned: use best hyperedge
    delta = state.alpha - best_inside;
    saved_hedges.push_back(best_edge);
    pruned --; // one more edge recovered
  }

  for (auto & edge : saved_hedges) {
    State *left = edge.left, *right = edge.right;
    if (!right) // edge
      left->logplus(edge.out + delta);
    else { // hyperedge
      left->logplus(edge.out + delta + right->alpha);
      right->logplus(edge.out + delta + left->alpha);      
    }      
  }
  saved += saved_hedges.size();
}
