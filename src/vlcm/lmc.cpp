//Copyright <2016> <Chu-Min Li & Hua Jiang>
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>

#include "lmc.h"

#define NONE -1
#define DELIMITER 0
#define PASSIVE false
#define ACTIVE true
#define P_TRUE 2
#define P_FALSE 0
#define NO_REASON -3
#define MAX_NODE 1000000 // TODO originally 80000000
#define max_expand_depth 100000
#define push(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer

#define CUR_CLQ_SIZE Clique_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]

#define SET_EDGE(row, col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) |= (1 << ((col) & 7)))

#define iMatrix(i) (Adj_Matrix+(i)*MATRIX_ROW_WIDTH)
#define Matrix(i, j) ((*((i) + ((j) >> 3))) & (1 << ((j) & 7)))

#define CORE_NO Vertex_UB


class LMC_Instance {

public:
    LMC_Instance();

    ~LMC_Instance();

    int NB_NODE, ADDED_NODE, NB_EDGE,
            MAX_CLQ_SIZE, MAX_ISET_SIZE, INIT_CLQ_SIZE, MATRIX_ROW_WIDTH, MAX_VERTEX_NO, K_CORE_G;
    int Max_Degree;
    int Node_Degree[MAX_NODE];
    char Node_State[MAX_NODE];
    int **Node_Neibors;

    int Candidate_Stack_fill_pointer;
    int Candidate_Stack[MAX_NODE * 2];
    int Vertex_UB[MAX_NODE * 2];
    int Clique_Stack_fill_pointer;
    int *Clique_Stack, *MaxCLQ_Stack;
    int Cursor_Stack[max_expand_depth];
    int Cursor_Stack_fill_pointer;

    int *Node_Reason;
    unsigned char *Adj_Matrix;

    int iSET_COUNT;
    int *iSET_Size;
    bool *iSET_State;
    bool *iSET_Used;
    bool *iSET_Tested;
    int *iSET_Index;
    bool *iSET_Involved;
    bool *Is_Tested;
    int **iSET;

    int *REASON_STACK;
    int REASON_STACK_fill_pointer;
    int *CONFLICT_ISET_STACK;
    int CONFLICT_ISET_STACK_fill_pointer;
    int *ADDED_NODE_iSET;
    int *REDUCED_iSET_STACK;
    int REDUCED_iSET_STACK_fill_pointer;
    int *PASSIVE_iSET_STACK;
    int PASSIVE_iSET_STACK_fill_pointer;
    int *FIXED_NODE_STACK;
    int FIXED_NODE_STACK_fill_pointer;
    int *UNIT_STACK;
    int UNIT_STACK_fill_pointer;
    int *NEW_UNIT_STACK;
    int NEW_UNIT_STACK_fill_pointer;

    int Rollback_Point;
    int Branching_Point;

    int *Old_Name;
    int *Second_Name;
    int NB_CANDIDATE, FIRST_INDEX;
    int START_MAXSAT_THD;

    int Extra_Node_Stack[100000];

    int LAST_IN;
    bool REBUILD_MATRIX;
    int CUR_MAX_NODE;
    int *Init_Adj_List;
    int BLOCK_COUNT;
    int *BLOCK_LIST[100];
    int *Adj_List;


    bool is_adjacent(int node1, int node2) const;

    void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge, int offset);

    void import_graph_instance(const std::vector<std::vector<int>> &adjacencyList);

    bool sort_by_degeneracy_ordering();

    bool re_number_adj(int node) const;

    bool re_number(int node) const;

    bool addIntoIsetTomitaBis_adj(int node);

    bool addIntoIsetTomitaBis(int node);

    bool cut_by_iset_less_vertices();

    bool cut_by_iset_last_renumber();

    int fix_newNode_for_iset(int fix_node, int fix_iset);

    int fix_oldNode_for_iset(int fix_node, int fix_iset);

    int fix_node_iset(int fix_iset);

    int unit_iset_process();

    int unit_iset_process_used_first();

    void identify_conflict_sets(int iset_idx);

    void enlarge_conflict_sets();

    void rollback_context_for_maxsatz(int start_fixed, int start_passive, int start_reduced);

    void reset_context_for_maxsatz();

    int further_test_reduced_iset(int start);

    int fix_anyNode_for_iset(int fix_node, int fix_iset);

    bool inc_maxsatz_lookahead_by_fl2();

    int inc_maxsatz_on_last_iset();

    int open_new_iset_old(int i);

    int simple_further_test_node(int start);

    int test_node_for_failed_nodes(int node, int iset);

    bool test_by_eliminate_failed_nodes();

    bool cut_by_inc_maxsat_eliminate_first();

    void compute_subgraph_degree(int start);

    void allocate_memory_for_maxsat();

    void store_maximum_clique(int node);

    void store_maximum_clique2();

    bool reduce_subgraph(int start);

    void rebuild_matrix(int start);

    bool cut_by_inc_ub();

    bool find_3_clique(int node);

    void init_for_search();

    void allocate_memory();

    void search_maxclique();

    void build_init_matrix();

    bool search_in_2_k_core_graph();

    void free_block();

    void reduce_instance();

    bool initialize();

    std::vector<int> getMaxCliqueVector() const;
};


LMC_Instance::LMC_Instance() {
    K_CORE_G = 0;
    Max_Degree = 0;
    Candidate_Stack_fill_pointer = 0;
    Cursor_Stack_fill_pointer = 0;
    iSET_COUNT = 0;
    REASON_STACK_fill_pointer = 0;
    REDUCED_iSET_STACK = Node_Degree;
    REDUCED_iSET_STACK_fill_pointer = 0;
    PASSIVE_iSET_STACK_fill_pointer = 0;
    FIXED_NODE_STACK_fill_pointer = 0;
    UNIT_STACK_fill_pointer = 0;
    NEW_UNIT_STACK_fill_pointer = 0;
    NB_CANDIDATE = 0;
    START_MAXSAT_THD = 15;
    REBUILD_MATRIX = false;
    BLOCK_COUNT = 0;
    MAX_VERTEX_NO = 0;

    Node_Reason = nullptr;
    ADDED_NODE_iSET = nullptr;
    FIXED_NODE_STACK = nullptr;
    iSET_State = nullptr;
    iSET_Used = nullptr;
    iSET_Tested = nullptr;
    UNIT_STACK = nullptr;
    NEW_UNIT_STACK = nullptr;
    PASSIVE_iSET_STACK = nullptr;
    iSET_Involved = nullptr;
    CONFLICT_ISET_STACK = nullptr;
    REASON_STACK = nullptr;
    Is_Tested = nullptr;
    Adj_Matrix = nullptr;
    Old_Name = nullptr;
    iSET_Index = nullptr;
    Second_Name = nullptr;
    iSET_Size = nullptr;
    MaxCLQ_Stack = nullptr;
    Clique_Stack = nullptr;
    Node_Neibors = nullptr;
    Adj_List = nullptr;
    iSET = nullptr;
}

LMC_Instance::~LMC_Instance() {
    free(Node_Reason);
    free(ADDED_NODE_iSET);
    free(FIXED_NODE_STACK);
    free(iSET_State);
    free(iSET_Used);
    free(iSET_Tested);
    free(UNIT_STACK);
    free(NEW_UNIT_STACK);
    free(PASSIVE_iSET_STACK);
    free(iSET_Involved);
    free(CONFLICT_ISET_STACK);
    free(REASON_STACK);
    free(Is_Tested);
    free(Adj_Matrix);
    free(Old_Name);
    free(iSET_Index);
    free(Second_Name);
    free(iSET_Size);
    free(MaxCLQ_Stack);
    free(Clique_Stack);
    free(Node_Neibors);
    free(Adj_List);
    if (iSET != nullptr)
        free(iSET[0]);
    free(iSET);
}

static int int_cmp(const void *a, const void *b) {
    return *((int *) b) - *((int *) a);
}

bool LMC_Instance::is_adjacent(int node1, int node2) const {
    int neibor, *neibors;
    neibors = Node_Neibors[node1];
    for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
        if (neibor == node2) {
            return true;
        }
    }
    return false;
}

void LMC_Instance::allcoate_memory_for_adjacency_list(int nb_node, int nb_edge, int offset) {
    int i, block_size = 40960000, free_size = 0;
    Init_Adj_List = (int *) malloc((2 * nb_edge + nb_node) * sizeof(int));
    if (Init_Adj_List == nullptr) {
        for (i = 1; i <= NB_NODE; i++) {
            if (Node_Degree[i - offset] + 1 > free_size) {
                Node_Neibors[i] = (int *) malloc(block_size * sizeof(int));
                BLOCK_LIST[BLOCK_COUNT++] = Node_Neibors[i];
                free_size = block_size - (Node_Degree[i - offset] + 1);
            } else {
                Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1 - offset] + 1;
                free_size = free_size - (Node_Degree[i - offset] + 1);
            }
        }
    } else {
        BLOCK_COUNT = 1;
        BLOCK_LIST[0] = Init_Adj_List;
        Node_Neibors[1] = Init_Adj_List;
        for (i = 2; i <= NB_NODE; i++) {
            Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1 - offset] + 1;
        }
    }
}

void LMC_Instance::import_graph_instance(const std::vector<std::vector<int>> &adjacencyList) {
    // node ids must be >= 1, thus we add 1 to every node (is subtracted at the end)
    NB_NODE = static_cast<int>(adjacencyList.size());  // number of nodes
    // get number of edges and node degrees
    int nb_edge = 0;
    for (int i = 1; i <= NB_NODE; i++) {
        int degree = static_cast<int>(adjacencyList[i - 1].size());
        nb_edge += degree;
        Node_Degree[i] = degree;
    }
    NB_EDGE = nb_edge / 2;
    REDUCED_iSET_STACK = Node_Degree;
    Node_Neibors = (int **) malloc((NB_NODE + 1) * sizeof(int *));
    allcoate_memory_for_adjacency_list(NB_NODE, NB_EDGE, 0);
    // copy adjacency list
    for (int i = 1; i <= NB_NODE; i++) {
        for (int j = 0; j < Node_Degree[i]; j++) {
            Node_Neibors[i][j] = adjacencyList[i - 1][j] + 1;
        }
    }

    Max_Degree = 0;
    for (int node = 1; node <= NB_NODE; node++) {
        Node_Neibors[node][Node_Degree[node]] = NONE;

        if (Node_Degree[node] > Max_Degree)
            Max_Degree = Node_Degree[node];
    }
}

bool LMC_Instance::sort_by_degeneracy_ordering() {
    int *degree_counter, *where;
    int p1, i, node = NONE, neibor, *neibors, t, j, h, k;
    int cur_degree;
    INIT_CLQ_SIZE = 0;
    //printf("I computing initial degeneracy ordering...\n");

    where = Candidate_Stack + NB_NODE + 1;
    degree_counter = Vertex_UB + NB_NODE + 1;
    memset(degree_counter, 0, (Max_Degree + 1) * sizeof(int));

    for (node = 1; node <= NB_NODE; node++) {
        degree_counter[Node_Degree[node]]++;
    }
    j = 0;
    for (i = 0; i <= Max_Degree; i++) {
        k = degree_counter[i];
        degree_counter[i] = j;
        j += k;
    }

    for (node = 1; node <= NB_NODE; node++) {
        Candidate_Stack[t = degree_counter[Node_Degree[node]]++] = node;
        where[node] = t;
    }
    for (i = Max_Degree; i > 0; i--) {
        degree_counter[i] = degree_counter[i - 1];
    }
    degree_counter[0] = 0;

    Candidate_Stack[NB_NODE] = DELIMITER;
    ptr(Candidate_Stack) = NB_NODE + 1;

    p1 = 0;
    cur_degree = Node_Degree[Candidate_Stack[p1]];
    while (p1 < NB_NODE) {
        node = Candidate_Stack[p1];

        if (Node_Degree[node] > cur_degree)
            cur_degree = Node_Degree[node];
        CORE_NO[p1] = cur_degree;

        if (cur_degree > K_CORE_G)
            K_CORE_G = cur_degree;

        if (p1 < NB_NODE - 1 && Node_Degree[node] == Node_Degree[Candidate_Stack[p1 + 1]]) {
            degree_counter[Node_Degree[node]] = p1 + 1;
        }
        if (Node_Degree[node] > MAX_VERTEX_NO)
            MAX_VERTEX_NO = Node_Degree[node];

        if (Node_Degree[node] == NB_NODE - p1 - 1) {
            INIT_CLQ_SIZE = NB_NODE - p1;
            //printf("I initial clique is %d\n", INIT_CLQ_SIZE);
            //printf("I maxcore number is %d\n", K_CORE_G);
            MaxCLQ_Stack = (int *) malloc((K_CORE_G + 2) * sizeof(int));
            Clique_Stack = (int *) malloc((K_CORE_G + 2) * sizeof(int));
            memcpy(MaxCLQ_Stack, Candidate_Stack + p1, INIT_CLQ_SIZE * sizeof(int));
            for (i = p1 + 1; i < NB_NODE; i++)
                CORE_NO[i] = cur_degree;
            break;
        }

        neibors = Node_Neibors[node];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (where[neibor] > p1) {
                t = where[neibor];
                h = degree_counter[Node_Degree[neibor]];

                k = Candidate_Stack[h];

                Candidate_Stack[h] = neibor;
                where[neibor] = h;

                Candidate_Stack[t] = k;
                where[k] = t;

                degree_counter[Node_Degree[neibor]]++;

                Node_Degree[neibor]--;
                if (Node_Degree[neibor] != Node_Degree[Candidate_Stack[h - 1]]) {
                    degree_counter[Node_Degree[neibor]] = h;
                }
            }
        }
        p1++;
    }

    if (MAX_VERTEX_NO + 1 == INIT_CLQ_SIZE) {
        MAX_CLQ_SIZE = INIT_CLQ_SIZE;
        //printf("I find the maximum clique in initial phase!\n");
        return true;
    } else {
        return false;
    }
}

bool LMC_Instance::re_number_adj(int node) const {
    int i, k, *neibors, *saved_neibors, neibor, one_neibor;
    for (i = 0; i < iSET_COUNT - 1; i++) {
        neibors = iSET[i];
        one_neibor = NONE;
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (is_adjacent(node, neibor)) {
                if (one_neibor == NONE) {
                    one_neibor = neibor;
                    saved_neibors = neibors;
                } else {
                    break;
                }
            }
        }
        if (one_neibor == NONE) {
            iSET[i][iSET_Size[i]] = node;
            iSET_Size[i]++;
            iSET[i][iSET_Size[i]] = NONE;
            return true;
        }
        if (neibor == NONE) {
            for (k = i + 1; k < iSET_COUNT; k++) {
                neibors = iSET[k];
                for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                    if (one_neibor > neibor) {
                        if (is_adjacent(one_neibor, neibor))
                            break;
                    } else {
                        if (is_adjacent(neibor, one_neibor))
                            break;
                    }
                }
                if (neibor == NONE) {
                    iSET[k][iSET_Size[k]] = one_neibor;
                    iSET_Size[k]++;
                    iSET[k][iSET_Size[k]] = NONE;
                    *saved_neibors = node;
                    return true;
                }
            }
        }
    }
    return false;
}

bool LMC_Instance::re_number(int node) const {
    int i, k, *neibors, *saved_neibors, neibor, one_neibor;
    unsigned char *adj_list1 = iMatrix(node), *adj_list2;
    for (i = 0; i < iSET_COUNT - 1; i++) {
        neibors = iSET[i];
        one_neibor = NONE;
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (Matrix(adj_list1, neibor) > 0) {
                if (one_neibor == NONE) {
                    one_neibor = neibor;
                    saved_neibors = neibors;
                } else {
                    break;
                }
            }
        }
        if (one_neibor == NONE) {
            iSET[i][iSET_Size[i]] = node;
            iSET_Size[i]++;
            iSET[i][iSET_Size[i]] = NONE;
            return true;
        }
        if (neibor == NONE) {
            adj_list2 = iMatrix(one_neibor);
            for (k = i + 1; k < iSET_COUNT; k++) {
                neibors = iSET[k];
                for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                    if (Matrix(adj_list2, neibor) > 0)
                        break;
                }
                if (neibor == NONE) {
                    iSET[k][iSET_Size[k]] = one_neibor;
                    iSET_Size[k]++;
                    iSET[k][iSET_Size[k]] = NONE;
                    iSET_Index[one_neibor] = k;
                    *saved_neibors = node;
                    iSET_Index[node] = i;
                    return true;
                }
            }
        }
    }
    return false;
}

bool LMC_Instance::addIntoIsetTomitaBis_adj(int node) {
    int j, *current_set, iset_node;
    for (j = 0; j < iSET_COUNT; j++) {
        current_set = iSET[j];
        for (iset_node = *current_set; iset_node != NONE; iset_node = *(++current_set)) {
            //assert(node > iset_node);
            if (is_adjacent(node, iset_node))
                break;
        }
        if (iset_node == NONE) {
            iSET_Size[j]++;
            *(current_set) = node;
            *(++current_set) = NONE;
            return true;
        }
    }
    if (iSET_COUNT < MAX_ISET_SIZE) {
        iSET_Size[j] = 1;
        iSET[j][0] = node;
        iSET[j][1] = NONE;
        iSET_COUNT++;
        return true;
    } else {
        return false;
    }
}

bool LMC_Instance::addIntoIsetTomitaBis(int node) {
    int j, *current_set, iset_node;
    unsigned char *adj_list = iMatrix(node);
    for (j = 0; j < iSET_COUNT; j++) {
        current_set = iSET[j];
        for (iset_node = *current_set; iset_node != NONE; iset_node = *(++current_set)) {
            if (Matrix(adj_list, iset_node) > 0)
                break;
        }
        if (iset_node == NONE) {
            iSET_Size[j]++;
            *(current_set) = node;
            *(++current_set) = NONE;
            iSET_Index[node] = j;
            return true;
        }
    }
    if (iSET_COUNT < MAX_ISET_SIZE) {
        iSET_Size[j] = 1;
        iSET[j][0] = node;
        iSET[j][1] = NONE;
        iSET_Index[node] = j;
        iSET_COUNT++;
        return true;
    } else {
        return false;
    }
}

bool LMC_Instance::cut_by_iset_less_vertices() {
    int i = ptr(Candidate_Stack) - 2, node;
    iSET_COUNT = 0;
    MAX_ISET_SIZE = MAX_CLQ_SIZE - CUR_CLQ_SIZE - 1;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[--i]) {
        if (!addIntoIsetTomitaBis_adj(node) && !re_number_adj(node)) {
            return false;
        }
    }
    return true;
}

bool LMC_Instance::cut_by_iset_last_renumber() {
    int i = ptr(Candidate_Stack) - 2, node;
    LAST_IN = INT_MAX;
    FIRST_INDEX = NONE;
    iSET_COUNT = 0;
    MAX_ISET_SIZE = MAX_CLQ_SIZE - CUR_CLQ_SIZE - 1;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[--i]) {
        if (!addIntoIsetTomitaBis(node) && !re_number(node)) {
            if (FIRST_INDEX == NONE)
                FIRST_INDEX = i;
            Candidate_Stack[i] = -node;
        } else {
            LAST_IN = i;
        }
    }
    if (FIRST_INDEX == NONE) {
        Vertex_UB[CURSOR] = iSET_COUNT + 1;
        return true;
    } else {
        Branching_Point = FIRST_INDEX + 1;
        if (MAX_CLQ_SIZE < START_MAXSAT_THD) {
            i = FIRST_INDEX;
            for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[--i]) {
                if (node < 0 && i > LAST_IN) {
                    Vertex_UB[i] = MAX_ISET_SIZE + 1;
                }
            }
        }
        return false;
    }
}

#define assign_node(node, value, reason) \
    Node_State[node] = value;\
    Node_Reason[node] = reason;\
    push(node, FIXED_NODE_STACK)

int LMC_Instance::fix_newNode_for_iset(int fix_node, int fix_iset) {
    int idx, iset_idx;
    iSET_State[fix_iset] = PASSIVE;
    push(fix_iset, PASSIVE_iSET_STACK);
    assign_node(fix_node, P_FALSE, fix_iset);
//iSET_Potential[fix_iset] -= Node_Potential[fix_node];
    idx = ADDED_NODE_iSET[fix_node];
    while ((iset_idx = CONFLICT_ISET_STACK[idx++]) != NONE) {
        if (iSET_State[iset_idx] == ACTIVE) {
            iSET_Size[iset_idx]--;
            //iSET_ADD_Size[iset_idx]--;
            push(iset_idx, REDUCED_iSET_STACK);
            //iSET_Potential[iset_idx] -= Node_Potential[fix_node];
            if (iSET_Size[iset_idx] == 1)
                push(iset_idx, NEW_UNIT_STACK);
            else if (iSET_Size[iset_idx] == 0)
                return iset_idx;
        }
    }
    return NONE;
}

int LMC_Instance::fix_oldNode_for_iset(int fix_node, int fix_iset) {
    int i, node, iset_idx;
    unsigned char *adj_list;
    assign_node(fix_node, P_TRUE, fix_iset);
    iSET_State[fix_iset] = PASSIVE;
    push(fix_iset, PASSIVE_iSET_STACK);

    adj_list = iMatrix(fix_node);
    i = ptr(Candidate_Stack) - 2;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[--i]) {
        if (node > 0 && Matrix(adj_list, node) == 0 && Node_State[node] == ACTIVE) {
            assign_node(node, P_FALSE, fix_iset);
            iset_idx = iSET_Index[node];
            if (iSET_State[iset_idx] == ACTIVE) {
                iSET_Size[iset_idx]--;
                push(iset_idx, REDUCED_iSET_STACK);
                if (iSET_Size[iset_idx] == 1) {
                    push(iset_idx, NEW_UNIT_STACK);
                } else if (iSET_Size[iset_idx] == 0) {
                    return iset_idx;
                }
            }
        }
    }

    return NONE;
}

#define fix_node(node, iset) (((node) > NB_NODE)? fix_newNode_for_iset(node, iset):fix_oldNode_for_iset(node, iset))

int LMC_Instance::fix_node_iset(int fix_iset) {
    int fix_node, *nodes = iSET[fix_iset];
    for (fix_node = *(nodes); fix_node != NONE; fix_node = *(++nodes)) {
        if (Node_State[fix_node] == ACTIVE) {
            if (fix_node > MAX_VERTEX_NO)
                return fix_newNode_for_iset(fix_node, fix_iset);
            else
                return fix_oldNode_for_iset(fix_node, fix_iset);
        }
    }
    nodes = iSET[fix_iset];
    for (fix_node = *(nodes); fix_node != NONE; fix_node = *(++nodes)) {
        printf("iset=%d,node=%d,active=%d\n", fix_iset, fix_node, Node_State[fix_node]);
    }
    printf("error in fix_node_iset\n");
    printf("iSET COUNT=%d\n", iSET_COUNT);
    exit(0);
}

int LMC_Instance::unit_iset_process() {
    int iset_idx, empty_iset;
    for (int i = 0; i < ptr(UNIT_STACK); i++) {
        iset_idx = UNIT_STACK[i];
        if (iSET_State[iset_idx] == ACTIVE && iSET_Size[iset_idx] == 1) {
            ptr(NEW_UNIT_STACK) = 0;
            if ((empty_iset = fix_node_iset(iset_idx)) > NONE) {
                return empty_iset;
            } else {
                for (int j = 0; j < ptr(NEW_UNIT_STACK); j++) {
                    iset_idx = NEW_UNIT_STACK[j];
                    if (iSET_State[iset_idx] == ACTIVE) {
                        if ((empty_iset = fix_node_iset(iset_idx)) > NONE)
                            return empty_iset;
                    }
                }
            }
        }
    }
    ptr(NEW_UNIT_STACK) = 0;
    return NONE;
}

int LMC_Instance::unit_iset_process_used_first() {
    int j, iset, iset_start = 0, used_iset_start = 0, my_iset;
    do {
        for (j = used_iset_start; j < ptr(NEW_UNIT_STACK); j++) {
            iset = NEW_UNIT_STACK[j];
            if (iSET_State[iset] == ACTIVE && iSET_Used[iset])
                if ((my_iset = fix_node_iset(iset)) != NONE)
                    return my_iset;
        }
        used_iset_start = j;
        for (j = iset_start; j < ptr(NEW_UNIT_STACK); j++) {
            iset = NEW_UNIT_STACK[j];
            if (iSET_State[iset] == ACTIVE) {
                if ((my_iset = fix_node_iset(iset)) != NONE)
                    return my_iset;
                iset_start = j + 1;
                break;
            }
        }
    } while (j < ptr(NEW_UNIT_STACK));
    return NONE;
}

void LMC_Instance::identify_conflict_sets(int iset_idx) {
    int i, reason_start = ptr(REASON_STACK), iset, *nodes, node, reason_iset;
    push(iset_idx, REASON_STACK);
    iSET_Involved[iset_idx] = true;
    for (i = reason_start; i < ptr(REASON_STACK); i++) {
        iset = REASON_STACK[i];
        nodes = iSET[iset];
        for (node = *nodes; node != NONE; node = *(++nodes))
            if (Node_State[node] == P_FALSE && Node_Reason[node] != NO_REASON
                && !iSET_Involved[Node_Reason[node]]) {
                reason_iset = Node_Reason[node];
                push(reason_iset, REASON_STACK);
                Node_Reason[node] = NO_REASON;
                iSET_Involved[reason_iset] = true;
            }
    }
    for (i = reason_start; i < ptr(REASON_STACK); i++) {
        iSET_Involved[REASON_STACK[i]] = false;
        iSET_Used[REASON_STACK[i]] = true;
    }
}

void LMC_Instance::enlarge_conflict_sets() {
    int i, iset;
    Node_State[ADDED_NODE] = ACTIVE;
    Node_Reason[ADDED_NODE] = NO_REASON;
//node_match_state[ADDED_NODE] = false;
    ADDED_NODE_iSET[ADDED_NODE] = ptr(CONFLICT_ISET_STACK);
    for (i = 0; i < ptr(REASON_STACK); i++) {
        iset = REASON_STACK[i];
        if (!iSET_Involved[iset]) {
            iSET_Involved[iset] = true;
            iSET[iset][iSET_Size[iset]++] = ADDED_NODE;
            iSET[iset][iSET_Size[iset]] = NONE;
            //assert(ptr(CONFLICT_ISET_STACK)<3*tab_node_size);
            push(iset, CONFLICT_ISET_STACK);
        }
    }
    push(NONE, CONFLICT_ISET_STACK);

    for (i = 0; i < ptr(REASON_STACK); i++) {
        iSET_Involved[REASON_STACK[i]] = false;
        iSET_Used[REASON_STACK[i]] = false;
    }
    ptr(REASON_STACK) = 0;
    ADDED_NODE++;
}

void LMC_Instance::rollback_context_for_maxsatz(int start_fixed, int start_passive, int start_reduced) {
    int i, node;
    for (i = start_fixed; i < ptr(FIXED_NODE_STACK); i++) {
        node = FIXED_NODE_STACK[i];
        Node_State[node] = ACTIVE;
        Node_Reason[node] = NO_REASON;
    }
    ptr(FIXED_NODE_STACK) = start_fixed;
    for (i = start_passive; i < ptr(PASSIVE_iSET_STACK); i++) {
        iSET_State[PASSIVE_iSET_STACK[i]] = ACTIVE;
    }
    ptr(PASSIVE_iSET_STACK) = start_passive;
    for (i = start_reduced; i < ptr(REDUCED_iSET_STACK); i++) {
        iSET_Size[REDUCED_iSET_STACK[i]]++;
    }
    ptr(REDUCED_iSET_STACK) = start_reduced;
    ptr(NEW_UNIT_STACK) = 0;
}

void LMC_Instance::reset_context_for_maxsatz() {
    int i, node;
    for (i = 0; i < ptr(FIXED_NODE_STACK); i++) {
        node = FIXED_NODE_STACK[i];
        Node_State[node] = ACTIVE;
        Node_Reason[node] = NO_REASON;
    }
    ptr(FIXED_NODE_STACK) = 0;
    for (i = 0; i < ptr(PASSIVE_iSET_STACK); i++) {
        iSET_State[PASSIVE_iSET_STACK[i]] = ACTIVE;
    }
    ptr(PASSIVE_iSET_STACK) = 0;
    for (i = 0; i < ptr(REDUCED_iSET_STACK); i++) {
        iSET_Size[REDUCED_iSET_STACK[i]]++;
    }
    ptr(REDUCED_iSET_STACK) = 0;
    ptr(NEW_UNIT_STACK) = 0;
}

int LMC_Instance::further_test_reduced_iset(int start) {
    int i, chosen_iset, empty_iset, node, *nodes;
    int saved_fixed = ptr(FIXED_NODE_STACK);
    int start_passive = ptr(PASSIVE_iSET_STACK);
    int start_reduced = ptr(REDUCED_iSET_STACK);
    for (i = start; i < ptr(REDUCED_iSET_STACK); i++)
        iSET_Tested[REDUCED_iSET_STACK[i]] = false;
    for (i = start; i < ptr(REDUCED_iSET_STACK); i++) {
        chosen_iset = REDUCED_iSET_STACK[i];
        if (!iSET_Tested[chosen_iset]
            && iSET_State[chosen_iset] == ACTIVE
            && iSET_Size[chosen_iset] == 2) {
            iSET_Tested[chosen_iset] = true;
            nodes = iSET[chosen_iset];
            for (node = *nodes; node != NONE; node = *(++nodes)) {
                if (Node_State[node] == ACTIVE && node > MAX_VERTEX_NO)
                    break;
            }
            if (node == NONE) {
                nodes = iSET[chosen_iset];
                for (node = *nodes; node != NONE; node = *(++nodes)) {
                    if (Node_State[node] == ACTIVE) {
                        ptr(NEW_UNIT_STACK) = 0;
                        if ((empty_iset = fix_oldNode_for_iset(node, chosen_iset)) != NONE
                            || (empty_iset = unit_iset_process_used_first()) != NONE) {
                            iSET_Involved[chosen_iset] = true;
                            identify_conflict_sets(empty_iset);
                            iSET_Involved[chosen_iset] = false;
                            rollback_context_for_maxsatz(saved_fixed, start_passive, start_reduced);
                        } else {
                            rollback_context_for_maxsatz(saved_fixed, start_passive, start_reduced);
                            break;
                        }
                    }
                }
                if (node == NONE)
                    return chosen_iset;
            }
        }
    }
    return NONE;
}

int LMC_Instance::fix_anyNode_for_iset(int fix_node, int fix_iset) {
    if (fix_node > MAX_VERTEX_NO)
        return fix_newNode_for_iset(fix_node, fix_iset);
    else
        return fix_oldNode_for_iset(fix_node, fix_iset);
}

bool LMC_Instance::inc_maxsatz_lookahead_by_fl2() {
    int i, j, empty_iset, iset_idx, *nodes, node;
    int sn = FIXED_NODE_STACK_fill_pointer;
    int sp = PASSIVE_iSET_STACK_fill_pointer;
    int sr = REDUCED_iSET_STACK_fill_pointer;
    int rs = REASON_STACK_fill_pointer;
    for (i = 0; i < sr; i++) {
        Is_Tested[REDUCED_iSET_STACK[i]] = false;
    }
    for (i = 0; i < sr; i++) {
        iset_idx = REDUCED_iSET_STACK[i];
        if (!Is_Tested[iset_idx] && iSET_State[iset_idx] == ACTIVE
            && iSET_Size[iset_idx] == 2) {
            nodes = iSET[iset_idx];
            ptr(REASON_STACK) = rs;
            Is_Tested[iset_idx] = true;
            for (node = *nodes; node != NONE; node = *(++nodes)) {
                if (Node_State[node] == ACTIVE) {
                    ptr(NEW_UNIT_STACK) = 0;
                    if ((empty_iset = fix_anyNode_for_iset(node, iset_idx)) != NONE
                        || (empty_iset = unit_iset_process_used_first()) != NONE
                        || (*(nodes + 1) == NONE && (empty_iset = further_test_reduced_iset(sr)) != NONE)) {
                        identify_conflict_sets(empty_iset);
                        rollback_context_for_maxsatz(sn, sp, sr);
                    } else {
                        rollback_context_for_maxsatz(sn, sp, sr);
                        break;
                    }
                }
            }
            if (node == NONE) {
                //reset_context_for_maxsatz();
                return true;
            } else {
                for (j = rs; j < ptr(REASON_STACK); j++) {
                    iSET_Involved[REASON_STACK[j]] = false;
                    iSET_Used[REASON_STACK[j]] = false;
                }
                ptr(REASON_STACK) = rs;
            }
        }
    }
//reset_context_for_maxsatz();
    return false;
}

int LMC_Instance::inc_maxsatz_on_last_iset() {
    int empty_iset, node, *nodes, iset_idx = iSET_COUNT - 1;
    ptr(REASON_STACK) = 0;
    ptr(NEW_UNIT_STACK) = 0;
    ptr(FIXED_NODE_STACK) = 0;
    ptr(PASSIVE_iSET_STACK) = 0;
    ptr(REDUCED_iSET_STACK) = 0;
    nodes = iSET[iset_idx];
    for (node = *nodes; node != NONE; node = *(++nodes)) {
        if (Node_State[node] == ACTIVE) {
            ptr(NEW_UNIT_STACK) = 0;
            if ((empty_iset = fix_oldNode_for_iset(node, iset_idx)) != NONE
                || (empty_iset = unit_iset_process_used_first()) != NONE
                || (empty_iset = unit_iset_process()) != NONE) {
                identify_conflict_sets(empty_iset);
                reset_context_for_maxsatz();
            } else if (inc_maxsatz_lookahead_by_fl2()) {
                reset_context_for_maxsatz();
            } else {
                reset_context_for_maxsatz();
                break;
            }
        }
    }
    if (node == NONE) {
        enlarge_conflict_sets();
    }
    return node;
}

int LMC_Instance::open_new_iset_old(int i) {
    int *current_set, node, iset_node, idx = 0;
    unsigned char *adj_list;
    iSET_Size[iSET_COUNT] = 0;
    iSET[iSET_COUNT][0] = NONE;
    iSET_Used[iSET_COUNT] = false;
    iSET_State[iSET_COUNT] = ACTIVE;

    while ((node = Candidate_Stack[i]) != DELIMITER) {
        if (node < 0) {
            node = -node;
            adj_list = iMatrix(node);
            current_set = iSET[iSET_COUNT];
            for (iset_node = *current_set; iset_node != NONE; iset_node = *(++current_set)) {
                if (Matrix(adj_list, iset_node) > 0)
                    break;
            }
            if (iset_node == NONE) {
                iSET_Size[iSET_COUNT]++;
                *(current_set) = node;
                *(++current_set) = NONE;
                iSET_Index[node] = iSET_COUNT;
                Node_State[node] = ACTIVE;
                Node_Reason[node] = NO_REASON;
                Candidate_Stack[i] = node;
                Extra_Node_Stack[idx++] = node;
                Extra_Node_Stack[idx++] = i;
            } else {
                break;
            }
        }
        i--;
    }
    if (iSET_Size[iSET_COUNT] == 1) {
        push(UNIT_STACK[0], UNIT_STACK);
        UNIT_STACK[0] = iSET_COUNT;
    }
    iSET_COUNT++;
    //assert(iSET_COUNT <= MAX_VERTEX_NO);
    Extra_Node_Stack[idx] = NONE;
    return i;
}

int LMC_Instance::simple_further_test_node(int start) {
    int my_iset, saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer, chosen_iset, node, *nodes, i, j;
    int my_saved_node_stack_fill_pointer,
            my_saved_passive_iset_stack_fill_pointer,
            my_saved_reduced_iset_stack_fill_pointer;
    bool conflict = false;

    saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
    saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
    saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;
    my_saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
    my_saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
    my_saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;

    for (i = start; i < ptr(REDUCED_iSET_STACK); i++)
        iSET_Tested[REDUCED_iSET_STACK[i]] = false;
    for (i = start; i < ptr(REDUCED_iSET_STACK); i++) {
        chosen_iset = REDUCED_iSET_STACK[i];
        if (iSET_State[chosen_iset] == ACTIVE && !iSET_Tested[chosen_iset] && iSET_Size[chosen_iset] <= 2) {
            nodes = iSET[chosen_iset];
            iSET_Tested[chosen_iset] = true;
            for (node = *nodes; node != NONE; node = *(++nodes)) {
                if (node <= MAX_VERTEX_NO && Node_State[node] == ACTIVE) {
                    ptr(NEW_UNIT_STACK) = 0;
                    my_iset = fix_oldNode_for_iset(node, chosen_iset);
                    if (my_iset == NONE)
                        my_iset = unit_iset_process_used_first();
                    rollback_context_for_maxsatz(my_saved_node_stack_fill_pointer,
                                                 my_saved_passive_iset_stack_fill_pointer,
                                                 my_saved_reduced_iset_stack_fill_pointer);
                    if (my_iset != NONE) {
                        assign_node(node, P_FALSE, NO_REASON);
                        iSET_Size[chosen_iset]--;
                        push(chosen_iset, REDUCED_iSET_STACK);
                        if (iSET_Size[chosen_iset] == 1) {
                            ptr(NEW_UNIT_STACK) = 0;
                            push(chosen_iset, NEW_UNIT_STACK);
                            if (unit_iset_process_used_first() != NONE) {
                                conflict = true;
                                break;
                            }
                            for (j = my_saved_reduced_iset_stack_fill_pointer; j < ptr(REDUCED_iSET_STACK); j++)
                                iSET_Tested[REDUCED_iSET_STACK[j]] = false;
                            my_saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
                            my_saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
                            my_saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;
                        }
                    }
                }
            }
            if (conflict)
                break;
        }
    }
    rollback_context_for_maxsatz(saved_node_stack_fill_pointer,
                                 saved_passive_iset_stack_fill_pointer,
                                 saved_reduced_iset_stack_fill_pointer);
    if (conflict)
        return chosen_iset;
    else
        return NONE;
}

int LMC_Instance::test_node_for_failed_nodes(int node, int iset) {
    int my_iset, saved_node_stack_fill_pointer,
            saved_passive_iset_stack_fill_pointer,
            saved_reduced_iset_stack_fill_pointer;

    saved_node_stack_fill_pointer = FIXED_NODE_STACK_fill_pointer;
    saved_passive_iset_stack_fill_pointer = PASSIVE_iSET_STACK_fill_pointer;
    saved_reduced_iset_stack_fill_pointer = REDUCED_iSET_STACK_fill_pointer;
    ptr(NEW_UNIT_STACK) = 0;
    if ((my_iset = fix_oldNode_for_iset(node, iset)) == NONE) {
        if ((my_iset = unit_iset_process_used_first()) == NONE) {
            my_iset = simple_further_test_node(saved_reduced_iset_stack_fill_pointer);
        }
    }
    rollback_context_for_maxsatz(saved_node_stack_fill_pointer,
                                 saved_passive_iset_stack_fill_pointer,
                                 saved_reduced_iset_stack_fill_pointer);
    return my_iset;
}

bool LMC_Instance::test_by_eliminate_failed_nodes() {
    int node, my_iset, *nodes, false_flag;
    bool conflict;
    do {
        false_flag = 0;
        for (my_iset = iSET_COUNT - 1; my_iset >= 0; my_iset--) {
            if (iSET_State[my_iset] == ACTIVE) {
                nodes = iSET[my_iset];
                conflict = false;
                ptr(NEW_UNIT_STACK) = 0;
                for (node = *nodes; node != NONE; node = *(++nodes)) {
                    if (node <= MAX_VERTEX_NO && Node_State[node] == ACTIVE
                        && test_node_for_failed_nodes(node, my_iset) != NONE) {
                        ptr(NEW_UNIT_STACK) = 0;
                        assign_node(node, P_FALSE, NO_REASON);
                        false_flag++;
                        iSET_Size[my_iset]--;
                        push(my_iset, REDUCED_iSET_STACK);
                        if (iSET_Size[my_iset] == 1) {
                            push(my_iset, NEW_UNIT_STACK);
                            break;
                        } else if (iSET_Size[my_iset] == 0) {
                            conflict = true;
                            break;
                        }
                    }
                }
                if (conflict)
                    break;
                else if (ptr(NEW_UNIT_STACK) > 0 && unit_iset_process_used_first() != NONE) {
                    conflict = true;
                    break;
                }
            }
        }
    } while (false_flag > 1 && !conflict);
    reset_context_for_maxsatz();
    return conflict;
}

bool LMC_Instance::cut_by_inc_maxsat_eliminate_first() {
    int i, j, k, node;
    ADDED_NODE = MAX_VERTEX_NO + 1;
    ptr(CONFLICT_ISET_STACK) = 0;
    ptr(UNIT_STACK) = 0;

    i = ptr(Candidate_Stack) - 2;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[--i]) {
        if (node > 0)
            Node_State[node] = ACTIVE;
        else
            Node_State[-node] = PASSIVE;
    }

    for (i = 0; i < iSET_COUNT; i++) {
        if (iSET_Size[i] == 1)
            push(i, UNIT_STACK);
    }

    while ((node = Candidate_Stack[FIRST_INDEX]) != DELIMITER) {
        j = open_new_iset_old(FIRST_INDEX);
        if (Candidate_Stack[j] == DELIMITER && test_by_eliminate_failed_nodes()) {
            FIRST_INDEX = j;
            break;
        } else if ((node = inc_maxsatz_on_last_iset()) == NONE) {
            FIRST_INDEX = j;
        } else {
            for (k = 0; Extra_Node_Stack[k] != node; k += 2);
            FIRST_INDEX = Extra_Node_Stack[k + 1];
            Branching_Point = FIRST_INDEX + 1;
            //SoMC
//			for (k = FIRST_INDEX; Candidate_Stack[k] != DELIMITER; k--) {
//				if (Candidate_Stack[k] > 0)
//					Candidate_Stack[k] = -Candidate_Stack[k];
//			}
            //DoMC
            for (; Extra_Node_Stack[k] != NONE; k += 2)
                Candidate_Stack[Extra_Node_Stack[k + 1]] = -Extra_Node_Stack[k];

            for (k = FIRST_INDEX; Candidate_Stack[k] != DELIMITER; k--) {
                if (Candidate_Stack[k] < 0 && k > LAST_IN)
                    Vertex_UB[k] = MAX_ISET_SIZE + 1;
            }
            break;
        }
    }

    i = ptr(Candidate_Stack) - 2;
    while ((node = Candidate_Stack[i--]) != DELIMITER) {
        if (node > 0)
            Node_State[node] = PASSIVE;
        else
            Node_State[-node] = PASSIVE;
    }

    for (node = MAX_VERTEX_NO; node <= ADDED_NODE; node++) {
        Node_State[node] = PASSIVE;
    }

    if (Candidate_Stack[FIRST_INDEX] == DELIMITER) {
        Vertex_UB[CURSOR] = MAX_CLQ_SIZE - CUR_CLQ_SIZE;
        return true;
    }
    return false;
}

void LMC_Instance::compute_subgraph_degree(int start) {
    int i, j = 0, node, neibor, *neibors;
    Max_Degree = 0;
    int nb_node = 0, nb_edge = 0;

    for (node = Candidate_Stack[i = start]; node != DELIMITER; node = Candidate_Stack[++i]) {
        Node_State[node] = ACTIVE;
        Node_Degree[node] = 0;
        iSET_Index[node] = j;
        iSET_Size[j++] = 0;
        nb_node++;
    }
    for (node = Candidate_Stack[i = start]; node != DELIMITER; node = Candidate_Stack[++i]) {
        neibors = Node_Neibors[node];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (Node_State[neibor] == ACTIVE) {
                nb_edge++;
                iSET[iSET_Index[node]][Node_Degree[node]++] = neibor;
                iSET[iSET_Index[neibor]][Node_Degree[neibor]++] = node;
            }
        }
    }

    for (node = Candidate_Stack[i = start]; node != DELIMITER; node = Candidate_Stack[++i]) {
        iSET[iSET_Index[node]][Node_Degree[node]] = NONE;
        Node_State[node] = PASSIVE;
        if (Node_Degree[node] > Max_Degree)
            Max_Degree = Node_Degree[node];
    }
}

void LMC_Instance::allocate_memory_for_maxsat() {
    Node_Reason = (int *) malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
    ADDED_NODE_iSET = (int *) malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));
    FIXED_NODE_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));

    iSET_State = (bool *) calloc((MAX_VERTEX_NO + 1), sizeof(bool));
    iSET_Used = (bool *) calloc((MAX_VERTEX_NO + 1), sizeof(bool));
    iSET_Tested = (bool *) malloc((MAX_VERTEX_NO + 1) * sizeof(bool));
    UNIT_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
    NEW_UNIT_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
    PASSIVE_iSET_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
    iSET_Involved = (bool *) calloc((MAX_VERTEX_NO + 1), sizeof(bool));
    CONFLICT_ISET_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
    REASON_STACK = (int *) malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
    Is_Tested = (bool *) malloc((MAX_VERTEX_NO + 1) * sizeof(bool));
}

void LMC_Instance::store_maximum_clique(int node) {
    if (CUR_CLQ_SIZE == 0)
        push(node, Clique_Stack);
    else
        push(Second_Name[node], Clique_Stack);
    MAX_CLQ_SIZE = ptr(Clique_Stack);
    memcpy(MaxCLQ_Stack, Clique_Stack, MAX_CLQ_SIZE * sizeof(int));
    ptr(Candidate_Stack) = NB_NODE + 1;
    ptr(Cursor_Stack) = 1;
    ptr(Clique_Stack) = 0;
    Rollback_Point = 0;
    Vertex_UB[CURSOR] = MAX_CLQ_SIZE;
    if (MAX_CLQ_SIZE == START_MAXSAT_THD)
        allocate_memory_for_maxsat();
}

void LMC_Instance::store_maximum_clique2() {
    MAX_CLQ_SIZE = ptr(Clique_Stack);
    memcpy(MaxCLQ_Stack, Clique_Stack, MAX_CLQ_SIZE * sizeof(int));
//ptr(Cursor_Stack) = 1;
    ptr(Clique_Stack) = 0;
    if (MAX_CLQ_SIZE == START_MAXSAT_THD)
        allocate_memory_for_maxsat();
}

bool LMC_Instance::reduce_subgraph(int start) {
    int *degree_counter, *where;
    int end, p1, i, node = NONE, neibor, *neibors, t, j, h, k;
    int max_degree = 0;
    compute_subgraph_degree(start);
    where = Candidate_Stack + ptr(Candidate_Stack) + 1;
    degree_counter = Vertex_UB + ptr(Candidate_Stack) + 1;
    memset(degree_counter, 0, (Max_Degree + 1) * sizeof(int));

    for (node = Candidate_Stack[i = start]; node != DELIMITER; node = Candidate_Stack[++i]) {
        degree_counter[Node_Degree[node]]++;
        Vertex_UB[i] = node;
    }
    Vertex_UB[i] = DELIMITER;

    end = i - 1;
    j = start;
    for (i = 0; i <= Max_Degree; i++) {
        k = degree_counter[i];
        degree_counter[i] = j;
        j += k;
    }

    for (node = Vertex_UB[i = start]; node != DELIMITER; node = Vertex_UB[++i]) {
        t = degree_counter[Node_Degree[node]]++;
        Candidate_Stack[t] = node;
        where[node] = t;
    }

    for (i = Max_Degree; i > 0; i--) {
        degree_counter[i] = degree_counter[i - 1];
    }
    degree_counter[0] = start;
    //return false;
    p1 = start;
    while ((node = Candidate_Stack[p1]) != DELIMITER) {
        if (Node_Degree[node] > max_degree)
            max_degree = Node_Degree[node];
        if (p1 < end
            && Node_Degree[node] == Node_Degree[Candidate_Stack[p1 + 1]]) {
            degree_counter[Node_Degree[node]] = p1 + 1;
        }
        if (Node_Degree[node] == end - p1) {
            if (end - p1 + 1 == MAX_CLQ_SIZE) {
                ptr(Clique_Stack) = 0;
                push(Candidate_Stack[CURSOR], Clique_Stack);
                while ((node = Candidate_Stack[p1++]) != DELIMITER)
                    push(node, Clique_Stack);
                store_maximum_clique2();
                return true;
            } else {
                while ((node = Candidate_Stack[++p1]) != DELIMITER)
                    Node_Degree[node] = end - p1;
            }
            break;
        }

        neibors = iSET[iSET_Index[node]];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors))
            if (where[neibor] > p1) {
                t = where[neibor];
                h = degree_counter[Node_Degree[neibor]];
                k = Candidate_Stack[h];

                Candidate_Stack[h] = neibor;
                where[neibor] = h;
                Candidate_Stack[t] = k;
                where[k] = t;

                degree_counter[Node_Degree[neibor]]++;
                Node_Degree[neibor]--;
                if (Node_Degree[neibor]
                    != Node_Degree[Candidate_Stack[h - 1]]) {
                    degree_counter[Node_Degree[neibor]] = h;
                }
            }

        p1++;
    }
    if (max_degree + 1 >= MAX_CLQ_SIZE) {
        CUR_MAX_NODE = 0;
        for (node = Candidate_Stack[i = start]; node != DELIMITER; node = Candidate_Stack[++i]) {
            if (Node_Degree[node] + 1 >= MAX_CLQ_SIZE)
                break;
        }
        j = start;
        for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[++i]) {
            Vertex_UB[j] = Node_Degree[node] + 1;
            Candidate_Stack[j++] = node;
            if (node > CUR_MAX_NODE)
                CUR_MAX_NODE = node;
        }
        Candidate_Stack[j] = DELIMITER;
        ptr(Candidate_Stack) = j + 1;
        return false;
    } else {
        Vertex_UB[CURSOR] = max_degree + 2;
        return true;
    }
}

void LMC_Instance::rebuild_matrix(int start) {
    int i = start, j = 1, node, neibor, *neibors;
//	for (node = 1; node <= NB_NODE; node++)
//		Node_State[node] = PASSIVE;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[++i]) {
        Candidate_Stack[i] = j;
        Second_Name[j] = node;
        Node_Degree[node] = j++;
        Node_State[node] = ACTIVE;
    }
    memset(Adj_Matrix, 0, (MAX_VERTEX_NO + 1) * (MATRIX_ROW_WIDTH) * sizeof(char));
    i = start;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[++i]) {
        neibors = Node_Neibors[Second_Name[node]];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (Node_State[neibor] == ACTIVE) {
                SET_EDGE(node, Node_Degree[neibor]);
                SET_EDGE(Node_Degree[neibor], node);
            }
        }
        //	Node_State[Second_Name[node]] = PASSIVE;
    }
    i = start;
    for (node = Candidate_Stack[i]; node != DELIMITER; node = Candidate_Stack[++i]) {
        Node_State[Second_Name[node]] = PASSIVE;
    }
}

bool LMC_Instance::cut_by_inc_ub() {
    int i = CURSOR, neibor, max = 0, *neibors;
    int node = Candidate_Stack[CURSOR];
    int start = ptr(Candidate_Stack);
    unsigned char *adj_list;
    NB_CANDIDATE = 0;
    CUR_MAX_NODE = 0;
    if (CUR_CLQ_SIZE == 0) {
        neibors = Node_Neibors[node];
        CUR_MAX_NODE = *(neibors);
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            i = NB_NODE - neibor;
            if (max < Vertex_UB[i])
                max = Vertex_UB[i];
            Vertex_UB[ptr(Candidate_Stack)] = Vertex_UB[i];
            push(neibor, Candidate_Stack);
            NB_CANDIDATE++;
        }
    } else {
        adj_list = iMatrix(node);
        while (Candidate_Stack[i] != DELIMITER)i--;
        for (neibor = Candidate_Stack[++i]; neibor != DELIMITER; neibor = Candidate_Stack[++i]) {
            if (neibor > 0 && Matrix(adj_list, neibor) > 0) {
                if (max < Vertex_UB[i])
                    max = Vertex_UB[i];
                Vertex_UB[ptr(Candidate_Stack)] = Vertex_UB[i];
                push(neibor, Candidate_Stack);
                NB_CANDIDATE++;
            }
        }
    }
    push(DELIMITER, Candidate_Stack);

    /*max=NB_CANDIDATE;//disable incremental upper bound*/
    if (NB_CANDIDATE < max) {
        max = NB_CANDIDATE;
    }
    if (Vertex_UB[CURSOR] - 1 < max) {
        max = Vertex_UB[CURSOR] - 1;
    }
    if (max < MAX_CLQ_SIZE - CUR_CLQ_SIZE) {
        Vertex_UB[CURSOR] = max + 1;
        return true;
    } else if (CUR_CLQ_SIZE == 0) {
        if (NB_CANDIDATE < 10 && cut_by_iset_less_vertices()) {
            Vertex_UB[CURSOR] = iSET_COUNT + 1;
            return true;
        } else if (reduce_subgraph(start)) {
            return true;
        } else if (REBUILD_MATRIX || CUR_MAX_NODE > MAX_VERTEX_NO) {
            rebuild_matrix(start);
            REBUILD_MATRIX = true;
            return false;
        }
        return false;
    } else {
        return false;
    }
}

bool LMC_Instance::find_3_clique(int node) {
    int neibor1, neibor2, neibor3, *neibors1, *neibors2, *neibors3;
    neibors1 = Node_Neibors[node];
    for (neibor1 = *neibors1; neibor1 != NONE; neibor1 = *(++neibors1)) {
        neibors2 = neibors1 + 1;
        for (neibor2 = *neibors2; neibor2 != NONE; neibor2 = *(++neibors2)) {
            if (is_adjacent(neibor1, neibor2)) {
                neibors3 = neibors2 + 1;
                for (neibor3 = *neibors3; neibor3 != NONE; neibor3 = *(++neibors3)) {
                    if (is_adjacent(neibor1, neibor3) && is_adjacent(neibor2, neibor3)) {
                        MaxCLQ_Stack[0] = Old_Name[node];
                        MaxCLQ_Stack[1] = Old_Name[neibor1];
                        MaxCLQ_Stack[2] = Old_Name[neibor2];
                        MaxCLQ_Stack[3] = Old_Name[neibor3];
                        MAX_CLQ_SIZE = 4;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void LMC_Instance::init_for_search() {
    int i, node;
    int neibor, neibor2, *neibors, *neibors2;

    MAX_CLQ_SIZE = 0;
    ptr(Clique_Stack) = 0;
    ptr(Cursor_Stack) = 0;
    Rollback_Point = 0;
    push(NB_NODE - INIT_CLQ_SIZE - 1, Cursor_Stack);
    MAX_CLQ_SIZE = INIT_CLQ_SIZE;
    for (i = 0; i < ptr(Candidate_Stack) - 1; i++) {
        node = Candidate_Stack[i];
        Vertex_UB[i] = Node_Degree[node] + 1;
        if (INIT_CLQ_SIZE == 3 && Vertex_UB[i] > 3) {
            if (find_3_clique(node)) {
                //printf("I improve the initial clique to 4\n");
                INIT_CLQ_SIZE = 4;
            } else {
                Vertex_UB[i] = 3;
            }
        }
        if (INIT_CLQ_SIZE == 4 && Vertex_UB[i] == 5) {
            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                neibors2 = neibors + 1;
                for (neibor2 = *neibors2; neibor2 != NONE; neibor2 = *(++neibors2)) {
                    if (!is_adjacent(neibor, neibor2)) {
                        Vertex_UB[i] = 4;
                        break;
                    }
                }
                if (Vertex_UB[i] == 4)
                    break;
            }
        }
    }
}

void LMC_Instance::allocate_memory() {
    int i;
    Second_Name = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
    iSET = (int **) malloc((MAX_VERTEX_NO + 1) * sizeof(int *));
    iSET[0] = (int *) malloc((MAX_VERTEX_NO + 1) * (MAX_VERTEX_NO + 1) * sizeof(int));
    for (i = 1; i < MAX_VERTEX_NO; i++) {
        iSET[i] = iSET[i - 1] + MAX_VERTEX_NO + 1;
    }
    iSET_Size = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
    iSET_Index = (int *) malloc((NB_NODE + 1) * sizeof(int));

    if (INIT_CLQ_SIZE >= START_MAXSAT_THD)
        allocate_memory_for_maxsat();
}

void LMC_Instance::search_maxclique() {
    int node;
    init_for_search();
    while (CURSOR > 0) {
        node = Candidate_Stack[--CURSOR];
        if (CUR_CLQ_SIZE > 0 && node > 0)
            continue;
        if (node == DELIMITER) {
            ptr(Candidate_Stack) = CURSOR + 1;
            ptr(Cursor_Stack)--;
            ptr(Clique_Stack)--;
            Vertex_UB[CURSOR] = MAX_CLQ_SIZE - CUR_CLQ_SIZE;
        } else {
            if (node < 0) {
                node = -node;
                Candidate_Stack[CURSOR] = -Candidate_Stack[CURSOR];
            }
            if (MAX_CLQ_SIZE == CUR_CLQ_SIZE) {
                store_maximum_clique(node);
            } else if (Vertex_UB[CURSOR] > MAX_CLQ_SIZE - CUR_CLQ_SIZE) {
                Rollback_Point = ptr(Candidate_Stack);

                if (cut_by_inc_ub()
                    || cut_by_iset_last_renumber()
                    || (MAX_CLQ_SIZE >= START_MAXSAT_THD && cut_by_inc_maxsat_eliminate_first())) {
                    ptr(Candidate_Stack) = Rollback_Point;
                } else {
                    if (CUR_CLQ_SIZE == 0)
                        push(node, Clique_Stack);
                    else
                        push(Second_Name[node], Clique_Stack);
                    push(Branching_Point, Cursor_Stack);
                }
            }
        }
    }
}

std::vector<int> LMC_Instance::getMaxCliqueVector() const {
    std::vector<int> maxClique(MAX_CLQ_SIZE);
    if (INIT_CLQ_SIZE < MAX_CLQ_SIZE) {
        for (int i = 0; i < MAX_CLQ_SIZE; i++)
            maxClique[i] = Old_Name[MaxCLQ_Stack[i]] - 1;
    } else {
        for (int i = 0; i < MAX_CLQ_SIZE; i++)
            maxClique[i] = MaxCLQ_Stack[i] - 1;
    }
    return maxClique;
}

void LMC_Instance::build_init_matrix() {
    int node, neibor, *neibors;
    MATRIX_ROW_WIDTH = MAX_VERTEX_NO / 8 + 1;
    Adj_Matrix = (unsigned char *) malloc((MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

    memset(Adj_Matrix, 0, (MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

    for (node = 1; node <= MAX_VERTEX_NO; node++) {
        Second_Name[node] = node;
        neibors = Node_Neibors[node];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            SET_EDGE(node, neibor);
            SET_EDGE(neibor, node);
        }
    }
}

#define New_Name Node_Degree

bool LMC_Instance::search_in_2_k_core_graph() {
    int i, node, neibor1, neibor2, neibor, *neibors;
    bool find = false;
    for (i = 0; i < NB_NODE; i++) {
        node = Candidate_Stack[i];
        if (Node_Degree[node] == 2) {
            neibor1 = Node_Neibors[node][0];
            neibor2 = Node_Neibors[node][1];

            neibors = Node_Neibors[neibor1];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                if (neibor == neibor2) {
                    find = true;
                    break;
                }
            }
            if (find) {
                MaxCLQ_Stack[0] = Old_Name[node];
                MaxCLQ_Stack[1] = Old_Name[neibor1];
                MaxCLQ_Stack[2] = Old_Name[neibor2];
                INIT_CLQ_SIZE = 3;
                MAX_CLQ_SIZE = 3;
                break;
            }
        } else if (Node_Degree[node] > 2) {
            break;
        }
    }
    if (find) {
        //printf("I find maximum clique in initial phase.");
        return true;
    } else if (i == NB_NODE) {
        //printf("I find maximum clique in initial phase.");
        MAX_CLQ_SIZE = 2;
        return true;
    } else {
        return false;
    }
}

void LMC_Instance::free_block() {
    for (int i = 0; i < BLOCK_COUNT; i++)
        free(BLOCK_LIST[i]);
}

void LMC_Instance::reduce_instance() {
    int i, j, nb, node, *neibors, *neibors2, *addr;
    MAX_VERTEX_NO = 0;
    j = 0;
    for (i = 0; i < NB_NODE; i++) {
        node = Candidate_Stack[i];
        if (CORE_NO[i] < INIT_CLQ_SIZE) {
            Node_State[node] = PASSIVE;
        } else {
            Candidate_Stack[j++] = node;
            Node_State[node] = ACTIVE;
        }
    }
    NB_NODE = j;
    Candidate_Stack[j] = DELIMITER;
    ptr(Candidate_Stack) = j + 1;
    Old_Name = (int *) malloc((NB_NODE + 1) * sizeof(int));
    for (i = 0; i < NB_NODE; i++) {
        Old_Name[NB_NODE - i] = Candidate_Stack[i];
        New_Name[Candidate_Stack[i]] = NB_NODE - i;
        Candidate_Stack[i] = NB_NODE - i;
    }

    NB_EDGE = 0;
    for (i = NB_NODE; i > 0; i--) {
        neibors = Node_Neibors[Old_Name[i]];
        neibors2 = neibors;
        nb = 0;
        for (node = *neibors; node != NONE; node = *(++neibors)) {
            if (Node_State[node] == ACTIVE && New_Name[node] < i) {
                (*neibors2) = New_Name[node];
                neibors2++;
                nb++;
            }
        }
        (*neibors2) = NONE;
        NB_EDGE += nb;
        qsort(Node_Neibors[Old_Name[i]], nb, sizeof(int), int_cmp);
    }

    Adj_List = (int *) malloc((NB_EDGE + NB_NODE) * sizeof(int));
    addr = Adj_List;

    if (Adj_List == nullptr) {
        printf("can't allocate memory for Adj_List!\n");
        exit(0);
    }

    for (i = NB_NODE; i > 0; i--) {
        Node_Degree[i] = 0;
        Node_State[i] = PASSIVE;
        neibors = Node_Neibors[Old_Name[i]];
        for (node = *neibors; node != NONE; node = *(++neibors)) {
            *(addr++) = node;
            Node_Degree[i]++;
        }
        *(addr++) = NONE;
        if (Node_Degree[i] > MAX_VERTEX_NO)
            MAX_VERTEX_NO = Node_Degree[i];
    }
    free_block();
    Node_Neibors[NB_NODE] = Adj_List;
    for (i = NB_NODE - 1; i > 0; i--) {
        Node_Neibors[i] = Node_Neibors[i + 1] + Node_Degree[i + 1] + 1;
    }
    //printf("I the reduced graph #node %d #edge %d #density %9.8f\n", NB_NODE, NB_EDGE, ((float) NB_EDGE * 2 / NB_NODE / (NB_NODE - 1)));
    //printf("I the largest subgraph is %d\n", MAX_VERTEX_NO);

}

bool LMC_Instance::initialize() {
    bool r = sort_by_degeneracy_ordering();
    if (!r) {
        reduce_instance();
        if (K_CORE_G <= 2) {
            r = search_in_2_k_core_graph();
        }
    } else {
        free_block();
    }
    return r;
}

std::vector<int> getMaximumClique(const std::vector<std::vector<int>> &adjacencyList) {
    auto lmc = new LMC_Instance();

    lmc->import_graph_instance(adjacencyList);
    if (!lmc->initialize()) {
        lmc->allocate_memory();
        lmc->build_init_matrix();
        lmc->search_maxclique();
    }

    auto maxClique = lmc->getMaxCliqueVector();
    delete lmc;
    return maxClique;
}

std::vector<int> getInitClique(const std::vector<std::vector<int>> &adjacencyList) {
    auto lmc = new LMC_Instance();

    lmc->import_graph_instance(adjacencyList);
    lmc->initialize();

    std::vector<int> initClique(lmc->INIT_CLQ_SIZE);
    for (int i = 0; i < lmc->INIT_CLQ_SIZE; i++)
        initClique[i] = lmc->MaxCLQ_Stack[i] - 1;

    delete lmc;
    return initClique;
}
