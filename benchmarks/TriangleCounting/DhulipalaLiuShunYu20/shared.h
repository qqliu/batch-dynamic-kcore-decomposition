#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
// #include "sparse_table.h"
// #include "gbbs/macros.h"

using namespace std;

#define EMPTYV numeric_limits<uintE>::max()
#define EMPTYVMAP numeric_limits<size_t>::max()
#define EMPTYKVB make_tuple(EMPTYV, 0)
#define EMPTYKV make_tuple(EMPTYV, (SetT *)NULL)
#define EMPTYWTV numeric_limits<uintE>::max()
#define UPDATET1 1
#define UPDATET2 2
#define UPDATET3 3
#define UPDATET4 4
#define UPDATET5 5
#define UPDATECLEAR 6
#define UPDATECLEANUP 7
#define REMOVE1 8

#define OLD_EDGE 0
#define NEW_EDGE 1
#define DEL_EDGE 2
#define NO_EDGE -1


namespace gbbs{
namespace DBTGraph{
    constexpr const size_t smallTasksForThreshold = kDefaultGranularity;
    
    using EdgeT = pair<uintE, uintE>;
    using SymGraph = symmetric_graph<symmetric_vertex, gbbs::empty>;
    using StaticEdgeT = tuple<uintE, gbbs::empty>;

    inline uintE getFirst(sequence<pair<EdgeT,bool>> &edges, size_t i){
    return edges[i].first.first;
    }

    inline uintE getSecond(sequence<pair<EdgeT,bool>> &edges, size_t i){
    return edges[i].first.second;
    }

    inline uintE getFirst(pbbslib::range<pair<EdgeT,bool>> edges, size_t i){
    return edges[i].first.first;
    }

    inline uintE getSecond(pbbslib::range<pair<EdgeT,bool>> edges, size_t i){
    return edges[i].first.second;
    }

    inline uintE getFirst(pair<EdgeT,bool> e){
    return e.first.first;
    }

    inline uintE getSecond(pair<EdgeT,bool> e){
    return e.first.second;
    }

    inline uintE getFirst(pair<uintE, int> e){
    return e.first;
    }

    inline uintE getSecond(pair<uintE, int> e){
    return e.second;
    }

    template <class W>
    inline uintE getFirst(tuple<uintE,W> e){
    return get<0>(e);
    }

    template <class W>
    inline uintE getSecond(tuple<uintE,W> e){
    return get<1>(e);
    }

    template<class SetT>
    struct UpperTable{
        sequence<SetT*> tb;//caution copy
        using K = uintE;
        using V = SetT*;
        using T = std::tuple<K, V>;
        using KT = K;

        UpperTable(size_t n){
            tb = sequence<SetT*>(n, (SetT*)nullptr); // must initialize
        }

        // UpperTable(size_t n, T _empty, K _tomb_key, KeyHash _key_hash, long inp_space_mult=-1){
        //     tb = sequence<SetT*>(n, (SetT*)nullptr); // must initialize
        // }
        
        inline size_t size() const {return tb.n;}

        // bool not_empty(K k){}

        // bool not_empty(T kv){}
        
        // static void clearA(T* A, long n, T kv) {
        //     parallel_for(0, n, [&] (size_t i) { A[i] = kv; });
        // }

        inline void clear() {
            parallel_for(0, tb.n, [&] (size_t i) { tb[i] = nullptr; });
        }

        inline void del() {
            tb.clear();
        }

        // undetermined behavior if delete + insert/update same key
        // at the same time
        inline bool deleteVal(K k) {
            V v = tb[k];
            if(v == nullptr) return false;
            return pbbslib::CAS(&tb[k], v, (SetT*) nullptr);
        }

        inline bool insert(std::tuple<K, V> kv){
            K k = get<0>(kv);
            V v = get<1>(kv);
            return pbbslib::CAS(&tb[k], (SetT*) nullptr, v);
        }

        // template <class F>
        // inline bool insert_f(std::tuple<K, V> kv, const F& f) {
        //     K k = get<0>(kv);
        //     V v = get<1>(kv);
        //     SetT** addr = &tb[k];
        //     if (pbbslib::CAS(addr, nullptr, v)) {
        //         // f(&std::get<1>(table[h]), kv);
        //         f(addr, kv);
        //         return true;
        //     }else{
        //         f(addr, kv);
        //         return false;
        //     }
            
        // }

        inline void updateSeq(K k, V val) {
            tb[k] = val;
        }

        inline void maybe_resize(size_t nt) {}

        inline void maybe_resize(size_t n_inc, size_t ne){}

        inline bool contains(K k) const {return tb[k] != nullptr;}

        inline V find(K k, V default_value) const {
            if(tb[k] == nullptr) return default_value;
            return tb[k];
        }



    };


    struct VtxUpdate{
        uintE id = -1;
        size_t degree = 0; // number of updates
        size_t insert_degree = 0; // number of insertion updates
        size_t insert_low_degree = 0; // number of insertion updates to low
        size_t offset = -1; // offsets in edges
        size_t delete_low_degree = 0;
        bool change_status = false;

        VtxUpdate(uintE a, size_t o):id(a), offset(o){}
        VtxUpdate(uintE a):id(a){}
        inline void setDeg(size_t a){degree = a;}
        inline void setInsDeg(size_t a){insert_degree = a;}
        inline size_t end() const {return offset + degree;}
        inline size_t insOffset() const {return offset + insert_degree;}
        inline size_t delDeg() const {return degree-insert_degree;}
        inline size_t newDeg(size_t oldDeg) const {return oldDeg + 2*insert_degree - degree;}
        inline size_t newLowDeg(size_t oldDeg) const {return oldDeg + insert_low_degree - delete_low_degree;}
        inline bool ins_low(){return insert_low_degree > 0;}
        inline bool ins_high(){return insert_low_degree < insert_degree;}

    };

    struct VtxUpdateInsDeg{
        sequence<VtxUpdate> &vtxNew;
        VtxUpdateInsDeg(sequence<VtxUpdate>& a):vtxNew(a){}

        size_t operator ()(size_t i)const {
            return vtxNew[i].insert_degree;
        }
    };

    template <class SetT>
    struct MakeEdge{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdge(uintE uu, T*t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        EdgeT operator ()(size_t i)const {
            return EdgeT(get<0>(table[i]),u);
        }
    };

    template <class SetT>
    struct MakeEdgeEntry{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeEntry(uintE uu, T*t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        std::pair<uintE, int> operator ()(size_t i) const {
            return make_pair(get<0>(table[i]),get<1>(table[i]));
        }
    };

    template <class SetT>
    struct MakeEdgeEntryMajor{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeEntryMajor(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        StaticEdgeT operator ()(size_t i)const {
            if (get<1>(table[i]) == DEL_EDGE) return make_pair(empty_key, gbbs::empty());
            return make_tuple(get<0>(table[i]), gbbs::empty());
        }
    };

    template <class SetT>
    struct MakeEdgeLtoH{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeLtoH(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        pair<EdgeT, bool> operator ()(size_t i)const {
            return make_pair(EdgeT(get<0>(table[i]),u), true);
        }
    };

    template <class SetT>
    struct MakeEdgeHtoL{
        using T = typename SetT::T; using K = typename SetT::KT;
        uintE u;T *table;K empty_key;
        MakeEdgeHtoL(uintE uu, T* t_table, K _empty):u(uu), table(t_table), empty_key(_empty){}

        pair<EdgeT, bool> operator ()(size_t i)const {
            return make_pair(EdgeT(get<0>(table[i]),u), false);
        }
    };

    struct VtxRbl{
        uintE id = -1;
        size_t LtoH = 0;
        size_t degree = 0;
        size_t offset = -1; // offsets in Ngh

        VtxRbl(uintE a, size_t o):id(a), offset(o){}
        VtxRbl(uintE a):id(a){}
        inline void setDeg(size_t a){degree = a;}
        inline void setInsDeg(size_t a){LtoH = a;}
        inline size_t newLowDeg(size_t oldDeg) const {
            return oldDeg + degree - LtoH - LtoH;
        }
        inline size_t getHtoL() const {return degree-LtoH;}
        inline size_t insOffset() const {return offset + LtoH;}
        inline size_t end() const {return offset + degree;}

    };

    struct TriangleCounts{ // cache line size 64
        sequence<size_t> c; //c1, c2, c3, c4, c5, c6 ... padding to 64; c1, c2, c3, c4, c5, c6; ....padding to 64;
        size_t P;
        static constexpr size_t eltsPerCacheLine = 128 /sizeof(size_t);
        
        //TriangleCounts():c1(0),c2(0),c3(0), c4(0), c5(0), c6(0){
        TriangleCounts(){
            P = num_workers();
            c = sequence<size_t>::uninitialized(P * eltsPerCacheLine);
            clear();
        }

        //flag \in {1,2,3}
        inline void increment(int flag, size_t val){
            c[eltsPerCacheLine * worker_id() + flag - 1 ] += val;
        }

        //flag \in {1,2,3}
        inline void decrement(int flag, size_t val){
            c[eltsPerCacheLine * worker_id() + 2 + flag] += val;
        }

        inline sequence<size_t> report(){
            sequence<size_t> result = sequence<size_t>(6, (size_t)0);
            par_for(0, 6, [&] (size_t j){
                for(size_t i=0;i<P;++i){
                    result[j] += c[eltsPerCacheLine * i + j];
                }
            });
            return result;
        }

        inline void clear(){
            // par_for(0, 6 * P, [&] (size_t i) {c[i] = 0;});
            par_for(0, P, [&] (size_t i){
                for(int j=0;j<6;++j){c[eltsPerCacheLine * i + j] = 0;}
            });
        }
    };

    

    struct WTV{
        volatile uintE c1, c2, c3, c4, c5;
        volatile bool changing = false;
        WTV():c1(0),c2(0),c3(0), c4(0), c5(0){}
        WTV(uintE a):c1(a),c2(a),c3(a), c4(a), c5(a){}

        // WTV(size_t cc1, size_t cc2, size_t cc3):c1(cc1),c2(cc2),c3(cc3){}
        WTV(size_t flag, size_t val):c1(flag),c2(val),c3(0), c4(0), c5(0){}

        // inline size_t getFlag(){return c1;}
        inline size_t getUpdateVal(){return c2;}
        inline void update(const  std::tuple<EdgeT, WTV>& kv){ //(edge key, flag, val, 0)
            size_t flag  = std::get<1>(kv).c1;
            switch(flag) {
            case UPDATET1:
                pbbslib::write_add(&c1, std::get<1>(kv).c2);
                break;
            case REMOVE1:
                pbbslib::write_minus(&c1, std::get<1>(kv).c2);
                break;
            case UPDATET2:
                pbbslib::write_add(&c2, std::get<1>(kv).c2);
                break;
            case UPDATET3:
                pbbslib::write_add(&c3, std::get<1>(kv).c2);
                break;
            case UPDATET4:
                pbbslib::write_add(&c4, std::get<1>(kv).c2);
                break;
            case UPDATET5:
                pbbslib::write_add(&c5, std::get<1>(kv).c2);
                break;
            case UPDATECLEANUP:
                cleanUp();
                break;
            case UPDATECLEAR:
                if(c1 !=0) c1 = 0;
                break;
            default:
                cout << "invalid update flag " << flag << endl;
                exit(1);
            }
        }

        inline size_t cleanUp(){
            if(pbbslib::atomic_compare_and_swap(&changing, false, true)){
            c1  = c1 + c2 + c3 -c4 -c5;
            c2 = 0;
            c3 = 0;
            c4 = 0;
            c5 = 0;
            changing = false;
            }
            return c1;
        }

        inline bool operator== (const WTV& j){
            return c1 == j.c1 && c2 == j.c2 && c3 == j.c3 && c4 == j.c4 && c5 == j.c5 ;
        }
        // bool operator == (const WTV& i, const WTV& j){
        //     return i.c1 == j.c1 && i.c2 == j.c2 && i.c3 == j.c3 && i.c4 == j.c4 && i.c5 == j.c5 ;
        // }
    };


    struct vertexHash { //TODO: check
        uint64_t operator ()(const uintE& v) const {return pbbslib::hash64_2(v);}
        int cmp(uintE v, uintE b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}
    };

    struct edgeHash { //TODO: check
        uint64_t operator ()(const EdgeT& v) const{
            return pbbslib::hash_combine(pbbslib::hash64_2(v.first), pbbslib::hash64_2(v.second));}
        int cmp(EdgeT v, EdgeT b) {return (v > b) ? 1 : ((v == b) ? 0 : -1);}

    }; 


    // inline pbbslib::sparse_table<uintE, int, vertexHash> make_vertex_set(size_t m, long space_mult=-1) {
    //     return pbbslib::sparse_table<uintE, int, vertexHash>(m, EMPTYKVB, vertexHash(), space_mult);
    // }
    
}}
