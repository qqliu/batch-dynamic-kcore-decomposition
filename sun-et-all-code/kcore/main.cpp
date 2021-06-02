#include <fstream>
#include <cassert>
#include <ctime>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "kcore.h"

using namespace std;

int main(int argc,char* argv[])
{

    bool batch = true;
    bool deletes = true;
    int space = 0;
    if (argc == 6)
        space = atoi(argv[5]);

    /*if(argc != 4) {
        fprintf(stderr,"Usage: %s <graph-file> <edge-file> <threads>\n",argv[0]);
        return -1;
    }*/
    Decomposition dec(atoi(argv[2]));
    if (!batch){
        cout << "started loading file" << endl;
        if(!dec.Load(argv[1],argv[2])) {
            fprintf(stderr,"Load file failed\n");
            return -1;
        }
        //cout << "started deletions"<< endl;
        //dec.Delete();
        //cout << "started insertions" << endl;
        //dec.Insert();
    } else {
        if (!deletes) {
            const char *graphfile = argv[1];
            ifstream update_edges;
            update_edges.open(graphfile, ios_base::in);
            if (update_edges.fail()) {
                return -1;
            }

            string content;
            content.assign((istreambuf_iterator<char>(update_edges)),
                    (istreambuf_iterator<char>()));
            stringstream ss(content);
            istream_iterator<string> begin(ss);
            istream_iterator<string> end;
            vector<string> vstrings(begin, end);
            assert(!(vstrings.size() & 1));
            assert(!(vstrings.size() == 0));

            int batch_size = atoi(argv[3]);

            for (int i = 0; i < vstrings.size() - 1; i += 2 * batch_size) {
                //cout << "iteration: " << i << endl;
                auto load_batch_time = dec.LoadBatch(vstrings, i, batch_size);
                dec.Insert(load_batch_time);
                //dec.Delete();
                if (space == 1)
                    cout << "Space: " << dec.GetSize() << endl;
            }
            update_edges.close();
        } else {
            if(!dec.LoadStatic(argv[1])) {
                fprintf(stderr,"Load file failed\n");
                return -1;
            }

            const char *updatefile = argv[3];
            ifstream update_edges;
            update_edges.open(updatefile, ios_base::in);
            if (update_edges.fail()) {
                return -1;
            }
            string content;
            content.assign((istreambuf_iterator<char>(update_edges)),
                    (istreambuf_iterator<char>()));
            stringstream ss(content);
            istream_iterator<string> begin(ss);
            istream_iterator<string> end;
            vector<string> vstrings(begin, end);
            assert(!(vstrings.size() & 1));
            assert(!(vstrings.size() == 0));

            int batch_size = atoi(argv[4]);

            //Timer total;
            for (int i = 0; i < vstrings.size() - 1; i += 2 * batch_size) {
                //total.reset();
                //cout << "iteration: " << i << endl;
                auto load_batch_time = dec.LoadBatchDeletes(vstrings, i, batch_size);
                //dec.Insert(load_batch_time);
                //cout << "LOADING TIME: " << total.elapsed() << endl;
                dec.Delete(load_batch_time);
                //cout << "TOTAL TIME: " << total.elapsed() << endl;
                if (space == 1)
                    cout << "Space: " << dec.GetSize() << endl;
            }
            update_edges.close();
        }
        /*std::cout.flush();
        std::ofstream out("new_all_inserts");
        std::cout.rdbuf(out.rdbuf());
        dec.dumpCores();*/
    }
}
