/**
 * Author: Yue Teng
 * Compile: g++ path_gener.cpp Vertex/Vertex.cpp -std=c++11 -o path_gener
 * 
 * Description: Find the shortest path between two nodes in an
 *              unweighted, undirected, graph. The input is the 
 *              output graph from constructor. The command to 
 *              invoke an output of shortest path:
 * 
 *              s 2 4
 * 
 *              after which a shortest path will output, like
 * 
 *              2-0-4
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <queue>
#include <algorithm>

using namespace std;

// Vertex class
class Vertex{
 public:
  // label
  int label_;
  //0 white, 1 grey, 2 black
  int color_;
  //distance
  int distance_;
  //pointer to predecessor
  Vertex *parent_;

  Vertex();
};

// Constructor for Vertex.
// Initialization cannot be specified
Vertex::Vertex() {
  // label specifying no vertex
  label_ = -1;
  // white
  color_ = 0;
  // -1 means infinity
  distance_ = -1;
  // no predecessor
  parent_ = NULL;
}

// generates a pointer to an adjacency
// list which is a vector containing
// vector pointers to a vector<int>
// specitying adjacent vertices
vector<vector<int> *> *GnrtAdj(int number, vector<int> edgepoint)
{
    // array of pointers to vectors for each vertex
    vector<vector<int> *> *vecarr = new vector<vector<int> *>(number);
    // vector<int> *vecarr[number];

    // initialize array with NULL pointers
    for (int i = 0; i < number; i++)
    {
        vecarr->at(i) = NULL;
    }

    // discard repeated end pointes of edges
    set<int> segpt(edgepoint.begin(), edgepoint.end());

    // construct a vector for each vertex as edge end point
    for (set<int>::iterator it = segpt.begin(); it != segpt.end(); it++)
    {
        // adjacency points
        vector<int> *adjp = new vector<int>;
        vecarr->at(*it) = adjp;
    }

    // connect the edge end
    // points to adjacent ones
    for (int i = 0; i < (edgepoint.size() / 2); i++)
    {
        int x = edgepoint[2 * i];
        int y = edgepoint[2 * i + 1];

        // push back if the element pushed
        // does not exist in the vector
        vector<int>::iterator it;
        vector<int> *p2vecx = vecarr->at(x);
        vector<int> *p2vecy = vecarr->at(y);
        it = find(p2vecx->begin(), p2vecx->end(), y);
        if (it == p2vecx->end())
        {
            p2vecx->push_back(y);
        }
        it = find(p2vecy->begin(), p2vecy->end(), x);
        if (it == p2vecy->end())
        {
            p2vecy->push_back(x);
        }
    }

    // return the pointer of
    // the array containing
    // pointer of vectors which
    // are adjacent points
    return vecarr;
}

// frees the allocated space
// for adjacency list in GnrtAdj
vector<vector<int> *> *FrVecArr(vector<vector<int> *> *vecarr)
{
    if (vecarr != NULL)
    {
        // free vectors for each
        // existing vertex in the array
        for (int i = 0; i < (vecarr->size()); i++)
        {
            if (vecarr->at(i) != NULL)
            {
                delete vecarr->at(i);
            }
        }
        // free vector array
        delete vecarr;
    }

    return NULL;
}

// generates a pointer to a vector
// containing pointers to Vertex
vector<Vertex *> *GnrtVtxArr(int number)
{
    vector<Vertex *> *vtxarr = new vector<Vertex *>(number);
    for (int i = 0; i < number; i++)
    {
        vtxarr->at(i) = new Vertex;
        vtxarr->at(i)->label_ = i;
    }

    return vtxarr;
}

// frees the allocated Vertex array in GnrtVtxArr
vector<Vertex *> *FrVtxArr(vector<Vertex *> *vtxarr)
{
    if (vtxarr != NULL)
    {
        // free the vertices
        for (int i = 0; i < (vtxarr->size()); i++)
        {
            delete vtxarr->at(i);
        }
        delete vtxarr;
    }

    return NULL;
}

// breadth-first search in Vertex array
// given adjacency list and source
void Bfs(vector<vector<int> *> *vecarr, vector<Vertex *> *vtxarr, int source)
{
    if (vecarr->at(source) != NULL && source < (vtxarr->size()))
    {
        Vertex *p2vtxscr = vtxarr->at(source);
        p2vtxscr->color_ = 1; // Grey
        p2vtxscr->distance_ = 0;

        queue<int> que;
        que.push(source);

        while (!que.empty())
        {
            // label of the first vertex in queue
            int u = que.front();
            Vertex *p2vtxu = vtxarr->at(u);

            // dequeue the first vertex label
            que.pop();
            int vecsize = (vecarr->at(u))->size();
            for (int i = 0; i < vecsize; i++)
            {
                // label for the vth vertex
                // stored in the ith position in vector
                int v = (vecarr->at(u))->at(i);
                Vertex *p2vtxv = vtxarr->at(v);
                if (p2vtxv->color_ == 0)
                {
                    p2vtxv->color_ = 1;
                    p2vtxv->distance_ = p2vtxu->distance_ + 1;
                    p2vtxv->parent_ = p2vtxu;
                    que.push(v);
                }
            }
            p2vtxu->color_ = 2;
        }
    }
}

// prints shortest path
// given source and target
void PrtPath(vector<Vertex *> *vtxarr, vector<int> *ppath, int source, int target)
{
    Vertex *p2vtxtg = vtxarr->at(target);
    if (source < (vtxarr->size()))
    {
        if (target == source)
        {
            ppath->push_back(source);
        }
        else if (p2vtxtg->parent_ == NULL)
        {
            cerr << "Error: No path exists.";
            exit(1);
        }
        else
        {
            PrtPath(vtxarr, ppath, source, (p2vtxtg->parent_)->label_);
            ppath->push_back(target);
        }
    }
}

// restores modified Vertex array
void Restore(vector<Vertex *> *vtxarr)
{
    for (int i = 0; i < (vtxarr->size()); i++)
    {
        Vertex *p2vtxi = vtxarr->at(i);
        p2vtxi->color_ = 0;
        p2vtxi->distance_ = -1;
        p2vtxi->parent_ = NULL;
    }
}

int main()
{
    // number of vertices
    int num = 0;

    // container for edge points
    vector<int> egpt;
    vector<vector<int> *> *adj = NULL;
    vector<Vertex *> *vtcs = NULL;

    // stores the input lines
    string lines[2] = {"", ""};
    
    // read from stdin until EOF
    while (!cin.eof())
    {
        // read a line of input until
        // EOL and store in a string
        string line;
        getline(cin, line);
        if (cin.eof())
        {
            break;
        }

        // create an input stream
        // based on the line
        // we will use the input
        // stream to parse the line
        istringstream input(line);

        if (!input.eof())
        {
            // parse the command
            char cmd;
            input >> cmd;
            if (input.fail())
            {
                cerr << "Error: Cannot read a command." << endl;
                exit(1);
            }
            else if (cmd == 'V')
            {
                int oldnum = num;
                lines[0] = "";
                lines[1] = "";
                num = 0;
                input >> num;
                if (num > 0)
                {
                    egpt.clear();
                    adj = FrVecArr(adj);
                    vtcs = FrVtxArr(vtcs);
                    vtcs = GnrtVtxArr(num);
                    lines[0] = line;
                }
                else
                {
                    num = oldnum;
                    cerr << "Error: Number following V "
                         << "must be bigger than 0." << endl;
                    exit(1);
                }
            }

            else if (cmd == 'E')
            {
                if (num == 0)
                {
                    cerr << "Error: Input 'V number (>0)' first." << endl;
                    exit(1);
                }
                else
                {
                    egpt.clear();
                    adj = FrVecArr(adj);
                    Restore(vtcs); // restoration
                    lines[1] = "";
                    char lfbr;  // left brace
                    char lfbra; // left bracket
                    char rgbra; // right bracket
                    char cm;    // comma
                    int end1;   // end point 1
                    int end2;   // end point 2
                    input >> lfbr;
                    do
                    {
                        input >> lfbra;
                        input >> end1 >> cm >> end2;
                        input >> rgbra;
                        if (end1 < num &&
                            end1 >= 0 &&
                            end2 < num &&
                            end2 >= 0)
                        {
                            egpt.push_back(end1);
                            egpt.push_back(end2);
                        }
                        else
                        {
                            cerr << "Error: Edge vertices out of range." << endl;
                            exit(1);
                            egpt.clear();
                            break;
                        }
                        input >> cm;
                    } while (cm == ',');
                    adj = GnrtAdj(num, egpt);

                    // output the input lines
                    lines[1] = line;
                    if (lines[0] != "" && lines[1] != "")
                    {
                        cout << lines[0] << endl;
                        cout << lines[1] << endl;
                    }
                }
            }

            else if (cmd == 's')
            {
                if (num == 0 || egpt.empty())
                {
                    cerr << "Error: No 'V number' "
                         << "or 'E {<edge>,...}' input." << endl;
                    exit(1);
                }
                else
                {
                    Restore(vtcs); // restoration
                    int st;        // starting point
                    int ed;        // end point
                    input >> st >> ed;
                    if (st < (vtcs->size()) && ed < (vtcs->size()))
                    {
                        Bfs(adj, vtcs, st);
                        vector<int> path;
                        PrtPath(vtcs, &path, st, ed);
                        for (int i = 0; i < path.size(); i++)
                        {
                            cout << path[i];
                            if (i != path.size() - 1)
                            {
                                cout << "-";
                            }
                        }
                        cout << endl;
                    }
                    else
                    {
                        cerr << "Error: Starting or "
                             << "ending point out of range." << endl;
                        exit(1);
                    }
                }
            }
            else
            {
                cerr << "Error: Wrong command." << endl;
                exit(1);
            }
        }
    }

    return 0;
}
