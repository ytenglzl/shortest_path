/**
 * Author: Yue Teng
 * Compile: g++ rgen.cpp -std=c++11 -o rgen
 * 
 * Description: A random generator. It randomly generates streets as input to 
 *              constructor.
 * 
 *              Example of a street:
 *              "street name" (2, 0) (2, 1) (5, 5)
 * 
 *              Command line arguments: 
 *              -s k --- k >= 2. Random number of streets in [2, k].
 *              -n k --- k >= 1. Random number of line segments per street in [1, k].
 *              -l k --- k >= 5. Random waiting time between two outputs in [5, k].
 *              -c k --- k >= 1. Random scope of graph point coordinates in [-k, k].
 *              
 */
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <unistd.h>
#include <random>

using namespace std;

// segmant containint linear coefficients
// range and name of the street
struct seg
{
    float A[2];
    float b;
    float x_range[2];
    float y_range[2];
    float xy0[2]; // end point
    float xy1[2]; // the other end point
    string name;
};

// generates a random int from uniform distr.
int RandUnifInt(int min, int max)
{
    if (max >= min)
    {
        uniform_int_distribution<int> d(min, max);
        random_device rd("/dev/urandom");
        return d(rd);
    }
}

void GnrtSeg(float x0, float y0, float x1, float y1, string name, seg &sg)
{
    sg.x_range[0] = min(x0, x1);
    sg.x_range[1] = max(x0, x1);
    sg.y_range[0] = min(y0, y1);
    sg.y_range[1] = max(y0, y1);
    sg.xy0[0] = x0;
    sg.xy0[1] = y0;
    sg.xy1[0] = x1;
    sg.xy1[1] = y1;
    sg.name = name;

    // compute coefficients
    if (x0 == x1)
    { // x = x0
        sg.A[0] = 1.;
        sg.A[1] = 0.;
        sg.b = x0;
    }
    else if (y0 == y1)
    { // y = y0
        sg.A[0] = 0.;
        sg.A[1] = 1.;
        sg.b = y0;
    }
    else
    {
        float dx = x1 - x0;
        float dy = y1 - y0;
        float ct = dy * x0 - dx * y0;
        sg.A[0] = dy;
        sg.A[1] = -dx;
        sg.b = ct;
    }
}

// calculates the intersection of two
// segments no matter whether existing or not
// but det is !0,
void CalItsc(float *A0, float b0, float *A1, float b1, float *itsc)
{
    float A00, A01, A10, A11, b;

    // find the segment with the first coefficient !0
    // and put it in the first position
    if (A0[0] == 0.)
    {
        A00 = A1[0];
        A01 = A1[1];
        A10 = A0[0];
        A11 = A0[1];
        b = b0;
        b0 = b1;
        b1 = b;
    }
    else
    {
        A00 = A0[0];
        A01 = A0[1];
        A10 = A1[0];
        A11 = A1[1];
    }

    // Gaussian elimination
    if (A10 != 0.)
    {
        float factor;
        factor = -A10 / A00;
        A10 = A10 + A00 * factor;
        A11 = A11 + A01 * factor;
        b1 = b1 + b0 * factor;
    }
    b1 = b1 / A11;
    A11 = A11 / A11;
    b0 = b0 - A01 * b1;
    A01 = A01 - A11 * A01;
    b0 = b0 / A00;
    A00 = A00 / A00;

    // store the intersection
    itsc[0] = b0;
    itsc[1] = b1;
}

// check the status between two segments
// return 0: no touch
// return 1: intersection
// return 2: overlap
int CheckSegStat(seg &sg0, seg &sg1)
{
    // determinant
    float A00 = sg0.A[0];
    float A01 = sg0.A[1];
    float A10 = sg1.A[0];
    float A11 = sg1.A[1];

    float det = A00 * A11 - A01 * A10;
    if (det != 0)
    { // intersection or no touch
        float itsc[2];
        CalItsc(sg0.A, sg0.b, sg1.A, sg1.b, itsc);
        float x = itsc[0];
        float y = itsc[1];
        if (x >= sg0.x_range[0] && x <= sg0.x_range[1] &&
            y >= sg0.y_range[0] && y <= sg0.y_range[1] &&
            x >= sg1.x_range[0] && x <= sg1.x_range[1] &&
            y >= sg1.y_range[0] && y <= sg1.y_range[1])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    { // overlap or no touch
        if (A01 == 0. && A11 == 0.)
        { // A = [[1, 0], [1, 0]]
            if (sg0.b == sg1.b)
            { // on the same line
                if (sg0.y_range[0] > sg1.y_range[1])
                {
                    return 0;
                }
                else if (sg1.y_range[0] > sg0.y_range[1])
                {
                    return 0;
                }
                else
                { // overlap
                    return 2;
                }
            }
            else
            {             // on different lines
                return 0; // no touch
            }
        }
        else if (A00 == 0. && A10 == 0.)
        { // A = [[0, 1], [0, 1]]
            if (sg0.b == sg1.b)
            { // on the same line
                if (sg0.x_range[0] > sg1.x_range[1])
                {
                    return 0;
                }
                else if (sg1.x_range[0] > sg0.x_range[1])
                {
                    return 0;
                }
                else
                {
                    return 2;
                }
            }
            else
            { // on different lines
                return 0;
            }
        }
        else
        { // common parallel
            if (sg0.b == 0 && sg1.b == 0)
            { // on the same line crossing the origin
                if (sg0.x_range[0] > sg1.x_range[1])
                { // no touch
                    return 0;
                }
                else if (sg1.x_range[0] > sg0.x_range[1])
                { // no touch
                    return 0;
                }
                else
                {
                    return 2;
                }
            }
            else if (sg0.b != 0 && sg1.b != 0 && // on the same line
                     A00 / A01 == sg0.b / sg1.b)
            { // not crossing the origin
                if (sg0.x_range[0] > sg1.x_range[1])
                {
                    return 0;
                }
                else if (sg1.x_range[0] > sg0.x_range[1])
                {
                    return 0;
                }
                else
                {
                    return 2;
                }
            }
            else
            { // on different lines
                return 0;
            }
        }
    }
}

// generates a random name
string RandName()
{
    // name length
    int len = RandUnifInt(5, 10);
    // name string
    string repo = "abcdefghijklmnopqrstuvwxyz";
    string name = "";
    for (int i = 0; i < len; ++i)
    {
        int idx = RandUnifInt(0, 25);
        name += repo[idx];
    }

    return name;
}

// generates a random street
// number of segments randomly sampled
// from [1, nk]
// all coordinates randomly sampled from
// [-ck, ck]*[-ck, ck]
// vector<int> streetcoords temporarily stores
// coordinates    of this street to be generated
string RandStreet(int nk, int ck)
{
    vector<int> streetcoords;
    string name = RandName();
    string street = "\"" + name + "\" ";
    int seg_num = RandUnifInt(1, nk);
    int repeat_pt = 0;    // two adjacent points are different
    int sgl_str_stat = 0; // single street status

    do
    { // no touch:0 , intersection: 1, overlap: 2
        repeat_pt = 0;
        sgl_str_stat = 0;

        // generate coordinates
        for (int i = 0; i <= seg_num; ++i)
        {
            int x = RandUnifInt(-ck, ck);
            int y = RandUnifInt(-ck, ck);
            streetcoords.push_back(x);
            streetcoords.push_back(y);
        }

        // generate segments
        vector<seg> street_segs(seg_num);
        for (int i = 0; i < seg_num; ++i)
        {
            float x0 = streetcoords.at(2 * i);
            float y0 = streetcoords.at(2 * i + 1);
            float x1 = streetcoords.at(2 * i + 2);
            float y1 = streetcoords.at(2 * i + 3);
            GnrtSeg(x0, y0, x1, y1, name, street_segs.at(i));
        }

        // check the rationality of this street
        // if there are repeated adjacent points,
        // intersection and overlap
        for (int i = 1; i <= seg_num; ++i)
        {
            if (streetcoords.at(2 * i) == streetcoords.at(2 * i - 2) &&
                streetcoords.at(2 * i + 1) == streetcoords.at(2 * i - 1))
            {
                repeat_pt = 1;
                break;
            }
        }

        if (repeat_pt == 1)
        {
            streetcoords.clear();
        }
        else
        {
            for (int i = 0; i < seg_num - 1; ++i)
            {
                for (int j = i + 1; j < seg_num; ++j)
                {
                    int seg_stat = CheckSegStat(street_segs.at(i), street_segs.at(j));
                    if (seg_stat == 2)
                    {
                        sgl_str_stat == 2;
                        break;
                    }
                    if ((j - i) > 1 && seg_stat == 1)
                    {
                        sgl_str_stat = 1;
                        break;
                    }
                }
                if (sgl_str_stat != 0)
                {
                    streetcoords.clear();
                    break;
                }
            }
        }
    } while (repeat_pt == 1 || sgl_str_stat != 0);

    // string containing name and coordinates
    stringstream coordss;
    for (int i = 0; i <= seg_num; ++i)
    {
        street += "(";
        coordss.str(""); // clear stringstream
        coordss << streetcoords.at(2 * i);
        street += coordss.str();
        street += ",";
        coordss.str("");
        coordss << streetcoords.at(2 * i + 1);
        street += coordss.str();
        street += ") ";
    }

    return street;
}

// check overlap between streets
// return 0: no touch
// return 1: intersection
// return 2: overlap
// return 3: same name
int CheckStrStat(string street0, string street1)
{
    int str_stat = 0;
    stringstream street0ss;
    street0ss << street0;
    stringstream street1ss;
    street1ss << street1;

    // discard street name
    string name0;
    string name1;
    street0ss >> name0;
    street1ss >> name1;
    if (name0.compare(name1) == 0)
    {
        str_stat = 3;
    }
    else
    {
        // parse coordinates
        char lfbr, rgbr, cma;
        float x, y;
        vector<float> street0_coords;
        while (!street0ss.eof())
        {
            street0ss >> lfbr;
            street0ss >> x;
            street0ss >> cma;
            street0ss >> y;
            street0ss >> rgbr;
            if (street0ss.eof())
                break;
            street0_coords.push_back(x);
            street0_coords.push_back(y);
        }

        // store segments
        vector<seg> street0_segs(street0_coords.size() / 2 - 1);
        for (int i = 0; i < street0_coords.size() / 2 - 1; ++i)
        {
            float x0 = street0_coords.at(2 * i);
            float y0 = street0_coords.at(2 * i + 1);
            float x1 = street0_coords.at(2 * i + 2);
            float y1 = street0_coords.at(2 * i + 3);
            GnrtSeg(x0, y0, x1, y1, name0, street0_segs.at(i));
        }
        vector<float> street1_coords;
        while (!street1ss.eof())
        {
            street1ss >> lfbr;
            street1ss >> x;
            street1ss >> cma;
            street1ss >> y;
            street1ss >> rgbr;
            if (street1ss.eof())
                break;
            street1_coords.push_back(x);
            street1_coords.push_back(y);
        }
        vector<seg> street1_segs(street1_coords.size() / 2 - 1);
        for (int i = 0; i < street1_coords.size() / 2 - 1; ++i)
        {
            float x0 = street1_coords.at(2 * i);
            float y0 = street1_coords.at(2 * i + 1);
            float x1 = street1_coords.at(2 * i + 2);
            float y1 = street1_coords.at(2 * i + 3);
            GnrtSeg(x0, y0, x1, y1, name1, street1_segs.at(i));
        }

        // check status of two streets
        for (int i = 0; i < street0_segs.size(); ++i)
        {
            for (int j = 0; j < street1_segs.size(); ++j)
            {
                seg street0seg = street0_segs.at(i);
                seg street1seg = street1_segs.at(j);
                int seg_stat = CheckSegStat(street0seg, street1seg);
                if (seg_stat != 0)
                {
                    if (seg_stat == 1)
                    {
                        str_stat = 1;
                    }
                    else
                    { // seg_stat==2
                        str_stat = 2;
                        break;
                    }
                }
            } // end of for
            if (str_stat == 2)
                break;
        } // end of for
    }

    return str_stat;
}

// input streets for ece650a1.py
// street number randomly sampled from [2, sk]
// streets being strings stored in vector<string>
void StreetsLayout(int sk, int nk, int ck,
                   vector<string> &streets,
                   int max_try = 25)
{
    // random street number
    int street_num = RandUnifInt(2, sk);

    // generate specification of streets
    // maximum: max_try trials otherwise exit all processes
    int trial_num = 0;
    int map_stat = 0;
    do
    {
        streets.clear();
        map_stat = 0;
        for (int i = 0; i < street_num; ++i)
        {
            string street = RandStreet(nk, ck);
            streets.push_back(street);
        } // finished generating map

        // start checking overlaps and intersections
        for (int i = 0; i < street_num - 1; ++i)
        {
            for (int j = i + 1; j < street_num; ++j)
            {
                // provide the status
                int str_stat = CheckStrStat(streets.at(i), streets.at(j));
                if (str_stat == 2)
                {
                    map_stat = 2;
                    break;
                }
                if (str_stat == 3)
                {
                    map_stat = 3;
                    break;
                }
                if (str_stat == 1)
                {
                    map_stat = 1;
                }
            }
            if (map_stat == 2 || map_stat == 3)
            {
                break;
            }
        } // finished checking
    } while (map_stat != 1 && ++trial_num < max_try);

    if (map_stat != 1)
    {
        cerr << "Error: failed to generate valid input for "
             << max_try << " simultaneous attempts" << endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    // check the number of args
    if (argc > 9)
    {
        cerr << "Error: too many options."
             << endl;
        exit(1);
    }
    else
    {
        int option;
        // set default option values for arguments
        int sk = 10, nk = 5, lk = 5, ck = 20;
        // read from command line arguments
        while ((option = getopt(argc, argv, "s:n:l:c:")) != -1)
        {
            switch (option)
            {
            case 's':
            {
                int k;
                istringstream(optarg) >> k;
                if (k >= 2)
                {
                    sk = k;
                }
                else
                {
                    cerr << "Error: s smaller than 2." << endl;
                    exit(1);
                }
                break;
            }
            case 'n':
            {
                int k;
                istringstream(optarg) >> k;
                if (k >= 1)
                {
                    nk = k;
                }
                else
                {
                    cerr << "Error: n smaller than 1." << endl;
                    exit(1);
                }
                break;
            }
            case 'l':
            {
                int k;
                istringstream(optarg) >> k;
                if (k >= 5)
                {
                    lk = k;
                }
                else
                {
                    cerr << "Error: l smaller than 5." << endl;
                    exit(1);
                }
                break;
            }
            case 'c':
            {
                int k;
                istringstream(optarg) >> k;
                if (k >= 1)
                {
                    ck = k;
                }
                else
                {
                    cerr << "Error: c smaller than 1." << endl;
                    exit(1);
                }
                break;
            }
            default:
            {
                cerr << "Error: wrong argument option." << endl;
                exit(1);
                break;
            }
            }
        } // successfully modified options

        // start outputing lines for ece650a1.py
        int slp_time = RandUnifInt(5, lk); // sleep time
        int max_try = 25;
        vector<string> streets;
        while (1)
        {
            streets.clear();
            StreetsLayout(sk, nk, ck, streets, max_try);
            string otpt = "r\n";
            for (int i = 0; i < streets.size(); ++i)
            {
                otpt = otpt + "a " + streets.at(i) + "\n";
            }
            otpt += "g\n";         // cout<<endl == cout<<”\n” <<flush
            cout << otpt << flush; // fresh buffer
            sleep(slp_time);
        }
    }

    return 0;
}
