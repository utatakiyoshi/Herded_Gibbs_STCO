#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <cassert>
using namespace std;

uint32_t draw_random_int(void) {
    return rand();
}
double draw_random_double(){
    return (double)rand()/RAND_MAX;
}

int main(int argc, char** argv){
    ifstream file_list;
    ifstream fin;
    ofstream fout;
    int H,W;
    string fin_name, fout_name;
    double p;
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-filelist")==0){i++; file_list.open(argv[i]);}
        if(strcmp(argv[i], "-p")==0){i++; p=atof(argv[i]);}
    }
    int N;
    file_list >> N;
    for(int i=0;i<N;i++){
        file_list >> H;
        file_list >> W;
        file_list >> fout_name;
        file_list >> fin_name;
        int HW=H * W;
        int original;
        fin.open(fin_name);
        fout.open(fout_name);
        for(int y=0;y<H;y++){
            for(int x=0;x<W;x++){
                fin >> original;
                if(draw_random_double() < p){
                    original = (original + 1) % 2;
                }
                fout << original << " ";
            }
            fout << endl;
        }
        fin.close();
        fout.close();
    }
}
