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
#include <boost/math/common_factor_rt.hpp>
using namespace std;

uint32_t draw_random_int(void) {
    return rand();
}
double draw_random_double(){
    return (double)rand()/RAND_MAX;
}

class settings{
  public: 
    string input_file_name;
    string observe_time_file_name;
    string output_header_file_name;
    string output_file_name;
    
    ifstream input_file;
    ifstream observe_time_file;
    ofstream output_header_file;
    ofstream output_file;
    
    int observe_time_num;
    vector<int> observe_time;
    int tmax;
    
    int H;
    int W;
    int HW;
    double B;
    double J;
    
    int num_trial;
    int input_num;
    string expe_name;
    
    int algorithm;
    int postprocess;
    int rng_length;
};

enum e_algorithm{
    HERDING,
    HERDING_SHARE,
    HERDING_EX_SHARE,
    GIBBS,
    QMCMC
};

enum e_postprocess{
    MEAN,
    SAMPLE
};

class variable{
  public:
    int spin;
    vector<double> weights;
    void init(settings& s, int weights_num){
        weights.resize(weights_num);
    }
    void reset(settings& s){
        for(int i=0;i<weights.size();i++){
            weights[i] = draw_random_double();
        }
        spin = draw_random_int() % 2;
    }
    void update_herding(double p, int state){
        weights[state] += p;
        if(weights[state] >= 1) { spin = 1; weights[state] -= 1; }
        else { spin = 0; }
    }
    void update_random(double p, double u){
        if(u<p){ spin = 1;}
        else { spin = 0; }
    }
};

class graph{
  public:
    vector<variable> vs;
    vector<int> y;
    vector<int> z;
    vector<int> order;
    vector<vector<int> > ns;
    vector<double> u_seq;
    int u_seq_head;
    int idx(settings& s, int x, int y){
        if(x < 0 || x >= s.W || y < 0 || y >= s.H){ return -1; }
        return x+s.W*y;
    }
    void init_graph(settings& s){
    }
    
    void load_graph(settings& s, ifstream& img_file, ifstream& answer_img_file){
        y.resize(s.HW);
        z.resize(s.HW);
        
        order.resize(s.HW);
        for(int i=0;i<s.HW;i++){
            order[i] = i;
        }
        random_shuffle(order.begin(), order.end());

        for(int i=0;i<s.HW;i++){
            img_file >> y[order[i]];
            answer_img_file >> z[order[i]];
        }
        vs.resize(s.HW);
        ns.resize(s.HW);
        for(int y=0;y<s.H;y++){
            for(int x=0;x<s.W;x++){
                int i = order[idx(s, x, y)];
                ns[i].clear();
                int n1 = idx(s, x-1, y);
                int n2 = idx(s, x+1, y);
                int n3 = idx(s, x, y-1);
                int n4 = idx(s, x, y+1);
                if(n1 >= 0){ ns[i].push_back(order[n1]); }
                if(n2 >= 0){ ns[i].push_back(order[n2]); }
                if(n3 >= 0){ ns[i].push_back(order[n3]); }
                if(n4 >= 0){ ns[i].push_back(order[n4]); }
            }
        }
        int weights_num;
        for(int i=0;i<s.HW;i++){
            switch(s.algorithm){
                case HERDING:
                    weights_num = 1 << (ns[i].size());
                    break;
                case HERDING_SHARE:
                    weights_num = 1 + (ns[i].size());
                    break;
                case HERDING_EX_SHARE:
                    weights_num = 1;
                    break;
                case GIBBS:
                case QMCMC:
                    weights_num = 0;
                    break;
            }
            vs[i].init(s, weights_num);
        }
    }
    void gen_u_seq(settings &s){
        //const int modulo = 15;
        //const int hoge[] = {1, 8, 4, 2, 9, 12, 6, 11, 5, 10, 13, 14, 15, 7, 3};
        const int modulo = 31;
        const int hoge[] = {1, 16, 8, 4, 18, 9, 20, 26, 13, 6, 19, 25, 28, 30, 31, 15, 7, 3, 17, 24, 12, 22, 27, 29, 14, 23, 11, 21, 10, 5, 2};
        // const int modulo = 63;
        // const int hoge[] = {1, 32, 16, 8, 4, 2, 33, 48, 24, 12, 6, 35, 17, 40, 20, 10, 37, 50, 57, 60, 30, 47, 23, 11, 5, 34, 49, 56, 28, 14, 39, 19, 9, 36, 18, 41, 52, 26, 45, 54, 59, 29, 46, 55, 27, 13, 38, 51, 25, 44, 22, 43, 21, 42, 53, 58, 61, 62, 63, 31, 15, 7, 3}
        s.rng_length = s.HW * modulo;
        //        assert(s.tmax == modulo);
        int b = boost::math::gcd(s.HW, modulo);
        int g = 3;
        assert(boost::math::gcd(g, modulo) == 1);
        //cerr<<b<<endl;
        u_seq.resize(s.rng_length);
        vector<double> offset(s.HW);
        for(int i=0;i<s.HW;i++){
            offset[i] = draw_random_double();
        }
        for(int j=0;j<b;j++){
            for(int i=0;i<s.rng_length/b;i++){
                int x = i+j*s.rng_length/b;
                int y = (i+j)%modulo;
                double u = (double)hoge[(y*g)%modulo]/modulo;
                u += offset[x%(s.HW)];
                if(u>=1.0){u-=1.0;}
                u_seq[x] = u;
            }
        }
        u_seq_head = 0;
    }

    void reset_graph(settings& s){
        for(int i=0;i<s.HW;i++){
            vs[i].reset(s);
            vs[i].spin = y[i];
        }
        if(s.algorithm == QMCMC){
            gen_u_seq(s);
        }
    }
    
    void run(settings& s, vector<double>& vec_rec_err){
        vector<int> sum(s.HW, 0);
        int observe_time_head = 0;
        for(int t=0;t<s.tmax;t++){
            for(int i=0;i<s.HW;i++){
                double p;
                double E = 0;
                E += -s.B * (y[i]*2-1);
                for(int j=0;j<ns[i].size();j++){
                    int spin_j = vs[ns[i][j]].spin;
                    E += -s.J * (spin_j*2-1);
                }
                p = 1.0/(1.0+exp(E));
                
                if(s.algorithm == HERDING || s.algorithm == HERDING_SHARE || s.algorithm == HERDING_EX_SHARE){
                    int state = 0;
                    for(int j=0;j<ns[i].size();j++){
                        int spin_j = vs[ns[i][j]].spin;
                        if(s.algorithm == HERDING){ state |= (spin_j << j); }
                        else if(s.algorithm == HERDING_SHARE){ state += spin_j; }
                        else { state = 0; }
                    }
                    vs[i].update_herding(p, state);
                } else {
                    double u;
                    if(s.algorithm == GIBBS){ u = draw_random_double(); }
                    else { u = u_seq[u_seq_head]; u_seq_head++; if(u_seq_head>=s.rng_length){u_seq_head=0;}}
                    vs[i].update_random(p, u);
                }
                
                sum[i] += vs[i].spin;
            }
            
            if(t+1 == s.observe_time[observe_time_head]){
                double rec_err = 0;
                for(int i=0;i<s.HW;i++){
                    int s1,s0;
                    if(s.postprocess==MEAN){
                        s1 = sum[i];
                        s0 = (t+1) - s1;
                    } else {
                        s1 = vs[i].spin;
                        s0 = 1-s1;
                    }
                    if(s0>s1 && z[i]==1){ rec_err += 1.0; }
                    if(s0==s1){  rec_err += 0.5; }
                    if(s0<s1 && z[i]==0){ rec_err += 1.0; }
                }
                vec_rec_err[observe_time_head] = rec_err;// / s.HW;
                observe_time_head++;
            }
        }
    }
};


int main(const int argc, const char** argv){
    settings s;
    s.output_header_file_name = "out_settings.txt";
    s.output_file_name = "result1.csv";
    s.observe_time_file_name = "observe_times.txt";
    s.input_file_name="input_file_list.txt";
    s.algorithm = GIBBS;
    s.postprocess = MEAN;
    for(int i=1;i<argc;i++){
//        if(strcmp(argv[i], "-num_trial")==0){i++; s.num_trial=atoi(argv[i]);}
//        if(strcmp(argv[i], "-H")==0){i++; s.H=atoi(argv[i]);}
//        if(strcmp(argv[i], "-W")==0){i++; s.W=atoi(argv[i]);}
//        if(strcmp(argv[i], "-height")==0){i++; s.H=atoi(argv[i]);}
//        if(strcmp(argv[i], "-width")==0){i++; s.W=atoi(argv[i]);}
//        if(strcmp(argv[i], "-B")==0){i++; s.B=atof(argv[i]);}
//        if(strcmp(argv[i], "-J")==0){i++; s.J=atof(argv[i]);}
        if(strcmp(argv[i], "-expe_name")==0){i++; s.expe_name = argv[i];}
        
        if(strcmp(argv[i], "-HG")==0){s.algorithm=HERDING;}
        if(strcmp(argv[i], "-sHG")==0){s.algorithm=HERDING_SHARE;}
        if(strcmp(argv[i], "-SHG")==0){s.algorithm=HERDING_SHARE;}
        if(strcmp(argv[i], "-exsHG")==0){s.algorithm=HERDING_EX_SHARE;}
        if(strcmp(argv[i], "-exSHG")==0){s.algorithm=HERDING_EX_SHARE;}
        if(strcmp(argv[i], "-Gibbs")==0){s.algorithm=GIBBS;}
        if(strcmp(argv[i], "-gibbs")==0){s.algorithm=GIBBS;}
        if(strcmp(argv[i], "-QMCMC")==0){s.algorithm=QMCMC;}
        if(strcmp(argv[i], "-qmcmc")==0){s.algorithm=QMCMC;}
        
        if(strcmp(argv[i], "-mean")==0){s.postprocess=MEAN;}
        if(strcmp(argv[i], "-sample")==0){s.postprocess=SAMPLE;}
        
        if(strcmp(argv[i], "-output_header_file")==0){i++; s.output_header_file_name=argv[i];}
        if(strcmp(argv[i], "-output_file")==0){i++; s.output_file_name=argv[i];}
        if(strcmp(argv[i], "-observe_time_file")==0){i++; s.observe_time_file_name=argv[i];}
        if(strcmp(argv[i], "-input_file")==0){i++; s.input_file_name=argv[i];}
    }
    
    //================
    
    s.observe_time_file.open(s.observe_time_file_name);
    
    int ot;
    char dummy;
    s.observe_time_file >> s.observe_time_num;
    s.observe_time.resize(s.observe_time_num);
    for(int i=0; i<s.observe_time_num; i++){
        s.observe_time_file >> ot >> dummy;
        s.observe_time[i] = ot;
    }
    s.tmax=s.observe_time[s.observe_time_num-1];
    s.observe_time_file.close();
    
    //================
    
    s.input_file.open(s.input_file_name);    
    
    s.input_file >> s.B;
    s.input_file >> s.J;
    s.input_file >> s.num_trial;
    cerr<<"num_trial: "<<s.num_trial<<endl;
    
    vector<double> rec_err(s.observe_time_num,0);
    vector<double> rec_err_sum(s.observe_time_num*s.num_trial,0);
    vector<double> rec_err_ssum(s.observe_time_num*s.num_trial,0);
    vector<double> rec_err_mean(s.observe_time_num,0);
    vector<double> rec_err_var(s.observe_time_num,0);
    vector<double> rec_err_std(s.observe_time_num,0);
    
    graph g;
    g.init_graph(s);
            
    int num_run=0;
    ifstream img_file, answer_img_file;
    string img_file_name, answer_img_file_name;
    s.input_file >> s.input_num;
    cerr<<"input_num: "<<s.input_num<<endl;
    double num_pixels = 0;
    for(int img_i = 0; img_i < s.input_num; img_i++){
        cerr<<endl;
        cerr<< img_i << " ";
        s.input_file >> s.H;
        s.input_file >> s.W;
        s.input_file >> img_file_name;
        s.input_file >> answer_img_file_name;
        cerr<< s.H <<" " << s.W << " ";
        s.HW = s.H * s.W;
        num_pixels += s.HW;
        img_file.open(img_file_name);
        answer_img_file.open(answer_img_file_name);
        g.load_graph(s, img_file, answer_img_file);
        img_file.close();
        answer_img_file.close();
        
        for(int k=0;k<s.num_trial;k++){
            cerr<< k << " ";
            g.reset_graph(s);
            g.run(s,rec_err);
            for(int i=0;i<s.observe_time_num;i++){
                rec_err_sum[i*s.num_trial + k] += rec_err[i];
            }
            num_run++;
        }
    }
    s.input_file.close();
    
    //================
    s.output_header_file.open(s.output_header_file_name, std::ofstream::app);
    s.output_file.open(s.output_file_name, std::ofstream::app);
    
    for(int k=0;k<s.num_trial;k++){
        for(int i=0;i<s.observe_time_num;i++){
            rec_err_sum[i*s.num_trial + k] /= num_pixels;
            rec_err_mean[i] += rec_err_sum[i*s.num_trial + k];
            rec_err_ssum[i] += rec_err_sum[i*s.num_trial + k] * rec_err_sum[i*s.num_trial + k];
        }
    }
    for(int i=0;i<s.observe_time_num;i++){
        rec_err_mean[i] /= s.num_trial;
        rec_err_var[i] = rec_err_ssum[i] / s.num_trial - rec_err_mean[i] * rec_err_mean[i];
        rec_err_std[i] = sqrt(rec_err_var[i]);
    }
    s.output_header_file << s.expe_name << "$mean ";
    for(double mean : rec_err_mean){ s.output_file << mean << ", "; } 
    s.output_file << endl;
    s.output_header_file << s.expe_name << "$std ";
    for(double std : rec_err_std){ s.output_file << std << ", ";  } 
    s.output_file << endl;
    
    s.output_file.close();
    s.output_header_file.close();
}
