#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <climits>
using namespace std;

uint32_t draw_random_int(void) { 
    return rand();
}
double draw_random_double(){
    return (double)rand()/RAND_MAX;
}
double draw_random_gaussian(double mu, double sigma){
    double alpha = draw_random_double();
    double beta = draw_random_double();
    return sqrt(-2.0*log(alpha))*sin(2*M_PI*beta)*sigma+mu;
}

static const int N=9;

struct variable;
variable* v[N];
vector<double> delta_s;

struct herding;
struct herding{
    double p;
    int algorithm;
    double w;
    
    void reset(double w_){
        w=w_;
    }
    void init(double p_, int algorithm_){
        p=p_; algorithm=algorithm_;
    }
    int ask(int& rng_drawed){
        int ret;
        if(algorithm==0){
            if(p>0){ rng_drawed=1; }
            ret = ((draw_random_double() < p) ? 1:0);
        } else if (algorithm==1) {
            ret = ((w < p) ? 1:0);
            w+=p; if(w>1.0){w-=1.0;}
        } 
        return ret;
    }
};
struct cond{
    double p;
    int type, algorithm; // [algorithm]  0: Random 1:  Herded
    herding* weight_variable;
    void reset(){
        weight_variable->reset(draw_random_double());
    }
    void init(int algorithm_, double p_){
        algorithm = algorithm_; p=p_;
        if(weight_variable==NULL){ weight_variable = (herding*)malloc(sizeof(herding)); }
        double p_new = ((p<=0.5) ? p : (1-p));
        type = ((p<=0.5) ? 0 : 1);
        weight_variable->init(p_new, algorithm);
    }
    int ask(int& rng_drawed){
        int ret;
        switch(algorithm){
            case 1: // (Incomplete) Herded
                ret = (type + weight_variable->ask(rng_drawed))%2; break;
            default: // Random
                rng_drawed=1; ret = ((draw_random_double() < p) ? 1 : 0); break;
        } 
        return ret;
    }
};
struct variable{
    static const int MAXNN = 20;
    static const int MAXD = 1024*1024;
    int spin, nn_size, d, algorithm;
    int nn[MAXNN];
    double K[MAXNN];
    double H;
    vector<cond> c, master_s;
    vector<int> group_id, type_s;
    vector<double> p_s;
    
    double** dp;
    int** dq;
    double* x;
    vector<pair<double, int> > vec;

    double get_E(){
        double E = spin==0?H:-H;
        for(int i=0;i<nn_size;i++){
            if(spin == v[nn[i]]->spin){
                E -= K[i]/2;
            } else {
                E += K[i]/2;
            }
        }
        return E;
    }

    void init(double H_){
        H = H_; nn_size = 0;
    }
    
    int check_add_nn(int nn_new){
        if(nn_size >= MAXNN){return 0;}
        for(int i=0;i<nn_size;i++){
            if(nn[i] == nn_new){return 0;}
        }
        return 1;
    }
    
    void add_nn(int nn_new, double K_){
        nn[nn_size] = nn_new;
        K[nn_size] = K_;
        nn_size++;
    }
    
    void reset(){
        for(int i=0;i<d;i++){
            c[i].reset();
        }
    }
    
    void init_final(int algorithm_){
        spin = draw_random_int()%2;
        algorithm = algorithm_; 
        d=(1<<nn_size);
        c.resize(d);
        p_s.resize(d);
        for(int i=0;i<d;i++){
            double H0 = -H, H1 = H;
            for(int k = 0; k < nn_size; k++){
                H0 += (((i & (1<<k))==0) ? +K[k] : -K[k]);
                H1 += (((i & (1<<k))==0) ? -K[k] : +K[k]);
            }
            double p = exp(H0), q = exp(H1);
            p_s[i]= q/(p+q);
        }
        
        for(int i=0;i<d;i++){
            c[i].init(algorithm, p_s[i]);
        }
    }
    
    int next_state(){
        int state=0;
        int rng_drawed = 0;
        for(int i=0;i<nn_size;i++){
            state |= ((v[nn[i]]->spin) << i);
        }
        
        spin = c[state].ask(rng_drawed);
        return rng_drawed;
    }
};
int add_nn_sub(int i, int j, double K){
    if(i==j){return 0;}
    if(v[i]->check_add_nn(j)==0){return 0;}
    if(v[j]->check_add_nn(i)==0){return 0;}
    v[i]->add_nn(j, K);
    v[j]->add_nn(i, K);
    return 1;
}

void modify_graph(int step){
    for(int i=0;i<N;i++){
        add_nn_sub(i, (i+step)%N, 0);
    }
}
void build_graph(){
    double H=0.2;
    double K=0.5;
    int head=0;
    double x;
    for(int i=0;i<N;i++){
        v[i]=(variable*)malloc(sizeof(variable));
        x = draw_random_gaussian(H, 0.05);
        v[i]->init(x);
    }
    
    for(int i=0;i<N;i++){
        x = draw_random_gaussian(K, 0.05);
        add_nn_sub(i, (i+1)%N, x);
    }
}
void init_graph(int algorithm){
    for(int i=0;i<N;i++){
        v[i]->init_final(algorithm);
    }
}

void reset_graph(){
    for(int i=0;i<N;i++){
        v[i]->reset();
    }
}


double run(int maxT, int* measureT, double* mu_out, double* sigma_out){
    reset_graph();

    double mu=0, sigma=0;
    
    int head = 0;
    
    for(int t=0; t <= maxT; t++){
        if(t==measureT[head]){
            mu_out[head] = mu/t;
            sigma_out[head] = sigma/t;
            head++;
        }
        int drawed_cnt=0;
        for(int i=0;i<N;i++){
            drawed_cnt += v[i]->next_state();
        }
        double x=0;
        for(int i=0;i<N;i++){
            x+=v[i]->spin;
        }
        mu+=x; sigma+=x*x;
    }
}

double calc_exact(double* mu_out, double* sigma_out){
    int M = (1<<N);
    double Z=0;
    double mu=0;
    double sigma=0;
    for(int i=0;i<M;i++){
        int S=0;
        double E=0;
        for(int j=0;j<N;j++){
            v[j]->spin = ((i>>j)&1);
            S+=v[j]->spin;
        }
        for(int j=0;j<N;j++){
            E+=v[j]->get_E();
        }
        double p = exp(-E);
        fprintf(stderr,"%lf %lf\n", p, E);

        Z+=p;
        mu+=p*S;
        sigma+=p*S*S;
    }
    *mu_out=mu/Z;
    *sigma_out=sigma/Z;
}

int hoge(int num_measure, int* measureT, int algorithm, double mu_approx, double sigma_approx){
    double mu_s[1000][50];
    double sigma_s[1000][50];
    int num_trial = 400;
    
    init_graph(algorithm);

    int maxT = measureT[num_measure-1];
    for(int trial=0;trial<num_trial;trial++){
        run(maxT, measureT, mu_s[trial], sigma_s[trial]);
        for(int i=0;i<num_measure; i++){
            mu_s[trial][i] = abs(mu_s[trial][i]-mu_approx);
            sigma_s[trial][i] = abs(sigma_s[trial][i]-sigma_approx);
        }
    }
    vector<double> mu_vec(num_trial);
    vector<double> sigma_vec(num_trial);
    int lower = num_trial*0.2;
    int middle = num_trial*0.5;
    int upper = num_trial*0.8;
    
    for(int i=0;i<num_measure; i++){
        for(int j=0;j<num_trial;j++){
            mu_vec[j]=mu_s[j][i];
            sigma_vec[j]=sigma_s[j][i];
        }
        sort(mu_vec.begin(), mu_vec.end());
        sort(sigma_vec.begin(), sigma_vec.end());
        printf("%5lf",mu_vec[middle]);
        printf("\n");
    }
    fflush(stdout);
}

int main(){
    srand((unsigned)time(NULL));
    
    int num_measure = 22;
    int measureT[50];
    
    measureT[0] = 2;
    for(int head = 1; head < num_measure; head++){
        measureT[head] = measureT[head-1] * 2;
    }
    
    build_graph();

    double mu_approx;
    double sigma_approx;
    int initial_run_maxT = 1024*1024*128;

    init_graph(0);

    calc_exact(&mu_approx, &sigma_approx);
    fprintf(stderr, "%lf %lf\n", mu_approx, sigma_approx);
    hoge(num_measure, measureT, 0, mu_approx, sigma_approx);
    hoge(num_measure, measureT, 1, mu_approx, sigma_approx);
    
    int M = (N-1)/2;
    for(int step=2;step<=M;step++){
        modify_graph(step);
        hoge(num_measure, measureT, 1, mu_approx, sigma_approx);
    }
    for(int i=0;i<num_measure; i++){
        printf("%d\n", measureT[i]);
    }
}