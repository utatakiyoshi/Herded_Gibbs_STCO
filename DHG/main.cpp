#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <climits>
#include <cassert>
using namespace std;

uint32_t draw_random_int(void) {
    return rand();
}
double draw_random_double(){
    return (double)rand()/RAND_MAX;
}

static const int N=8;

struct variable_settings;

struct variable_settings{
    int algorithm; // 0: Random 1:Herded
    int algorithm_y; // 0: HG 1:RDHG
    int divide_method; // =0: no divide, =2: linear =4: linear+1
    double divide_B;
    double range; // =0: herding, >0: EBR
};
void print_for_debug_vs(variable_settings* vs){
    fprintf(stderr,"algorithm: %d algorithm_y: %d divide_method: %d divide_B: %lf range: %lf \n", vs->algorithm, vs->algorithm_y, vs->divide_method, vs->divide_B, vs->range);
}

struct settings{
    int num_trial;
    int algorithm;
    bool RB;
    int DHG_type;
    double divide_B;
    long long maxT; // configured
    int num_measure;
    long long* measureT;
    int* order;
    double range;

    variable_settings* vs;

    bool Tdiff;
    bool lin;

    bool fullyconnected;
    double Kmin;
    double Kmax;
    double Hmin;
    double Hmax;
    bool random;
    int seed;

    int N;
    int M;
    double* H_list;
    int* Ei_list;
    int* Ej_list;
    double* K_list;

    void print_for_debug(){
        fprintf(stderr, "num_trial: %d algorithm: %d DHG_type: %d divide_B: %lf range: %lf maxT: %lld num_measure: %d\n", num_trial, algorithm, DHG_type, divide_B, range, maxT, num_measure);
        fprintf(stderr, "seed: %d fullyconnected: %c K: [%lf %lf] H: [%lf %lf] random: %c\n", seed, fullyconnected?'T':'F',Kmin,Kmax,Hmin,Hmax,random?'T':'F');
        for(int i=0;i<N;i++){
            print_for_debug_vs(&(vs[i]));
        }
    }

    void configure(){
        maxT = measureT[num_measure-1];
        vs = (struct variable_settings*)malloc(sizeof(struct variable_settings)*N);
        for(int i=0;i<N;i++){
            switch(DHG_type){
                case 0: vs[i].divide_method=0; break;
                case 1: vs[i].divide_method=2; break;
                case 5: vs[i].divide_method=4; break;
            }
            vs[i].divide_B = divide_B;
            switch(algorithm){
                case 1: vs[i].algorithm=0;              break; // All random
                case 2: vs[i].algorithm=(i==(N-1))?1:0; break; // single herding
                case 3: vs[i].algorithm=1;              break; // All herding
                case 4: vs[i].algorithm=1;              break; // Bounded-Error-Gibbs
                case 5: vs[i].algorithm=(i==(N-1))?1:0; break; // single BEG
            }

            switch(DHG_type){
                case 0: vs[i].algorithm_y=0; break;
                case 1: vs[i].algorithm_y=0; break;
                case 5: vs[i].algorithm_y=1; break; // RDHG
            }

            vs[i].range=range;
        }
        print_for_debug();
    }
};


struct cond{
    double w[2];
    double u;
    double u_0;
    long long u_cnt;
    double range;
    int ask_random(double p){
        return (draw_random_double()<p)?1:0;
    }

    int ask(double p){
//      fprintf(stderr, "%.2lf [%.2lf %.2lf] ", p, w[0], w[1]);
        int ret;
        int argmax_w = (w[0]>w[1])?0:1;
        if(w[argmax_w] < range){
            ret = ask_random(p);
        } else {
            ret = argmax_w;
        }
        return ret;
    }
    void update_w_p(double p){
        w[1]+=p;
        w[0]+=1.0-p;
    }
    void update_w_x(int x){
        w[0]-=(x==0)?1:0;
        w[1]-=(x==1)?1:0;
//        fprintf(stderr, "-> [%.2lf %.2lf] ", w[0], w[1]);
    }
    void init(double range_){
        range = range_;
    }
    void reset(double p){
        w[0]=draw_random_double()-p;
        w[1]=-w[0];
        u_0 = draw_random_double();
        u_cnt = 0;
        u=u_0;
//        fprintf(stderr, "%.2lf\n", p);
    }
};

struct variable;
variable* v[N];

struct variable{
    static const int MAXNN = 20;
    static const int MAXD = 1024*1024;
    int spin, nn_size, d, y_card;
    int nn[MAXNN];
    double K[MAXNN];
    double H;
    double range;
    double divide_B;
    int algorithm;
    int algorithm_y;
    vector<cond> c;
    vector<double> p_s;
    vector<int> y_s;
    int z;
    int y;
    int divide_method;
    double p;

    void add(double H_){
        nn_size = 0; H = H_;
    }

    void add_nn(int nn_new, double K_){
        nn[nn_size] = nn_new;
        K[nn_size] = K_;
        nn_size++;
    }

    void init(variable_settings* vs){
//        print_for_debug_vs(vs);
        spin = draw_random_int()%2;

        d=(1<<nn_size);

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

        y_s.resize(d);
        algorithm = vs->algorithm;
        algorithm_y = vs->algorithm_y;
        divide_method = vs->divide_method;
        if(algorithm!=0 &&(vs->divide_method==2 || vs->divide_method==4)){
            divide_B = vs->divide_B;
            switch(vs->divide_method){
                case 2: y_card=1.0/divide_B; break;
                case 4: y_card=(1.0/divide_B)+1; break;
            }
            for(int i=0;i<d;i++){
                y_s[i] = int(p_s[i]/divide_B);
                if(y_s[i]==int(1.0/divide_B)){y_s[i]=y_card-1;}
            }
        } else {
            for(int i=0;i<d;i++){
                y_s[i]=i;
            }
            y_card=d;
        }
        c.resize(y_card);
        for(int i=0;i<y_card;i++){ c[i].init(vs->range); }
    }

    void reset(){
        spin = draw_random_int()%2;
        for(int i=0;i<y_card;i++){
            if(divide_method==2){
              c[i].reset(double(i*2+1)*divide_B/2.0);
            } else if(divide_method==4){
              c[i].reset(double(i)*divide_B);
            } else {
              c[i].reset(p_s[i]);
            }
        }
    }

    int calc_z(){
        int z_t=0;
        for(int i=0;i<nn_size;i++){
            z_t |= ((v[nn[i]]->spin) << i);
        }
        return z_t;
    }

    int calc_z_from_packed(int x){
        int z_t=0;
        for(int i=0;i<nn_size;i++){
            z_t |= (((x>>(nn[i]))&1) << i);
        }
        return z_t;
    }

    int calc_y_from_packed(int x){
        return y_s[calc_z_from_packed(x)];
    }

    void update_z(){
        z=calc_z();
    }
    void draw_y(){
      y=y_s[z];
      p=p_s[z];
      if(algorithm_y==1){
          double lt = (double)(y)*divide_B;
          double ut = (double)(y+1)*divide_B;
          double q = (ut-p)/(ut-lt);
          //fprintf(stderr, "%d %.2lf %.2lf %.2lf %.2lf %.2lf->", y, p, lt, ut, q, q*lt+(1.0-q)*ut);
          y = (draw_random_double()<q)?(y):(y+1);
          p=double(y)*divide_B;
          //fprintf(stderr, "%d %.2lf\n", y, p);
      }
    }
    void next_state(){ spin = c[y].ask(p); }
    void next_state_random(){ spin = (draw_random_double()<p)?1:0; }
    void update_w_p(){ c[y].update_w_p(p); }
    void update_w_x(){ c[y].update_w_x(spin); }

    void update(){
        update_z();
        draw_y();
//        fprintf(stderr, "(%.2lf -> %.2lf) ", p_s[z], p);
        switch(algorithm){
            case 0: next_state_random(); break;
            case 1: next_state(); update_w_p(); update_w_x(); break;
        }
    }


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
};

void add_nn_sub(settings* s, int i, int Ei, int Ej, double K){
    s->Ei_list[i]=Ei;
    s->Ej_list[i]=Ej;
    s->K_list[i]=K;
}

void draw_graph(settings* s){
    srand(s->seed);
    draw_random_double();
    draw_random_double();
    draw_random_double();
    double H=0;
    double K=0;
    int head=0;
    double x;
    s->N=N;
    s->H_list = (double*)malloc(sizeof(double)*N);
    for(int i=0;i<N;i++){
        if(s->random){
            x = draw_random_double()*((s->Hmax)-(s->Hmin))+(s->Hmin);
        } else {
            x = s->Hmax;
        }
        s->H_list[i]=x;
    }

    if(s->fullyconnected){
        s->M=(s->N)*((s->N)-1)/2;
    } else {
        s->M=s->N;
    }
    s->Ei_list = (int*)malloc(sizeof(int)*(s->M));
    s->Ej_list = (int*)malloc(sizeof(int)*(s->M));
    s->K_list = (double*)malloc(sizeof(double)*(s->M));
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if((!(s->fullyconnected) && (i+1)%N ==j)||(s->fullyconnected && (i<j) )){
                if(s->random){
                    x = draw_random_double()*((s->Kmax)-(s->Kmin))+(s->Kmin);
                } else {
                    x = s->Kmax;
                }
                add_nn_sub(s, head, i, j, x);
                head++;
            }
        }
    }
}

void save_graph(settings* s, FILE* fout){
    fprintf(fout, "%d ", s->N);
    fprintf(fout, "%d\n", s->M);
    for(int i=0;i<(s->N);i++){
        fprintf(fout, "%lf\n", s->H_list[i]);
    }
    for(int i=0;i<(s->M);i++){
        fprintf(fout, "%d %d %lf\n", s->Ei_list[i], s->Ej_list[i], s->K_list[i]);
    }
}

void load_graph(settings* s, FILE* fin){
    fscanf(fin, "%d", &(s->N));
    fscanf(fin, "%d", &(s->M));
    s->H_list = (double*)malloc(sizeof(double)*(s->N));
    s->Ei_list = (int*)malloc(sizeof(int)*(s->M));
    s->Ej_list = (int*)malloc(sizeof(int)*(s->M));
    s->K_list = (double*)malloc(sizeof(double)*(s->M));
    for(int i=0;i<(s->N);i++){
        fscanf(fin, "%lf", &(s->H_list[i]));
    }
    for(int i=0;i<(s->M);i++){
        fscanf(fin, "%d %d %lf", &(s->Ei_list[i]), &(s->Ej_list[i]), &(s->K_list[i]));
    }
}

void build_graph(settings* s){
    for(int i=0;i<N;i++){
        v[i]=(variable*)malloc(sizeof(variable));
//        fprintf(stderr,"H_%d = %lf\n", i, s->H_list[i]);
        v[i]->add(s->H_list[i]);
    }
    for(int i=0;i<s->M;i++){
//        fprintf(stderr,"K_%d,%d = %lf\n", s->Ei_list[i], s->Ei_list[j], s->K_list[i]);
        v[s->Ei_list[i]]->add_nn(s->Ej_list[i], s->K_list[i]);
        v[s->Ej_list[i]]->add_nn(s->Ei_list[i], s->K_list[i]);
    }
}

void init_variables(settings* s){
    for(int i=0;i<N;i++){
        v[i]->init(&(s->vs[i]));
    }
}

void reset_variable_states(){
    for(int i=0;i<N;i++){
        v[i]->reset();
    }
}

struct dist{
    int x;
    int* zw;

    int y_card;
    int x_size;
    int z_size;
    int zw_size;
    int zwx_size;
    int x_card;
    int z_card;
    int zw_card;
    int zwx_card;
    int card;

    double* p;
    double* py;
    double* pz;
    double** px_on_zw;
    double** py_on_zw;
    double** pzwx_on_y;
    double** pzw_on_y;
    double** px_on_y;
    double** pzwpx_on_y;
    double* px;
    double* pzw;
    double* pzwx;

    int calc_full_state(){
        int y_t = v[x]->y;
        int x_t = v[x]->spin;
        int state=0;
        for(int i=0;i<zw_size;i++){
            state |= ((v[zw[i]]->spin)<<i);
        }
        state = state * x_card + x_t;
        state = state * y_card + y_t;
        return state;
    }

    void free_all(){
        free(pzwx);
        free(pzw);
        free(pz);
        free(px);
        for(int i=0;i<y_card;i++){
            free(pzwpx_on_y[i]);
            free(pzwx_on_y[i]);
            free(pzw_on_y[i]);
            free(px_on_y[i]);
        }
        for(int i=0;i<zw_size;i++){
          free(px_on_zw[i]);
          free(py_on_zw[i]);
        }
        free(px_on_zw);
        free(py_on_zw);
        free(py);
        free(pzwpx_on_y);
        free(pzwx_on_y);
        free(pzw_on_y);
        free(px_on_y);
        free(p);
        free(zw);
    }

    void allocate(int y_card_, int x_card_, int z_size_, int zw_size_){
        y_card = y_card_;
        x_size = 1;
        x_card = x_card_;
        z_size = z_size_;
        z_card = (1<<z_size);
        zw_size = zw_size_;
        zw_card = (1<<zw_size);
        zwx_size = zw_size+x_size;
        zwx_card = zw_card*x_card;
        card = y_card * zwx_card;

        zw=(int*)malloc(sizeof(int)*zw_size);

        p=(double*)malloc(sizeof(double)*card);

        px_on_y   =(double**)malloc(sizeof(double*)*y_card);
        pzw_on_y  =(double**)malloc(sizeof(double*)*y_card);
        pzwx_on_y =(double**)malloc(sizeof(double*)*y_card);
        pzwpx_on_y =(double**)malloc(sizeof(double*)*y_card);

        px_on_zw =(double**)malloc(sizeof(double*)*zw_card);
        py_on_zw =(double**)malloc(sizeof(double*)*zw_card);

        py =(double*)malloc(sizeof(double)*y_card);
        for(int i=0;i<y_card;i++){
            px_on_y[i]   = (double*)malloc(sizeof(double)*x_card);
            pzw_on_y[i]  = (double*)malloc(sizeof(double)*zw_card);
            pzwx_on_y[i] = (double*)malloc(sizeof(double)*zwx_card);
            pzwpx_on_y[i] = (double*)malloc(sizeof(double)*zwx_card);
        }
        for(int i=0;i<zw_card;i++){
          px_on_zw[i] = (double*)malloc(sizeof(double)*x_card);
          py_on_zw[i] = (double*)malloc(sizeof(double)*y_card);
        }
        px   =(double*)malloc(sizeof(double)*x_card);
        pz   =(double*)malloc(sizeof(double)*z_card);
        pzw  =(double*)malloc(sizeof(double)*zw_card);
        pzwx =(double*)malloc(sizeof(double)*zwx_card);
    }
    void print_for_debug(){
        fprintf(stderr,"========\n");
        for(int y_t=0;y_t<y_card;y_t++){
            fprintf(stderr, "[%lf]\t", py[y_t]);
            for(int x_t=0;x_t<x_card;x_t++){
                fprintf(stderr, "%lf\t", px_on_y[y_t][x_t]);
            }
            fprintf(stderr,"\n");
            for(int zw_t=0;zw_t<zw_card;zw_t++){
                fprintf(stderr, "%lf\t", pzw_on_y[y_t][zw_t]);
                for(int x_t=0;x_t<x_card;x_t++){
                    int zwx_t = zw_t*x_card+x_t;
                }
            }
        }
    }
    void calc_conditionals(){
        memset(py, 0, sizeof(double)*y_card);
        for(int i=0;i<y_card;i++){
            memset(px_on_y[i], 0, sizeof(double)*x_card);
            memset(pzw_on_y[i], 0, sizeof(double)*zw_card);
            memset(pzwx_on_y[i], 0, sizeof(double)*zwx_card);
        }
        for(int i=0;i<zw_card;i++){
          memset(px_on_zw[i], 0, sizeof(double)*x_card);
          memset(py_on_zw[i], 0, sizeof(double)*y_card);
        }
        memset(pz, 0, sizeof(double)*z_card);
        memset(px, 0, sizeof(double)*x_card);
        memset(pzw, 0, sizeof(double)*zw_card);
        memset(pzwx, 0, sizeof(double)*zwx_card);

        for(int i=0;i<card;i++){
            int y_t = i % y_card;
            int zwx_t = i / y_card;
            int x_t = zwx_t % x_card;
            int zw_t =zwx_t / x_card;
            int z_t = zw_t % z_card;
            px[x_t] += p[i];
            pz[z_t] += p[i];
            pzw[zw_t] += p[i];
            pzwx[zwx_t] += p[i];
            py[y_t] += p[i];
            px_on_zw[zw_t][x_t] += p[i];
            py_on_zw[zw_t][y_t] += p[i];
            px_on_y[y_t][x_t] += p[i];
            pzw_on_y[y_t][zw_t] += p[i];
            pzwx_on_y[y_t][zwx_t] += p[i];
        }

        for(int zw_t=0;zw_t<zw_card;zw_t++){
          if(pzw[zw_t]>0){
            for(int j=0;j<x_card  ;j++){ px_on_zw[zw_t][j] /= pzw[zw_t]; }
            for(int j=0;j<y_card  ;j++){ py_on_zw[zw_t][j] /= pzw[zw_t]; }
          }
        }
        for(int y_t=0;y_t<y_card;y_t++){
            if(py[y_t]>0){
                for(int j=0;j<x_card  ;j++){ px_on_y[y_t][j] /= py[y_t]; }
                for(int j=0;j<zw_card ;j++){ pzw_on_y[y_t][j] /= py[y_t]; }
                for(int j=0;j<zwx_card;j++){ pzwx_on_y[y_t][j] /= py[y_t]; }
            }
            for(int j=0;j<zwx_card;j++){
                int x_t = j % x_card;
                int zw_t = j / x_card;
                pzwpx_on_y[y_t][j] = px_on_y[y_t][x_t] * pzw_on_y[y_t][zw_t];
            }
        }
    }

    void calc_qx_on_y(dist* pi, double*** qx_on_y_out, bool flag_RDHG=false, double divide_B=0){
      double** qx_on_y = (double**)malloc(sizeof(double*)*y_card);
      for(int i=0;i<y_card;i++){
        qx_on_y[i]   = (double*)malloc(sizeof(double)*x_card);
      }
      for(int i=0;i<y_card;i++){
        memset(qx_on_y[i], 0, sizeof(double)*x_card);
      }
      if(flag_RDHG){
        for(int y_t=0;y_t<y_card;y_t++){
            qx_on_y[y_t][1] = y_t*divide_B;
            qx_on_y[y_t][0] = 1.0-y_t*divide_B;
        }
      } else {
        for(int i=0;i<card;i++){
            int y_t = i % y_card;
            int zwx_t = i / y_card;
            int x_t = zwx_t % x_card;
              int zw_t =zwx_t / x_card;
            qx_on_y[y_t][x_t] += (pi->px_on_zw[zw_t][x_t])*pzw_on_y[y_t][zw_t];
        }
      }
      *qx_on_y_out = qx_on_y;
    }
    void free_qx_on_y(double*** qx_on_y_out){
      double** qx_on_y = *qx_on_y_out;
      for(int i=0;i<y_card;i++){
        free(qx_on_y[i]);
      }
      free(qx_on_y);
    }

    void reset_p(){
        memset(p, 0, sizeof(double)*card);
    }


    void configure_and_allocate(int x_){
        x=x_;
        allocate(v[x]->y_card, 2, v[x]->nn_size, N-1);

        int head=0;
        bool flag[N];
        for(int i=0;i<N;i++){
            flag[i]=true;
        }
        flag[x] = false;
        for(int i=0;i<(v[x]->nn_size);i++){
            zw[head++] = v[x]->nn[i];
            flag[v[x]->nn[i]] = false;
        }
        for(int i=0;i<N;i++){
            if(flag[i]){ zw[head++] = i;}
        }
    }

    void init_from_dist_p(double** q){
        reset_p();
        for(int i=0;i<(1<<N);i++){
            int x_t=(i>>x)&1;
            int zw_t=0;
            for(int j=0;j<zw_size;j++){
                zw_t |= ((i>>(zw[j]))&1)<<j;
            }
            int state_zwx = zw_t * x_card+x_t;
            int z_t = v[x]->calc_z_from_packed(i);
            for(int y_t=0;y_t<y_card;y_t++){
              p[state_zwx*y_card+y_t] += q[i][y_t];
            }
        }
        calc_conditionals();
    }

    void init_from_dist_q(double* q, bool flag_RDHG=false){
        reset_p();
        for(int i=0;i<(1<<N);i++){
            int x_t=(i>>x)&1;
            if(x_t==0){
              int i0=i;
              int i1=i|(1<<x);
              int zw_t=0;
              for(int j=0;j<zw_size;j++){
                  zw_t |= ((i>>(zw[j]))&1)<<j;
              }
              int state_zwx0 = zw_t * x_card+0;
              int state_zwx1 = zw_t * x_card+1;
              int y_t = v[x]->calc_y_from_packed(i);
              int z_t = v[x]->calc_z_from_packed(i);
              if(flag_RDHG){
                double divide_B=v[x]->divide_B;
                  double lt = double(y_t)*divide_B;
                  double ut = double(y_t+1)*divide_B;
                  double q_t = (ut-v[x]->p_s[z_t])/(ut-lt);
                  p[state_zwx1 * y_card+y_t] += (q[i0]+q[i1])*(lt)*q_t;
                  p[state_zwx0 * y_card+y_t] += (q[i0]+q[i1])*(1.0-lt)*q_t;
                  p[state_zwx1 * y_card+(y_t+1)] += (q[i0]+q[i1])*(ut)*(1.0-q_t);
                  p[state_zwx0 * y_card+(y_t+1)] += (q[i0]+q[i1])*(1.0-ut)*(1.0-q_t);
              } else {
                  p[state_zwx0 * y_card+y_t] += q[i0];
                  p[state_zwx1 * y_card+y_t] += q[i1];
              }
            }
          }
        calc_conditionals();
    }

};

void calc_exact(double* q_out){
    int M = (1<<N);
    double Z=0;
    double mu=0;
    double sigma=0;
    for(int i=0;i<M;i++){
        double E=0;
        for(int j=0;j<N;j++){
            v[j]->spin = ((i>>j)&1);
        }
        for(int j=0;j<N;j++){
            E+=v[j]->get_E();
        }
        double p = exp(-E);
        q_out[i]=p;
        Z+=p;
    }
    for(int i=0;i<M;i++){
        q_out[i]/=Z;
    }
}

double calc_D(int N_, double* p, double* q){
    double D=0;
    for(int i=0;i<(1<<N_);i++){
        D+=abs(p[i]-q[i]);
    }
    return D;
}

double calc_D_all(dist* p, dist* q){
    return calc_D(p->zwx_size, p->pzwx, q->pzwx);
}

double calc_D_prev(dist* p, dist* q){
    return calc_D(p->zw_size, p->pzw, q->pzw);
}

double calc_D_herding(dist* p, dist* q, double** qx_on_y){
    double D=0;
    for(int y_t=0;y_t<(p->y_card);y_t++){
        D += (p->py[y_t]) * calc_D(p->x_size, p->px_on_y[y_t], qx_on_y[y_t]);
    }
    return D;
}

double calc_D_bin(dist* p, dist* q, int y_card){
    double D=0;
    for(int zw_t=0;zw_t<(p->zw_card);zw_t++){
        for(int y=0;y<y_card;y++){
          D += (p->pzw[zw_t]) * abs(p->py_on_zw[zw_t][y]-q->py_on_zw[zw_t][y]);
        }
    }
    return D;
}

double calc_D_cor(dist* p){
    double D=0;
    for(int y_t=0;y_t<(p->y_card);y_t++){
        D += (p->py[y_t]) * calc_D(p->zwx_size, p->pzwx_on_y[y_t], p->pzwpx_on_y[y_t]);
    }
    return D;
}

double calc_D_approx(dist* p, dist* q, double** qx_on_y){
    double D=0;
    for(int zw_t=0;zw_t<(p->zw_card);zw_t++){
      for(int x_t=0;x_t<(p->x_card);x_t++){
        double tmp=0;
        for(int y_t=0;y_t<(p->y_card);y_t++){
          tmp += qx_on_y[y_t][x_t]*(q->py_on_zw[zw_t][y_t]);
        }
        D += (p->pzw[zw_t]) * abs(tmp - q->px_on_zw[zw_t][x_t]);
      }
    }
    return D;
}

int calc_full_state(){
    int state=0;
    for(int i=0;i<N;i++){
        state |= ((v[i]->spin) << i);
    }
    return state;
}

void run(settings* s, int* states_out, int* y_out){
    reset_variable_states();
    fprintf(stderr,"================\n");
    for(int t=0; t < (s->maxT)*N; t++){
        int i=s -> order[t%N];
        v[i]->update();
        int state = calc_full_state();
        states_out[t]=state;
        y_out[t]=v[i]->y;
        if(t%N==N-1 && t/N<10){
          for(int j=0;j<N;j++){
            fprintf(stderr, "%d ", v[j]->spin);
          }
          fprintf(stderr,"\n");
        }
    }
}

void lin_main(settings* s, int x, double* q_approx){
    init_variables(s);

    dist q;
    q.configure_and_allocate(x);
    q.init_from_dist_q(q_approx, false);
    double ans=q.px[1];

    reset_variable_states();
    double* sum_abs_diff=(double*)malloc(sizeof(double)*(s->num_measure));
    memset(sum_abs_diff,0,sizeof(double)*s->num_measure);
    for(int trial=0;trial<s->num_trial;trial++){
      fprintf(stderr, "trial: %d\n", trial);
//?    reset_variable_states();
      double sum=0;
      int cnt=0;
      int head=0;
      for(int t=0; t < (s->maxT)*N; t++){
          int i=s -> order[t%N];
          v[i]->update();
          if(i==x){
              int x_t = v[i]->spin;
              if(s->RB){
                int z = v[i]->z;
                double p = v[i]->p_s[z];
                sum += p;
              } else {
                sum += x_t;
              }
              cnt++;
              if(cnt==s->measureT[head]){
                  sum_abs_diff[head] += abs(sum/cnt-ans);
                  head++;
              }
          }
      }
    }
    for(int i=0;i<s->num_measure;i++){
      printf("%.15lf\n", sum_abs_diff[i]/s->num_trial);
    }
}

void tdiff_main(settings* s, int x){
    init_variables(s);

    reset_variable_states();
    double** cnt = (double**)malloc(sizeof(double*)*(v[x]->y_card));
    double Z = 0;
    double*** Tdist = (double***)malloc(sizeof(double**)*(v[x]->y_card));

    for(int y_t=0;y_t<(v[x]->y_card);y_t++){
        cnt[y_t] = (double*)malloc(sizeof(double)*2);
        memset(cnt[y_t],0,sizeof(double)*2);
        Tdist[y_t] = (double**)malloc(sizeof(double*)*2);
        for(int x_t=0;x_t<2;x_t++){
            Tdist[y_t][x_t] = (double*)malloc(sizeof(double)*(1<<N));
            memset(Tdist[y_t][x_t], 0, sizeof(double)*(1<<N));
        }
    }
    
    int* prev_x = (int*)malloc(sizeof(int)*(v[x]->y_card));
    for(int t=0; t < (s->maxT)*N; t++){
        int i=s -> order[t%N];
        v[i]->update();
        if(i==x){
            int state = calc_full_state();
            state = state & (~(1<<x));
            int y_t = v[i]->y;
            int x_t = v[i]->spin;
            if(prev_x[y_t]>=0){
                cnt[y_t][prev_x[y_t]] += 1;
                Tdist[y_t][prev_x[y_t]][state] += 1;
                Z += 1;
            }
            prev_x[y_t]=x_t;
        }
    }
    
    for(int y_t=0;y_t<(v[x]->y_card);y_t++){
        for(int x_t=0;x_t<2;x_t++){
            for(int i=0;i<(1<<N);i++){
                if(cnt[y_t][x_t]>0){
                    Tdist[y_t][x_t][i] /= cnt[y_t][x_t];
                }
            }
        }
    }

    for(int y_t=0;y_t<(v[x]->y_card);y_t++){
        fprintf(stdout, "%d ", y_t);
        double tt = min(1.0-(v[x]->p_s[y_t]), v[x]->p_s[y_t]);
        double pyt = (cnt[y_t][0]+cnt[y_t][1])/Z;
        double sum_diff = 0;
        double sum_var = 0;
        for(int i=0;i<(1<<N);i++){
            double ty0=Tdist[y_t][0][i];
            double ty1=Tdist[y_t][1][i];
            sum_diff += abs(ty0-ty1);
        }
        fprintf(stdout, "%.15f %.15lf %.15lf ", pyt, tt, sum_diff);
        fprintf(stdout, "\n");
    }
}


void D_main(settings* s, int x, double* q_approx, double** stats_out){
    double stats[6][101][50];

    init_variables(s);

    dist q;
    q.configure_and_allocate(x);
    bool flag_RDHG = false;
    if(s->DHG_type==5){flag_RDHG=true;}
    q.init_from_dist_q(q_approx, flag_RDHG);

    int* states=(int*)malloc(sizeof(int)*((s->maxT)*N));
    int* y_s=(int*)malloc(sizeof(int)*((s->maxT)*N));
    double** p_cnts=(double**)malloc(sizeof(double*)*((1<<N)));
    double** p_tmp=(double**)malloc(sizeof(double*)*((1<<N)));
    int y_card = v[x]->y_card;
    for(int i=0;i<(1<<N);i++){
      p_cnts[i] = (double*)malloc(sizeof(double)*y_card);
      p_tmp[i] = (double*)malloc(sizeof(double)*y_card);
    }

    dist p;
    p.configure_and_allocate(x);
    for(int trial=0;trial<s->num_trial;trial++){
        run(s, states, y_s);

        int head=0;
        for(int j=0;j<(1<<N);j++){memset(p_cnts[j],0,sizeof(double)*y_card);}

        for(int t=0; t < (s->maxT)*N; t++){
            // p_cnt++;
            int state = states[t];
            int y_t = y_s[t];
            if(t%N==x){
              p_cnts[state][y_t]+=1;
              if((t/N)+1==s->measureT[head]){
                  long long Tsample = s->measureT[head];
                    for(int i=0;i<(1<<N);i++){
                      for(int j=0;j<y_card;j++){
                        p_tmp[i][j] = p_cnts[i][j]/Tsample;
                      }
                     }
                  p.init_from_dist_p(p_tmp);
                  double** qx_on_y;
                  p.calc_qx_on_y(&q, &qx_on_y, flag_RDHG, v[x]->divide_B);
                  /* D */ stats[0][trial][head] = calc_D_all(&(p), &(q));
                  /* D_prev */ stats[1][trial][head] = calc_D_prev(&(p), &(q));
                  /* D_herding */  stats[2][trial][head]= calc_D_herding(&(p), &(q), (qx_on_y));
                  /* D_approx */ stats[3][trial][head] = calc_D_approx(&(p),&(q), (qx_on_y));
                  /* D_cor */ stats[4][trial][head] = calc_D_cor(&(p));
                  /* D_bin */ stats[5][trial][head] = calc_D_bin(&(p),&(q),v[x]->y_card);
                  p.free_qx_on_y(&qx_on_y);
                  if((t+1)%N==0){
                      head++;
                  }
              }
          }
        }
//        fprintf(stderr,"DONE.", trial+1);
    }
    p.free_all();
    q.free_all();

    for(int c=0;c<6;c++){
        for(int i=0;i<s->num_measure; i++){
            double sum=0;
            for(int j=0;j<s->num_trial;j++){
                sum+=stats[c][j][i];
            }
            stats_out[c][i]=sum/(s->num_trial);
        }
    }

    free(states);
    for(int i=0;i<N+1;i++){
        free(p_cnts[i]);
    }
    free(p_cnts);
}

int main(int argc, char** argv){
    bool flag_draw_and_save=false;
    settings s;

    s.algorithm = 1;
    s.RB = false;
    s.DHG_type = 0;
    s.divide_B = 1.0/8;
    s.range = 0;
    s.num_trial = 10;
    s.num_measure = 24;
    s.fullyconnected=true;
    s.Kmin=0;
    s.Kmax=1;
    s.Hmin=0;
    s.Hmax=1;
    s.random = false;

    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-draw")==0){flag_draw_and_save=true;}
        if(strcmp(argv[i], "-num_trial")==0){i++; s.num_trial=atoi(argv[i]);}
        if(strcmp(argv[i], "-num_measure")==0){i++; s.num_measure=atoi(argv[i]);}
        if(strcmp(argv[i], "-algorithm")==0){i++; s.algorithm=atoi(argv[i]);}
        if(strcmp(argv[i], "-RB")==0){s.RB=true;}
        if(strcmp(argv[i], "-DHG_type")==0){i++; s.DHG_type=atoi(argv[i]);}
        if(strcmp(argv[i], "-divide_B")==0){i++; s.divide_B=atof(argv[i]);}
        if(strcmp(argv[i], "-divide_B_inv")==0){i++; double j = double(atoi(argv[i])); s.divide_B=1.0/j;}
        if(strcmp(argv[i], "-range")==0){i++; s.range=atof(argv[i]);}
        if(strcmp(argv[i], "-fullyconnected")==0){s.fullyconnected = true;}
        if(strcmp(argv[i], "-ring")==0){s.fullyconnected = false;}
        if(strcmp(argv[i], "-random")==0){s.random = true;}
        if(strcmp(argv[i], "-fixed")==0){s.random = false;}
        if(strcmp(argv[i], "-Kmin")==0){i++; s.Kmin = atof(argv[i]);}
        if(strcmp(argv[i], "-Kmax")==0){i++; s.Kmax = atof(argv[i]);}
        if(strcmp(argv[i], "-Hmin")==0){i++; s.Hmin = atof(argv[i]);}
        if(strcmp(argv[i], "-Hmax")==0){i++; s.Hmax = atof(argv[i]);}
        if(strcmp(argv[i], "-seed")==0){i++; s.seed = atoi(argv[i]);}
        if(strcmp(argv[i], "-tdiff")==0){s.Tdiff=true;}
        if(strcmp(argv[i], "-lin")==0){s.lin=true;}
    }

    if(flag_draw_and_save){
        s.print_for_debug();
        draw_graph(&s);
        save_graph(&s, stdout);
        return 0;
    }
    
    srand((unsigned)time(NULL));
    
    long long measureT[50];
    measureT[0] = 1;
    for(int head = 1; head < s.num_measure; head++){ measureT[head] = measureT[head-1] * 2; }
    s.measureT = measureT;
    
    int order[N];
    for(int i=0;i<N;i++){order[i]=i;}
    s.order = order;

    load_graph(&s, stdin);

    s.configure();

    build_graph(&s);

    double* q_approx=(double*)malloc(sizeof(double)*(1<<N));
    calc_exact(q_approx);

    if(s.Tdiff){
        tdiff_main(&s, N-1);
    }else if (s.lin){
        lin_main(&s, N-1, q_approx);
    } else {
        double* stats[6];
        for(int i=0;i<6;i++){
            stats[i] = (double*)malloc(sizeof(double)*(s.num_measure));
        }
        D_main(&s, N-1, q_approx, stats);
        for(int j=0;j<(s.num_measure);j++){
            fprintf(stdout, "2^%d\t", j+1);
            for(int c=0;c<6;c++){
                fprintf(stdout, "%.15lf\t", stats[c][j]);
            }
            fprintf(stdout,"\n");
        }
    }
}
