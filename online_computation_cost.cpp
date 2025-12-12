#include "Tree.h"
#include "PRGTest.h"
#include <iostream>
#include <random>
#include <limits>
#include <bitset>
#include <cassert>
#include <utility>
#include <thread>
#include <iomanip>
#include <bitset>
#include <algorithm>
using namespace std;
ull seed = 1UL;
const ull STV_SIZE = 258;
const ull CW_SIZE = 32508;
// bit length for elements
const ull L = 64;
const ull lambda = 128;
// average latency(ns)
const ull trans_delay = 75422000;
// 0.168820736 bit/ns
const double bandwidth = 161.0 * 1024 * 1024 /1000 /1000/1000;
ull get_random()
{
    // 随机数生成引擎
    // std::random_device rd;     // 随机数种子
    // std::mt19937_64 gen(rd()); // 64位 Mersenne Twister 引擎

    std::mt19937_64 gen(seed++);
    // 定义无符号长整型的范围
    std::uniform_int_distribution<unsigned long long> dis(0, std::numeric_limits<unsigned long long>::max());

    // 生成随机 unsigned long long 数
    unsigned long long random_number = dis(gen);

    // 输出生成的随机数
    // std::cout << "随机生成的 unsigned long long 数: " << random_number << std::endl;
    return random_number;
}

void genMatrixTriples(vvu& A, vu& b, vu& c, vvu& M, vu& x)
{
    for(int i = 0; i < A.size(); i++)
    {
        for(int j = 0; j < A[i].size(); j++)
        {
            A[i][j] = get_random();
            M[i][j] = get_random();
            if(i == 0)
            {
                b[j] = get_random();
                x[j] = get_random();
            }
        }
        c[i] = get_random();
    }
}

vu matrixProduct(vvu matrix, vu vect)
{
    vu result;
    ull temp = 0;
    for(int i = 0; i < matrix.size(); i++)
    {
        temp = 0;
        for(int j = 0; j < matrix[i].size(); j++)
        {
            temp += (matrix[i][j] * vect[j]);
        }
        result.push_back(temp);
    }
    return result;
}

void getOursShare(vvu A, vu b, vu c, vvu Delta_M, vu Delta_x)
{
    vu mx = matrixProduct(Delta_M, Delta_x);
    vu mb = matrixProduct(Delta_M, b);
    vu ax = matrixProduct(A, Delta_x);
    vu result;
    for(int i = 0; i < mx.size(); i++)
    {
        result.push_back(mx[i]-mb[i]-ax[i]+c[i]);
    }
}

// need open m*n+n elements
void getZ22Share(vvu A, vu b, vu c, vvu Delta_M, vu Delta_x)
{
    for(int i = 0; i < Delta_M.size(); i++)
    {
        for(int j = 0; j < Delta_M[i].size(); j++)
        {
            Delta_M[i][j] -= A[i][j];
        }
    }
    for(int i = 0; i < Delta_x.size(); i++)
    {
        Delta_x[i] -= b[i];
    }
    vu mx = matrixProduct(Delta_M, Delta_x);
    vu mb = matrixProduct(Delta_M, b);
    vu ax = matrixProduct(A, Delta_x);
    vu result;
    for(int i = 0; i < mx.size(); i++)
    {
        result.push_back(mx[i]-mb[i]-ax[i]+c[i]);
    }
}

void getZ23Share(vu x, ull Ij0)
{
    ull rm = get_random();
    ull sm = get_random();
    vu p(x.size());
    for(ull i = 0; i < x.size(); i++)
    {
        ull ii = (i + sm) % x.size();
        p[i] = x[ii] + rm;
    }
    ull r = get_random();
    ull Ijp0 = Ij0 + r;
    ull s3 = Ij0 + r + sm;
    // send Ijp0 to other and receive Ijp0+Ijp1+s1+r
    ull Ip1 = (get_random() - r) % x.size();

    ull p1pip1 = r;
    ull p0xip0 = sm;
    ull share = p0xip0 + p1pip1 + -rm; 

}

void testFeatureSelection(vector<pair<int, int> > params)
{
    vector<pair<int, int> > z23ot_cost;
    z23ot_cost.push_back(make_pair(2,7));
    z23ot_cost.push_back(make_pair(3,7));
    z23ot_cost.push_back(make_pair(4,8));
    z23ot_cost.push_back(make_pair(8,9));
    z23ot_cost.push_back(make_pair(13,13));
    z23ot_cost.push_back(make_pair(17,80));
    for(int i = 0; i < params.size(); i++)
    {
        int m = pow(2, params[i].first-1) - 1;
        int n = params[i].second;
        vvu A(m, vu(n, 0)), Delta_M(m, vu(n, 0));
        vu b(n), c(m), Delta_x(n);
        ull tx = get_random();
        genMatrixTriples(A, b, c, Delta_M, Delta_x);
        clockrep our = measureTime(getOursShare, A, b, c, Delta_M, Delta_x);
        clockrep z22 = measureTime(getZ22Share, A, b, c, Delta_M, Delta_x);
        clockrep z23 = measureTime(getZ23Share,Delta_x, tx);
        // getShare(A, b, c, Delta_M, Delta_x);
        
        cout << "fs:z22," << params[i].first <<"," << n << "," <<fixed<<setprecision(0)<< z22 + (m*n+n)*L/bandwidth + trans_delay << endl;
        cout << "fs:z23," << params[i].first <<"," << n << "," << z23 + (z23ot_cost[i].second*1000000) + (6*trans_delay) + (2*L/bandwidth) << endl;
        cout << "fs:ours," << params[i].first <<"," << n << "," << our << endl;


    }
}

bool getBitProduct(bool x, bool y){
    bool a1 = true;
    bool a2 = true;
    bool b1 = false;
    bool b2 = true;
    bool c1 = false;
    bool c2 = false;
    bool e1 = x ^ a1;
    bool f1 = y ^ b1;
    //get e2 f2
    bool e2 = false;
    bool f2 = true;

    bool e = e1^e2;
    bool f = f1^f2;

    bool ef = e & f;

    return ef ^ (b1&e) ^(a1&f) ^ c1;

}
ull getRingPrroduct(ull x, ull y)
{
    ull a1=14490808261858112199,a2=2177743953316042629;
    ull b1=12415856028556828342,b2=16338301252824554741;
    ull c1=14315882575126838720,c2=18193727490463377156;
    ull e1 = x * a1;
    ull f1 = y * b1;
    //get e2 f2
    ull e2 = 23587943239643466562;
    ull f2 = 18553232343782334634;

    ull e = e1+e2;
    ull f = f1+f2;

    bool ef = e * f;

    return ef - (b1*e) -(a1*f) + c1;
}

void z22DecisionNodeEvl(ull yji, ull xsigmi)
{
    ull ai = yji - xsigmi;
    // p aka p_0
    bitset<L> p(ai);
    bitset<L> q(0);
    bitset<L> w(ai);
    vector<bool> c,d,e;
    bool temp = getBitProduct(p[0], q[0]);
    c.push_back(temp);
    c.push_back(temp);
    for(int k = 2; k < L; k++)
    {
        d.push_back(getBitProduct(p[k], q[k]) ^ 1);
        e.push_back(getBitProduct(w[k], c[k-1]) ^ 1);
        c.push_back(getBitProduct(e[k], d[k]) ^ 1);
    }
    
    bool a63 = w[63] + c[63];
    ull t10 = a63, t20 = 0;
    ull t11 = 0;
    ull al0 = t10 + t20 - 2 * getBitProduct(t10, t20);
}

pair<bool, bool> operat(bool g2, bool p2, bool g1, bool p1)
{
    bool g1p2 = getBitProduct(g1, p2);
    bool g = g2 ^ g1p2;
    bool p = getBitProduct(p2, p1);
    return make_pair(g, p);
}

void z23DecisionNodeEvl(ull yji, ull xsigmi)
{
    ull ai = yji - xsigmi;
    // p aka p_0
    bitset<L> a(ai);
    bitset<L> b(0);
    bitset<L> w(ai);
    vector<bool> g, p;

    for(int i = 0; i < L; i++)
    {
        g.push_back(getBitProduct(a[i], b[i]));
        p.push_back(a[i] ^ b[i]);
    }

    vector<bool> g1, p1;
    g1.push_back(g[0]);
    p1.push_back(p[0]);

    pair<bool, bool> pp;
    for(int k = 1; k < (L >> 1); k++)
    {
        pp = operat(g[2*k], p[2*k], g[2*k-1], p[2*k-1]);
        g1.push_back(pp.first);
        p1.push_back(pp.second);
    }

    vector<bool> g2, p2;
    for(int k = 0; k < (L >> 2); k++)
    {
        pp = operat(g1[2*k+1], p1[2*k+1], g1[2*k], p1[2*k]);
        g2.push_back(pp.first);
        p2.push_back(pp.second);
    }

    vector<bool> g3, p3;
    for(int k = 0; k < (L >> 3); k++)
    {
        pp = operat(g2[2*k+1], p2[2*k+1], g2[2*k], p2[2*k]);
        g3.push_back(pp.first);
        p3.push_back(pp.second);
    }

    vector<bool> g4, p4;
    for(int k = 0; k < (L >> 4); k++)
    {
        pp = operat(g3[2*k+1], p3[2*k+1], g3[2*k], p3[2*k]);
        g4.push_back(pp.first);
        p4.push_back(pp.second);
    }

    vector<bool> g5, p5;
    for(int k = 0; k < (L >> 5); k++)
    {
        pp = operat(g4[2*k+1], p4[2*k+1], g4[2*k], p4[2*k]);
        g5.push_back(pp.first);
        p5.push_back(pp.second);
    }
    bool g6 = g5[1] ^ getBitProduct(g5[0], p5[1]);
    bool c = g6;
    bool v = a[L-1] ^ c;
}

// vector<bool> PRG(vector<bool> seed)
// {
//     vector<bool> res(258);

//     return res;
// }

vector<vector<bool> > split2vec(vector<bool> src)
{
    vector<vector<bool> > res;
    vector<bool> s0(src.begin(), src.begin()+64);
    vector<bool> s1(src.begin()+64, src.begin()+128);
    vector<bool> t0(1), t1(1);
    t0[0] = src[128];
    t1[0] = src[129];
    vector<bool> v0(src.begin()+130, src.begin()+194);
    vector<bool> v1(src.begin()+194, src.begin()+258);
    res.push_back(s0);
    res.push_back(s1);
    res.push_back(t0);
    res.push_back(t1);
    res.push_back(v0);
    res.push_back(v1);
    return res;
}

void xorVectors(const std::vector<bool>& a, const std::vector<bool>& b, vector<bool>& result) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }

    std::transform(a.begin(), a.end(), b.begin(), result.begin(),
                   [](bool x, bool y) { return x ^ y; });

}

vector<bool> getvSum(const std::vector<bool>& a, const std::vector<bool>& b)
{
    //64bit
    vector<bool> res(64);
    bool value;
    ull w = 1, aval = 0, bval = 0;
    for(int i = 0; i < L; i++)
    {
        aval += (a[i] * w);
        bval += (b[i] * w);
        w *= 2;
    }
    ull result = aval + bval;
    bitset<64> ans(result);
    for(int i = 0; i < L; i++)
    {
        res[i] = ans[i];
    }
    return res;

}

void oursDecisionNodeEvl(pair<vector<bool>, vector<bool> > key, ull x)
{
    // cout <<"++++++++++"<<endl;

    bitset<L> xbit(x);
    vector<bool> s0(key.first.begin(), key.first.begin()+64);
    vector<bool> s1(key.first.begin()+64, key.first.begin()+128);
    bool t0 = key.first[128];
    bool t1 = key.first[129];
    vector<bool> v0(key.first.begin()+130, key.first.begin()+194);
    vector<bool> v1(key.first.begin()+194, key.first.begin()+258);
    
    
    vector<vector<vector<bool> > > cw(2, vector<vector<bool> >(L));
    int step = 0;
    for(int i = 0; i < 2; i++)
    {
        for(int j = 1; j < L; j++)
        {
            cw[i][j].assign(key.second.begin() + step, key.second.begin() + step + STV_SIZE);
            step += STV_SIZE;
        }
    }

    vector<bool> s = (xbit[0] == 0 ? s0: s1);
    bool t = (t == 0 ? t0 : t1);
    vector<bool> v = (xbit[0] == 0 ? v0 : v1);
    // cout << "s:" <<"\n";
    // cout << "s:" <<"\n";
    // cout << "s:" <<"\n";
    // cout << "s:" <<"\n";
    // for(int i = 0; i < 64; i++)
    // {
    //     cout << s[i] << " ";
    // }
    // cout << endl;
    // cout <<"=========" << endl;
    for(int i = 2; i < L; i++)
    {
        vector<bool> compose = PPRG(s);
        // vector<bool> compose(264);
        // cout << "comp111ose:";
        // for(int j = 0; j < 64; j++)
        // {
        //     cout <<compose[j] << " ";
        // }
        // cout << endl;
        vector<vector<bool> > gs = split2vec(compose);
        vector<vector<bool> > cwords = split2vec(cw[t][i-1]);
        xorVectors(gs[xbit[i]], cwords[t], s);
        t = gs[xbit[i] == 0 ? 2 : 3][0] ^ cwords[t == 0 ? 2 : 3][0];
        vector<bool> tmp = getvSum(gs[xbit[i] == 0 ? 4 : 5], cwords[t == 0 ? 4 : 5]);
        v = getvSum(v, tmp);
    }

}


void testDecisionNodeEvl(vector<pair<int, int> > params)
{
    std::random_device rd;  // 随机种子
    std::mt19937 gen(rd()); // Mersenne Twister 引擎
    std::uniform_int_distribution<> distrib(0, 1); // 生成 0 或 1
    ull x = get_random();
    // 258=64+64+1+1+64+64
    // 32508=(64-1)*2*258 CW_SIZE
    vector<bool> head(STV_SIZE), tail(CW_SIZE);
    for(int i = 0; i < STV_SIZE; i++)
    {
        head[i] = (distrib(gen) == 1);
    }

    for(int i = 0; i < CW_SIZE; i++)
    {
        tail[i] = (distrib(gen) == 1);
    }

    for(int i = 0; i < params.size(); i++)
    {
        ull yj = get_random();
        ull xsigm = get_random();
        // ull dn = pow(2, params[i].first-1)-1;
        ull dn = 1;
        clockrep temp;
        temp = measureTime(z22DecisionNodeEvl, yj, xsigm);
        cout << "z22DecisionNodeEvl,"<< std::fixed<<params[i].first << ","<<params[i].second << "," << temp <<", "
        << (ull)(temp + (2+3*(64-2))*trans_delay + (2+3*(64-2))*2/bandwidth) *  dn<< endl;
        temp = measureTime(z23DecisionNodeEvl,yj, xsigm);
        cout << "z23DecisionNodeEvl:" << params[i].first << ","<<params[i].second << ","<< temp <<", "<< 
        (ull)(temp + 7*trans_delay + 63*2*2/bandwidth) *  dn << endl;
        
        temp = measureTime(oursDecisionNodeEvl, make_pair(head, tail), x);
        // oursDecisionNodeEvl(make_pair(head, tail), x);
        cout << "oursDecisionNodeEvl:" <<params[i].first << ","<<params[i].second << ","<< temp <<", "<< temp *  dn << endl;
    }
    
}

void resGenz22(int d, vu dnodes, vu leaves, vu edges)
{
    int lowst_eg_st = pow(2, d-1)-2;
    int lowst_eg_ed = lowst_eg_st + pow(2, d-1);
    ull ans = 0;
    for(int j = lowst_eg_st; j < lowst_eg_ed; j++)
    {
        if(j == 0)break;
        ull res = 1;
        if(j & 1 == 0)
        {
            int k = (j-2) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-2) / 2)
            {
                res = getRingPrroduct(res, edges[k]);
            }
        }else if(j&1 != 0)
        {
            int k = (j-3) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-3) / 2)
            {
                res = getRingPrroduct(res, edges[k]);
            }
        }
        res = getRingPrroduct(res, leaves[j-lowst_eg_st]);
        ans += res;
    }
}

void resGenz23(int d, vu dnodes, vu leaves, vu edges)
{
    int lowst_eg_st = pow(2, d-1)-2;
    int lowst_eg_ed = lowst_eg_st + pow(2, d-1);
    ull ans = 0;
    for(int j = lowst_eg_st; j < lowst_eg_ed; j++)
    {
        if(j == 0)break;
        ull res = 1;
        if(j & 1 == 0)
        {
            int k = (j-2) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-2) / 2)
            {
                res = getRingPrroduct(res, edges[k]);
            }
        }else if(j&1 != 0)
        {
            int k = (j-3) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-3) / 2)
            {
                res = getRingPrroduct(res, edges[k]);
            }
        }
        res = getRingPrroduct(res, leaves[j-lowst_eg_st]);
        ans += res;
    }
}

vu permute(vu &pc){
    int len = pc.size();
    int offset = 8;
    vu npc(len);
    for(int i = 0; i < len; i++)
    {
        npc[i] = pc[(i+offset)%len];
    }
    return npc;

}


void resGenours(int d, vu dnodes, vu leaves, vu edges, vu r, vu rr)
{
    int lowst_eg_st = pow(2, d-1)-2;
    int lowst_eg_ed = lowst_eg_st + pow(2, d-1);
    ull ans = 0;
    vu pc(leaves.size());
    for(int j = lowst_eg_st; j < lowst_eg_ed; j++)
    {
        if(j == 0)break;
        ull res = 1;
        if(j & 1 == 0)
        {
            int k = (j-2) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-2) / 2)
            {
                res = res + edges[k];
            }
        }else if(j&1 != 0)
        {
            int k = (j-3) / 2;
            res = edges[j];
            for(;k >= 0; k = (k-3) / 2)
            {
                res = res + edges[k];
            }
        }

        pc[j-lowst_eg_st] = res * r[j-lowst_eg_st];
        leaves[j-lowst_eg_st] = leaves[j-lowst_eg_st] * rr[j-lowst_eg_st];
        
    }
    pc = permute(pc);
    leaves = permute(leaves);
}


void testInferenceResultGen(vector<pair<int, int> > params)
{
    for(int i = 0; i < params.size(); i++)
    {
        int d = params[i].first;
        int dn = pow(2, d-1) -1;
        vu dnodes(dn), edges1(dn*2+1),edges2(dn*2+1), leaves(dn+1), r(dn+1);

        for(int i = 0; i < dn; i++)
        {
            dnodes[i] = get_random();
            leaves[i] = get_random();
            r[i] = get_random();
            edges1[i*2] = 1-dnodes[i];
            edges1[i*2+1] = dnodes[i]; 
            edges2[i*2] = dnodes[i];
            edges2[i*2+1] = 1- dnodes[i];
        }
        leaves[dn] = get_random();
        r[dn] = get_random();
        ull ln = pow(2, d-1);
        ull tty = ln*(d-1) +1;
        clockrep temp;
        temp = measureTime(resGenz22, d, dnodes, leaves, edges1);
        cout << "z22," << params[i].first <<"," << params[i].second << "," 
        << temp + tty*trans_delay + tty*2*64/bandwidth<< endl;
        temp = measureTime(resGenz23, d, dnodes, leaves, edges1);
        cout << "z23," << params[i].first <<"," << params[i].second << "," 
        << temp + tty*trans_delay + tty*2*64/bandwidth << endl;
        temp = measureTime(resGenours, d, dnodes, leaves, edges2, r, r);
        cout << "ours," << params[i].first <<"," << params[i].second << "," 
        << temp + 2*trans_delay + (ln*4+1)/bandwidth<< endl;
        

    }
}

void testComputCost(vector<pair<int, int> > params)
{
    // testFeatureSelection(params);
    // testDecisionNodeEvl(params);
    testInferenceResultGen(params);

}

int main()
{
    // freopen("onlinecomp.txt", "w", stdout);
    vector<pair<int, int> > params;
    params.push_back(make_pair(2, 19));
    params.push_back(make_pair(3, 13));
    params.push_back(make_pair(4, 15));
    params.push_back(make_pair(8, 9));
    params.push_back(make_pair(13, 13));
    params.push_back(make_pair(17, 57));
    cout << "article,d,n,cost"<<fixed<<setprecision(0)<<endl;
    testComputCost(params);
    // vector<bool> seed(64);
    // for(int i = 0; i < 64; i++)
    // {
    //     seed[i] = (i & 1);
    // }
    // vector<bool> res = PPRG(seed);
    // cout << res.size()<<endl;
    // for(int i = 0; i < res.size(); i++)
    // {
    //     cout << res[i] << " ";
    // }
    // cout << "\n";
}