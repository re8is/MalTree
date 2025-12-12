#include "Tree.h"
#include <iostream>
#include <random>
#include <limits>
#include <bitset>
#include <cassert>
// using namespace std;
Tree *t, *t0, *t1;
// bit length
const int L = 64;
// nof nodes
const int N = 9;
// nof labels
const int LENY = 3;
// nof features
const int LENF = 4;
// nof depth of tree
const int DEPTH = 3;
// random seed for dev only.
ull seed = 1LL;
const int LENBT = 1000;
// Beaver Triples
vvu bt(3), bt0(3), bt1(3);
// BTriples index
int btIdx = 0;
// bool BT index
int bbtIdx = 0;
// arithmetic index
int atIdx = 0;
//bt max value
const int BOOL_TRIPLES_MAX = 100000;
const int ARITH_TROPLES_MAX = 10000;
void fore(vu v)
{
    for (ull num : v)
    {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

void printTree(Tree *t)
{
    std::cout << "left_children: ";
    fore(t->left_children);
    std::cout << "\n"
              << "right_children: ";
    fore(t->right_children);
    std::cout << "\n"
              << "feature: ";
    fore(t->feature);
    std::cout << "\n"
              << "threshold: ";
    fore(t->threshold);
    std::cout << "\n"
              << "n_node_samples: ";
    fore(t->n_node_samples);
    std::cout << "\n"
              << "value: \n";
    for (vu v : t->value)
        fore(v);
    std::cout << "\ndepth: " << t->depth << std::endl;
}

bool lt(ull a, ull b)
{
    bool am = a >> 63;
    bool bm = b >> 63;
    return am == bm ? (a < b ? true : false) : (am > bm ? true : false);
}

void predict(Tree *t, ull s[])
{
    ull idx = 0;
    ull dep = 0;
    while (dep < t->depth)
    {
        if (lt(s[t->feature[idx]], t->threshold[idx]))
            idx = t->left_children[idx];
        else
            idx = t->right_children[idx];

        if ((t->left_children[idx] >> 63) == 1)
        {
            vu value = t->value[idx];
            fore(value);
        }
        dep++;
    }
}

void init(Tree *t)
{
    freopen("in.txt", "r", stdin);
    int n = N, tn;
    while (n--)
    {
        std::cin >> tn;
        t->left_children.push_back(tn);
    }

    n = N;
    while (n--)
    {
        std::cin >> tn;
        t->right_children.push_back(tn);
    }

    n = N;
    while (n--)
    {
        std::cin >> tn;
        t->feature.push_back(tn);
    }

    n = N;
    while (n--)
    {
        std::cin >> tn;
        t->threshold.push_back(tn);
    }

    n = N;
    while (n--)
    {
        std::cin >> tn;
        t->n_node_samples.push_back(tn);
    }

    n = N;
    while (n--)
    {
        int tm = 3;
        vu v;
        while (tm--)
        {
            std::cin >> tn;
            v.push_back(tn);
        }
        t->value.push_back(v);
    }
    t->depth = DEPTH;
    printTree(t);
}

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

ull get_random2()
{
    std::srand(std::time(0)); // 用系统时间初始化种子
    ull random_number = rand();
    return random_number;
}

void gen_tree_share(Tree *t, Tree *t0, Tree *t1)
{
    int n = N;

    for (int i = 0; i < N; i++)
    {
        t0->left_children.push_back(get_random());
        t1->left_children.push_back(t->left_children[i] - t0->left_children[i]);

        t0->right_children.push_back(get_random());
        t1->right_children.push_back(t->right_children[i] - t0->right_children[i]);

        t0->feature.push_back(get_random());
        t1->feature.push_back(t->feature[i] - t0->feature[i]);

        t0->threshold.push_back(get_random());
        t1->threshold.push_back(t->threshold[i] - t0->threshold[i]);

        t0->n_node_samples.push_back(get_random());
        t1->n_node_samples.push_back(t->n_node_samples[i] - t0->n_node_samples[i]);

        vu v, u;
        for (int j = 0; j < LENY; j++)
        {
            v.push_back(get_random());
            u.push_back(t->value[i][j] - v[j]);
        }
        t0->value.push_back(v);
        t1->value.push_back(u);
    }

    t0->depth = DEPTH;
    t1->depth = DEPTH;
    std::cout << "Random Tree t0 as follows:" << std::endl;
    printTree(t0);
    std::cout << "Random Tree t1 as follows:" << std::endl;
    printTree(t1);
}

void gen_matrix_share(vvu &m, vvu &m0, vvu &m1)
{
    int tmp;
    for (int i = 0; i < N; i++)
    {
        std::cin >> tmp;
        m[i][tmp] = 1;
        vu v1;
        vu v2;
        for (int j = 0; j < LENF; j++)
        {
            ull tp = get_random();
            v1.push_back(tp);
            v2.push_back(m[i][j] - tp);
        }
        m0.push_back(v1);
        m1.push_back(v2);
    }
}

// provider
void input_prepare(Tree *t, Tree *t0, Tree *t1, vvu &m, vvu &m0, vvu &m1)
{
    gen_tree_share(t, t0, t1);
    gen_matrix_share(m, m0, m1);
}

// client
void send_data(vu &s, vu &s0, vu &s1)
{
    int tv;

    for (int i = 0; i < LENF; i++)
    {
        std::cin >> tv;
        s.push_back(tv);
        s0.push_back(get_random());
        s1.push_back(s[i] - s0[i]);
    }
}

void gen_btriples(vvu &bt, vvu &bt0, vvu &bt1)
{
    ull a, b, c, a0, b0, c0;
    int n = LENBT;
    int cnt = 0;
    while (n--)
    {
        a = get_random();
        if (cnt++ % N == 0)
            b = get_random();
        a0 = get_random();
        b0 = get_random();
        c0 = get_random();
        c = a * b;
        bt[0].push_back(a);
        bt[1].push_back(b);
        bt[2].push_back(c);

        bt0[0].push_back(a0);
        bt0[1].push_back(b0);
        bt0[2].push_back(c0);

        bt1[0].push_back(a - a0);
        bt1[1].push_back(b - b0);
        bt1[2].push_back(c - c0);
    }
}

void gen_tripl_matrix(vvu &G1, vvu &G10, vvu &G11,
                      vu &g2, vu &g20, vu &g21,
                      vu &g3, vu &g30, vu &g31)
{
    int cnt = 0;
    // G1
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < LENF; j++)
        {
            G1[i][j] = get_random();
            G10[i][j] = get_random();
            G11[i][j] = G1[i][j] - G10[i][j];
        }
    }
    // g2
    for (int i = 0; i < LENF; i++)
    {
        g2[i] = get_random();
        g20[i] = get_random();
        g21[i] = g2[i] - g20[i];
    }
    // g3
    for (int i = 0; i < N; i++)
    {
        ull tmp = 0;
        for (int j = 0; j < LENF; j++)
        {
            tmp += G1[i][j] * g2[j];
        }
        g3[i] = tmp;
        g30[i] = get_random();
        g31[i] = tmp - g30[i];
    }
    // std::cout << G
}
// 1 2 3     7
// 4 5 6     5
// 5 2 1     4
// 6 5 4

vvu eval_e(vvu m, vvu G1)
{
    vvu e;
    for (int i = 0; i < N; i++)
    {
        vu tmp;
        for (int j = 0; j < LENF; j++)
        {
            tmp.push_back(m[i][j] - G1[i][j]);
        }
        e.push_back(tmp);
    }
    return e;
}

vu eval_f(vu x, vu g2)
{
    vu f;
    for (int i = 0; i < LENF; i++)
    {
        f.push_back(x[i] - g2[i]);
    }
    return f;
}

vu get_mult(vvu x, vu y)
{
    vu res;
    for (int i = 0; i < x.size(); i++)
    {
        ull tmp = 0;
        for (int j = 0; j < y.size(); j++)
        {
            tmp += x[i][j] * y[j];
        }
        res.push_back(tmp);
    }
    return res;
}

vu get_add(vu x, vu y)
{
    vu res;
    for (int i = 0; i < x.size(); i++)
    {
        res.push_back(x[i] + y[i]);
    }
    return res;
}

vu eval_Mx(int id, vvu E, vu f, vvu G1, vu g2, vu g3)
{
    vu eg = get_mult(E, g2);
    vu gf = get_mult(G1, f);
    vu res = get_add(get_add(eg, gf), g3);
    if (id == 1)
    {
        res = get_add(res, get_mult(E, f));
    }

    return res;
}

// parties' m, s, bts
void input_select(vvu m0, vvu m1,
                  vu s0, vu s1,
                  vu &mx0, vu &mx1)
{
    vvu G1(N, vu(LENF, 0)), G10(N, vu(LENF, 0)), G11(N, vu(LENF, 0));
    vu g2(LENF), g20(LENF), g21(LENF);
    vu g3(N), g30(N), g31(N);
    gen_tripl_matrix(G1, G10, G11, g2, g20, g21, g3, g30, g31);
    vvu E0 = eval_e(m0, G10);
    vvu E1 = eval_e(m1, G11);
    vu f0 = eval_f(s0, g20);
    vu f1 = eval_f(s1, g21);
    vvu E;
    for (int i = 0; i < E0.size(); i++)
    {
        vu tmp;
        for (int j = 0; j < E0[i].size(); j++)
        {
            tmp.push_back(E0[i][j] + E1[i][j]);
        }
        E.push_back(tmp);
    }

    vu f;
    for (int i = 0; i < f0.size(); i++)
    {
        f.push_back(f0[i] + f1[i]);
    }
    mx0 = eval_Mx(0, E, f, G10, g20, g30);
    mx1 = eval_Mx(1, E, f, G11, g21, g31);

    // vu real = get_add(mx0, mx1);
    // vu target = get_mult(m, get_add(s0, s1));
}

void gen_bool_bts(vvb &bbt, vvb &bbt0, vvb &bbt1)
{

    vb t1, t2, t3, t4, t5, t6, t7, t8, t9;
    for (int i = 0; i < BOOL_TRIPLES_MAX; i++)
    {
        bool b1 = get_random() % 2;
        bool b2 = get_random() % 2;
        bool b3 = get_random() % 2;
        t1.push_back(b1);
        t2.push_back(b2);
        t3.push_back(b1 & b2);
        b1 = get_random() % 2;
        b2 = get_random() % 2;
        t4.push_back(b1);
        t5.push_back(b2);
        t6.push_back(b3);
        t7.push_back(t4[i] ^ t1[i]);
        t8.push_back(t5[i] ^ t2[i]);
        t9.push_back(t6[i] ^ t3[i]);
    }
    bbt.push_back(t1);
    bbt.push_back(t2);
    bbt.push_back(t3);
    bbt0.push_back(t4);
    bbt0.push_back(t5);
    bbt0.push_back(t6);
    bbt1.push_back(t7);
    bbt1.push_back(t8);
    bbt1.push_back(t9);
}
vvb bbt, bbt0, bbt1;
pbb mult_bool(bool x0, bool x1, bool y0, bool y1)
{
    assert(bbtIdx < BOOL_TRIPLES_MAX);
    bool e = (x0 ^ bbt0[0][bbtIdx]) ^ (x1 ^ bbt1[0][bbtIdx]);
    bool f = (y0 ^ bbt0[1][bbtIdx]) ^ (y1 ^ bbt1[1][bbtIdx]);
    bool p1 = e & f ^ bbt0[1][bbtIdx] & e ^ bbt0[0][bbtIdx] & f ^ bbt0[2][bbtIdx];
    bool p2 = bbt1[1][bbtIdx] & e ^ bbt1[0][bbtIdx] & f ^ bbt1[2][bbtIdx];
    bbtIdx++;
    pbb res(p1, p2);
    return res;
}

// process rhombus operator
std::pair<pbb, pbb> up_op(std::pair<pbb, pbb> gp2, std::pair<pbb, pbb> gp1)
{
    pbb g1p2 = mult_bool(gp1.first.first, gp1.first.second, gp2.second.first, gp2.second.second);
    pbb g_new(gp2.first.first ^ g1p2.first, gp2.first.second ^ g1p2.second);
    pbb p_new = mult_bool(gp2.second.first, gp2.second.second, gp1.second.first, gp1.second.second);
    return std::make_pair(g_new, p_new);
}

void get_MSB(ull num)
{

    // gen_bool_bts(bbt, bbt0, bbt1);
    pbb p;
}

pbb dnode_eval(ull x0, ull x1, ull y0, ull y1)
{
    ull delta0 = x0 - y0;
    ull delta1 = x1 - y1;

    std::bitset<64> a(delta0);
    std::bitset<64> b(delta1);

    vb w0, w1;
    std::vector<pbb> w;
    std::vector<pbb> g;
    std::vector<pbb> p;
    std::vector<std::pair<pbb, pbb>> node;

    for (int i = 0; i < L; i++)
    {
        // process w
        pbb wi;
        wi.first = a[i];
        wi.second = b[i];
        w.push_back(wi);

        // process g p
        pbb gi = mult_bool(a[i], 0, 0, b[i]);
        pbb pi(wi);
        node.push_back(std::make_pair(gi, pi));
    }
#ifdef DEBUG
    //
    std::cout << "x0 x1 y0 y1:" << x0 << " " << x1 << " " << y0 << " " << y1 << std::endl;
    std::cout << "a:";
    for (int k = L - 1; k >= 0; k--)
        std::cout << w[k].first;
    std::cout << std::endl;
    std::cout << "b:";
    for (int k = L - 1; k >= 0; k--)
        std::cout << w[k].second;
    std::cout << std::endl;

    std::cout << "G:";
    for (int k = L - 1; k >= 0; k--)
    {
        bool tg = node[k].first.first ^ node[k].first.second;
        // bool tp = node[k].second.first ^ node[k].second.second;
        std::cout << tg;
    }
    std::cout << "\nP:";
    for (int k = L - 1; k >= 0; k--)
    {
        // bool tg = node[k].first.first ^ node[k].first.second;
        bool tp = node[k].second.first ^ node[k].second.second;
        std::cout << tp;
    }
    std::cout << std::endl;

    std::cout << "-------------------" << std::endl;
#endif
    for (int curL = L >> 1; curL > 0; curL >>= 1)
    {
        if (curL == 32)
            curL--;
        for (int i = 0; i < curL; i++)
        {
            node[i] = up_op(node[2 * i + 1], node[2 * i]);
        }
        if (curL == 31)
        {
            node[curL] = node[curL * 2];
            curL++;
        }
#ifdef DEBUG
        std::cout << curL << "G: ";
        for (int k = curL - 1; k >= 0; k--)
        {
            bool tg = node[k].first.first ^ node[k].first.second;
            // bool tp = node[k].second.first ^ node[k].second.second;
            std::cout << tg;
        }
        std::cout << "\nP:";
        for (int k = curL - 1; k >= 0; k--)
        {
            // bool tg = node[k].first.first ^ node[k].first.second;
            bool tp = node[k].second.first ^ node[k].second.second;
            std::cout << tp;
        }
        std::cout << std::endl;

        std::cout << "-------------------" << std::endl;
#endif
    }

    return std::make_pair(node[0].first.first ^ w[L - 1].first, node[0].first.second ^ w[L - 1].second);
}

bool judge_bound(ull a, ull b)
{
    ull lower_bound = 1UL << (L - 2);
    ull upper_bound = (-1) - lower_bound + 1;
    ull sum = a + b;
    return !(sum >= lower_bound && sum < upper_bound);
}

vvu atriples0, atriples1, atriples;
puu mult_arth(ull x0, ull x1, ull y0, ull y1)
{
    ull e = x0 - atriples0[0][atIdx] + x1 - atriples1[0][atIdx];
    ull f = y0 - atriples0[1][atIdx] + y1 - atriples1[1][atIdx];
    ull p1 = e*f + atriples0[1][atIdx]*e + atriples0[0][atIdx]*f + atriples0[2][atIdx];
    ull p2 = atriples1[1][atIdx]*e + atriples1[0][atIdx]*f + atriples1[2][atIdx];
    return std::make_pair(p1, p2);
}

vpu b2a(vpb list)
{
    vpu list_a;
    for(int i = 0; i < N; i++)
    {
        bool t1 = list[i].first;
        bool t2 = list[i].second;
        mult_bool(list[i].first, 0, 0, list[i].second);

    }
    return list_a;
}

void gen_arith_triples()
{
    
}

puu sec_inference_gen(vpb list, vpu leaves)
{
    vpu cr = b2a(list);
}

int main()
{
    // Fake Model t
    // Tree tree;
    // t = &tree;
    // init(t);
    // Tree tree0;
    // t0 = &tree0;
    // Tree tree1;
    // t1 = &tree1;

    // // matrix
    // vvu m(N, vu(LENF, 0)), m0, m1;

    // input_prepare(t, t0, t1, m, m0, m1);

    // // Fake x
    // vu s, s0, s1;
    // send_data(s, s0, s1);

    // // Beaver Triples

    // gen_btriples(bt, bt0, bt1);

    // vvu e(N, vu(LENF, 0)), e0(N, vu(LENF, 0)), e1(N, vu(LENF, 0));
    // vu f(LENF), f0(LENF), f1(LENF);
    // vu mx0, mx1;
    // input_select(m0, m1, s0, s1, mx0, mx1);

    // // pbb res = dnode_eval(20, 3,10, 5);
    // gen_bool_bts(bbt, bbt0, bbt1);
    // vvu v(5, vu(200, 0));
    // ull mm = 6358486341003018661;
    // int ccnt = 0;
    // for (int i = 0; i < 200; i++)
    // {
    //     do
    //     {
    //         v[0][i] = get_random();
    //         v[1][i] = get_random();
    //         v[2][i] = get_random();
    //         v[3][i] = get_random();
    //     } while (!(judge_bound(v[0][i], v[1][i]) && judge_bound(v[2][i], v[3][i])));
    //     // v[0][i] = get_random() % mm;
    //     // v[1][i] = get_random() % mm;
    //     // v[2][i] = get_random() % mm;
    //     // v[3][i] = get_random() % mm;

    //     ull curr1 = v[0][i] + v[1][i];
    //     ull curr2 = curr1 - v[2][i];
    //     ull curr3 = curr2 - v[3][i];
    //     std::bitset<64> rre(curr3);
    //     v[4][i] = curr3 >> 63;
    //     pbb rr = dnode_eval(v[0][i], v[1][i], v[2][i], v[3][i]);
    //     // std::cout << (((rr.first ^ rr.second) == v[4][i]) ? true : false) << std::endl;
    //     if ((rr.first ^ rr.second) != v[4][i])
    //     {
    //         std::cout << ++ccnt << "-" << i << "-" << v[0][i] << " " << v[1][i] << " " << v[2][i] << " " << v[3][i] << std::endl;
    //     }
    // }

    // pbb rr = dnode_eval(12, 11, 19, 7);

    // puu res = sec_inference_gen();
    std::vector<clockrep >s(5);
    std::cout << s[0] <<s[1] <<s[2] <<s[3] <<s[4] <<std::endl;
    for(int i = 0; i < 10; i++)
    {
        measureTime(testClient, 1000);
        clockrep a15 = measureTime(testClient, 15);
        s[0] += a15;
        clockrep a131 = measureTime(testClient, 13);
        s[1] += a131;
        clockrep a9 = measureTime(testClient, 9);
        s[2] += a9;
        clockrep a132 = measureTime(testClient, 13);
        s[3] += a132;
        clockrep a57 = measureTime(testClient, 57);
        s[4] += a57;
    }
    for(int i = 0; i < 5; i++)std::cout << s[i] /10<< std::endl;
    // for(int i = 0; i < 1; i++)
    // {
    //     measureTime(testModelProvider, 10, 20);
    //     clockrep a15 = measureTime(testModelProvider, 7, 15);
    //     s[0] += a15;
    //     clockrep a131 = measureTime(testModelProvider, 3, 13);
    //     s[1] += a131;
    //     clockrep a9 = measureTime(testModelProvider, 127, 9);
    //     s[2] += a9;
    //     clockrep a132 = measureTime(testModelProvider, 4095, 13);
    //     s[3] += a132;
    //     clockrep a57 = measureTime(testModelProvider, 65535, 57);
    //     s[4] += a57;
    // }
    // for(int i = 0; i < 5; i++)std::cout << s[i] /1<< std::endl;
    // std::cout << get_random2() << std::endl;
    return 0;
}

void testClient(int n)
{ 
    // 1024  15
    // 966   13
    // 885   9 
    // 970   13
    // 2183  57
    vu x;
    for(int i = 0; i < n; i++) 
    {
        x.push_back(1);
        get_random();
        x[i] = x[i]+i;
    }
}
// 3   13   3
// 4   15   7
// 8   9    127
// 13  13   4095
// 17  57   65535
void testModelProvider(int dn, int n)
{
    vvu sm1(dn, vu(n, 0));
    vvu sm2(dn, vu(n, 0));
    for(int i = 0; i < dn; i++)
    {
        for(int j = 0; j < n; j++)
        {
            sm1[i].push_back(j);
            sm1[i][j] += get_random2();
            sm2[i][j] = sm1[i][j] - i;
        }
    }


}
