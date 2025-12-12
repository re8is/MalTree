#include <iostream>
#include <typeinfo>
#include <vector>
#include <bitset>
#include <random>
#include <limits>
#include <cassert>
#include <utility>
#include <iomanip>
#include <format>
#include <openssl/bn.h>
#include "Tree.h"
#include<cstring>
// using ull = unsigned long long;
using namespace std;
int a;
void gen_matrix_share(std::vector<std::vector<ull>> &m)
{
    int tmp;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            ull tp = 12345678910111213;
            m.at(i).push_back(tp);
        }
    }
}

void fun(bool a, bool b, bool c)
{
    bool res1 = a ^ b & c;
    bool res2 = (a ^ b) & c;
    std::cout << res1 << res2 << std::endl;
}

pbb mult_bool(bool x0, bool x1, bool y0, bool y1)
{

    return std::make_pair(0, 1);
}

ull dnode_eval(ull y0, ull y1, ull x0, ull x1)
{
    ull delta0 = x0 - y0;
    ull delta1 = x1 - y1;

    std::bitset<64> a(delta0);
    std::bitset<64> b(delta1);
    pbb gi = mult_bool(a[0], 0, b[0], 0);
    return 1UL;
}

ull seed = 2;
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

int main()
{
    // ull x = 15373252921830533746;
    // ull y = 2203844724969914387;

    // ull a = 1UL<<(64-2);
    // ull b = (-1) - a + 1;

    // bool a = get_random() & 1;
    // bool b = get_random() & 1;
    // bool c = a & b;

    // bool a1 = get_random() & 1;
    // bool a2 = a ^ a1;

    // bool b1 = get_random() & 1;
    // bool b2 = b ^ b1;

    // bool c1 = get_random() & 1;
    // bool c2 = c ^ c1;

    // double a = 1.23;
    // double b = 4.56;
    // cout << std::format("{:.0f}", a) << endl;
    // cout << b << endl;

    // BIGNUM* bn = BN_new();

    

    // BN_free(bn);
    
    // vector<bool> bb = {1,0,1,1,0};

    // vector<bool> c;
    // c.assign(bb.begin(), bb.begin()+3);

    // for(int i = 0; i < c.size(); i++)
    // {
    //     cout << c[i] << endl;
    // }
    // unsigned char t[100],key[20], iv[30];
    // for(int i = 0; i < 100; i++)t[i] = i;

    // memcpy(key, t, 20);
	// memcpy(iv, t+20, 30);
    // for(int i = 0; i < 10; i++)
    // {
    //     cout << (int) key[i] << " ";
    // }
    // cout << endl;
    // for(int i = 0; i < 10; i++)
    // {
    //     cout << (int) iv[i] << " ";
    // }
    // return 0;
    // cout << endl;

    ull a = get_random();
    ull b = get_random();
    ull c = a * b;

    ull a1 = get_random();
    ull a2 = a - a1;

    ull b1 = get_random();
    ull b2 = b - b1;

    ull c1 = get_random();
    ull c2 = c - c1;
    cout <<"ull a1="<< a1 << ",a2=" << a2 << "\n"
         <<"ull b1="<<  b1 << ",b2=" << b2 << "\n"
         <<"ull c1="<<  c1 << ",c2=" << c2 << "\n";


}
