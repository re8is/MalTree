#include "Tree.h"
#include<iostream>
#include<gmp.h>
#include<ctime>
using namespace std;

void fore(vector<int> v)
{
    for(int num : v){
        std::cout << num << " ";
    }
    cout << "\n";
}

void printTree(Tree t) {
    std::cout << "left_children: ";
    fore(t.left_children);
    std::cout << "\n" << "right_children: ";
    fore(t.right_children);
    std::cout << "\n" << "feature: ";
    fore(t.feature);
    std::cout << "\n" << "threshold: ";
    fore(t.threshold);
    std::cout << "\n" << "n_node_samples: ";
    fore(t.n_node_samples);
    std::cout << "\n" << "value: \n";
    for(vector<int> v : t.value) fore(v);
    std::cout << "\ndepth: " << t.depth << std::endl;

}

void predict(Tree t, int s[])
{
    int idx = 0;
    int dep = 0;
    while (dep < t.depth)
    {
        if(s[t.feature[idx]] <= t.threshold[idx])
            idx = t.left_children[idx];
        else
            idx = t.right_children[idx];
        
        if(t.left_children[idx] < 0)
        {
            vector<int> value = t.value[idx];
            fore(value);
        }
        dep++;
    }
    
}

//provider: Input Preparation
// struct Tree* prepare()
// {

// }

void maincode()
{
    //Fake Model t
    Tree t;
    t.left_children = {1, -1, 3, 4, -1, -1,  7, -1, -1};
    t.right_children = {2, -1, 6, 5, -1, -1,  8, -1, -1};
    t.feature = {3, -2,  3,  2, -2, -2,  2, -2, -2};
    //times 100
    t.threshold={80, -200, 175, 495, -200, -200, 485, -200, -200};
    t.n_node_samples = {150, 50, 100, 54, 48, 6, 46, 3, 43};
    t.depth = 3;
    t.value = {{33,33,33},
                {100,0,0},
                {0,50,50},
                {0,91,9},
                {0,98,2},
                {0,33,67},
                {0,2,98},
                {0,33,67},
                {0,0,100}};

    printTree(t);

    int s[4];
    while (cin >> s[0] >> s[1] >> s[2] >> s[3])
    {
        predict(t, s);
    }
}

void testRandom()
{
     // 初始化大整数类型
    mpz_t upper_bound, random_num;
    mpz_init(upper_bound);
    mpz_init(random_num);

    // 设置上限值（比如 1000）
    mpz_set_ui(upper_bound, 1000);  // 上限为 1000

    // 获取当前时间作为种子
    gmp_randstate_t state;
    gmp_randinit_mt(state);  // 初始化随机数生成器为 Mersenne Twister（MT）
    gmp_randseed_ui(state, static_cast<unsigned int>(time(nullptr)));  // 设置随机数种子

    // 生成一个在 [0, 1000) 范围内的随机整数
    mpz_urandomm(random_num, state, upper_bound);

    // 输出随机整数
    std::cout << "生成的随机整数：";
    mpz_out_str(stdout, 10, random_num);  // 输出到控制台，10表示十进制
    std::cout << std::endl;

    // 清理资源
    mpz_clear(upper_bound);
    mpz_clear(random_num);
    gmp_randclear(state);
}

void testASMDM(){
    mpz_t a, b, mod, res;
    mpz_init(res);
    mpz_init(a);
    mpz_set_str(a, "1234", 10);
    // mpz_set_si(a, 12345);
    mpz_out_str(stdout, 10, a);
    cout << endl;
    // cout << mpz_out_str(stdout, 10, b) << endl;
    // cout << mpz_out_str(stdout, 10, res) << endl;

    // mpz_add(res, a, b);
    // cout << mpz_out_str(stdout, 10, res) << endl;

    // mpz_sub(res, a, b);
    // cout << mpz_out_str(stdout, 10, res) << endl;

    // mpz_mul(res, a, b);
    // cout << mpz_out_str(stdout, 10, res) << endl;

    // mpz_div(res, a, b);
    // cout << mpz_out_str(stdout, 10, res) << endl;

    // mpz_mod(res, a, b);
    // cout << mpz_out_str(stdout, 10, res) << endl;
    mpz_clear(a);
    // mpz_clear(b);
    // mpz_clear(mod);
    // mpz_clear(res);
    
}

int main()
{
    // // maincode();
    // mpz_t upper_bound, random_num;
    // // mpz_inits(upper_bound, random_num);
    // mpz_init(upper_bound);
    // mpz_init(random_num);
    // // mpz_set_ui(a, 7);
    
    // mpz_set_ui(upper_bound, 1000);
    
    // // mpz_out_str(stdout, 10, a);

    // gmp_randstate_t state;
    // gmp_randinit_mt(state);
    // gmp_randseed_ui(state, static_cast<unsigned int>(time(nullptr)));

    // mpz_urandomm(random_num, state, upper_bound);

    // mpz_out_str(stdout, 10, random_num);
    // std::cout << std::endl;

    // mpz_clears(upper_bound, random_num);
    // // mpz_clear(upper_bound);
    // // mpz_clear(random_num);
    // gmp_randclear(state);
    //testRandom();
    testASMDM();
}
