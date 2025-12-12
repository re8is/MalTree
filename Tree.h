#pragma once
#include <vector>
#include <vector>
#include <chrono>
#include <iostream>
#include <cmath>
using ull = unsigned long long;
using vu = std::vector<ull>;
using vvu = std::vector<std::vector<ull> >;
using vb = std::vector<bool>;
using vvb = std::vector<std::vector<bool> >;
using pbb = std::pair<bool, bool>;
using puu = std::pair<ull, ull>;
using vpb = std::vector<pbb>;
using vpu = std::vector<puu>;
using clockrep = std::chrono::_V2::system_clock::rep;

struct Tree
{
    vu left_children;
    vu right_children;
    vu feature;
    vu threshold;
    vu impurity;
    vu n_node_samples;
    vvu value;
    int depth;
};

void testClient(int n);
void testModelProvider(int dn, int n);

template <typename Func, typename... Args>
std::chrono::_V2::system_clock::rep measureTime(Func func, Args&&... args)
{
    auto start = std::chrono::high_resolution_clock::now();

    func(std::forward<Args>(args)...);

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    auto res = duration.count();
    // std::cout << res << " " << res+1<<std::endl;
#ifdef DEBUG
    (std::cout << "args="<< ... << args) <<" Execution time: " << duration.count() << " ns" << std::endl;
#endif
    return res;
}



