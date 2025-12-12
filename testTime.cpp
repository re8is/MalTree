#include<iostream>
#include<chrono>



int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // auto ns_since_epoch = start.time_since_epoch().count();
    // std::cout << ns_since_epoch << std::endl;
    // 待测试的代码
    // int y = fun();
    int x;
    for (int i = 0; i < 10; ++i) {
        // 模拟耗时操作
        x = i * i; 
    }
    std::cin >> x;
    // std::cout << x << std::endl;
    // [[assume(x != 0)]]
    // int y = x;
    // std::cout << x << std::endl;
    // std::chrono::nanoseconds	纳秒（1e-9 秒）
    // std::chrono::microseconds	微秒（1e-6 秒）
    // std::chrono::milliseconds	毫秒（1e-3 秒）
    // std::chrono::seconds	秒
    // std::chrono::minutes	分钟
    // std::chrono::hours	小时
    auto end = std::chrono::high_resolution_clock::now();
    // auto ee = end.time_since_epoch().count();
    // std::cout << ee << std::endl;
    // 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ns" << std::endl;
    // std::cout << "====="<< (ee -  ns_since_epoch)<< std::endl;
    
}
