#include "Tree.h"
#include <iostream>
#include <random>
#include <limits>
#include <bitset>
#include <cassert>
#include <utility>
using namespace std;

// random seed for dev only.
ull seed = 1LL;
const int TOTAL = 10;
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

void client_mine(vu features ,vu masks)
{
    
    for(int j = 0; j < features.size(); j++){
        // m_features[i].push_back(get_random());
        masks.push_back(get_random());
    }

    vu result;
    for(int i = 0; i < features.size(); i++)
    {
        result.push_back(features[i] + masks[i]);
    }
}

void client_z22(vu features ,vu masks)
{
    for(int j = 0; j < features.size(); j++){
        // m_features[i].push_back(get_random());
        masks.push_back(get_random());
    }
    vu result;
    for(int i = 0; i < features.size(); i++)
    {
        result.push_back(features[i] - masks[i]);
    }
}

void clientFeatureEncZ22(vu features, vu masks)
{
    vu result;
    ull rv = get_random();
    for(int i = 0; i < features.size(); i++){
        // masks.push_back(get_random());
        // masks.push_back(rv);
        result.push_back(features[i] - masks[i]);
    }
}

void clientFeatureEncZ23(vu features, vu masks)
{
    vu result;
    ull rv = get_random();
    for(int i = 0; i < features.size(); i++){
        // masks.push_back(get_random());
        // masks.push_back(rv);
        result.push_back(features[i] - masks[i]);
    }
}

void clientFeatureEncOurs(vu features, vu masks)
{
    vu result;
    ull rv = get_random();
    for(int i = 0; i < features.size(); i++){
        // masks.push_back(get_random());
        // masks.push_back(rv);
        result.push_back(features[i] - masks[i]);
    }
}

void client_z23(vu features ,vu masks)
{
    for(int j = 0; j < features.size(); j++){
        masks.push_back(get_random());
    }
    vu result;
    for(int i = 0; i < features.size(); i++)
    {
        result.push_back(features[i] - masks[i]);
    }
}

void rec_z22(int d, vu pc0, vu pc1, vu v0, vu v1)
{
    int leaves = pow(2, d-1);
    vu res;
    for(int i = 0; i < pc0.size(); i++)
    {
        res.push_back(pc0[i]+pc1[i]);
    }
    ull ress = v0[0] + v1[0];
}
void rec_z23(ull x0, ull x1)
{
    ull res = x0 + x1;
}
void recOurs(ull x0, ull x1)
{
    ull res = x0 + x1;
}

void clientFeatureEnc(vector<pair<int, int> > params)
{
    for(int i = 0; i < params.size(); i++)
    {
        vu features, masks;
        for(int j = 0; j < params[i].second; j++)
        {
            features.push_back(get_random());
            masks.push_back(get_random());
        }

        ull x0 = get_random();
        ull x1 = get_random();
        vu pc0, pc1, v0, v1;
        for(int j = 0; j < pow(2, params[i].first-1); j++)
        {
            pc0.push_back(get_random());
            pc1.push_back(get_random());
            v0.push_back(get_random());
            v1.push_back(get_random());
        }

        clockrep temp_enc = 0, temp_dec = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp_enc += measureTime(clientFeatureEncZ22, features, masks);
            temp_dec += measureTime(rec_z22, params[i].first, pc0, pc1, v0, v1);
        }
        cout << "z22," << params[i].first << "," << params[i].second << "," << temp_enc / TOTAL <<"," << temp_dec / TOTAL << endl;

        temp_enc = 0, temp_dec = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp_enc += measureTime(clientFeatureEncZ23, features, masks);
            temp_dec += measureTime(rec_z23, x0, x1);
        }
        cout << "z23," << params[i].first << "," << params[i].second << "," << temp_enc / TOTAL <<"," << temp_dec / TOTAL << endl;

        temp_enc = 0, temp_dec = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp_enc += measureTime(clientFeatureEncOurs, features, masks);
            temp_dec += measureTime(recOurs, x0, x1);
        }
        cout << "Ours," << params[i].first << "," << params[i].second << "," << temp_enc / TOTAL <<"," << temp_dec / TOTAL << endl;
    }
}


void testClient(vector<pair<int, int> > params)
{

    clientFeatureEnc(params);
    // clientResDec(params);
}

void mpMatrixEncZ22(vvu matrix)
{
    int m = matrix.size();
    int n = matrix[0].size();
    vvu masks(m, vu(n, 0)),result(m, vu(n, 0)), result2;
    for(int i = 0; i < m; i++)
    {
        // get_random();
        vu tmp;
        for(int j = 0; j < n; j++)
        {
            // masks[i][j] = get_random();
            // masks[i][j] = j;
            result[i][j] = matrix[i][j] - masks[i][j];
            tmp.push_back(result[i][j]);
            // result2[i][j] = masks[i][j];
        }
        result2.push_back(tmp);
    }
}

void mpIndexArray(vu idx)
{
    vu masks, result, result2;
    for(int i = 0; i < idx.size(); i++)
    {
        // get_random();
        // masks.push_back(get_random());
        masks.push_back(i);
        result.push_back(idx[i] - masks[i]);
        result2.push_back(masks[i]);

    }
}

void mpMatrixEncOurs(vvu matrix)
{
    int m = matrix.size();
    int n = matrix[0].size();
    vvu masks(m, vu(n, 0)),result(m, vu(n, 0)), result2;
    for(int i = 0; i < m; i++)
    {
        // get_random();
        vu tmp;
        for(int j = 0; j < n; j++)
        {
            // masks[i][j] = get_random();
            masks[i][j] = j;
            result[i][j] = matrix[i][j] + masks[i][j];
            tmp.push_back(result[i][j]);
            // result2[i][j] = result[i][j];
        }
        result2.push_back(tmp);
    }
}

void mpNodeEncZ22(vu decision_nodes, vu leaf_nodes)
{
    vu d_masks, l_masks, decision_res, leaf_res;
    int i = 0;
    for(; i < decision_nodes.size(); i++)
    {
        // d_masks.push_back(get_random());
        // l_masks.push_back(get_random());
        d_masks.push_back(i);
        l_masks.push_back(i);
        decision_res.push_back(decision_nodes[i] - d_masks[i]);
        leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
    }
    l_masks.push_back(get_random());
    leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
}

void mpNodeEncZ23(vu decision_nodes, vu leaf_nodes)
{
    vu d_masks, l_masks, decision_res, leaf_res;
    int i = 0;
    for(; i < decision_nodes.size(); i++)
    {
        // d_masks.push_back(get_random());
        // l_masks.push_back(get_random());
        d_masks.push_back(i);
        l_masks.push_back(i);
        decision_res.push_back(decision_nodes[i] - d_masks[i]);
        leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
    }
    l_masks.push_back(get_random());
    leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
}

void mpNodeEncOurs(vu decision_nodes, vu leaf_nodes)
{
    vu d_masks, l_masks, decision_res, leaf_res;
    int i = 0;
    for(; i < decision_nodes.size(); i++)
    {
        // d_masks.push_back(i);
        // l_masks.push_back(i);
        d_masks.push_back(i);
        l_masks.push_back(i);
        decision_res.push_back(decision_nodes[i] - d_masks[i]);
        leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
    }
    l_masks.push_back(i);
    leaf_res.push_back(leaf_nodes[i] - l_masks[i]);
}



void testMatrixEnc(vector<pair<int, int> > params)
{
    clockrep temp = 0;
    for(int i = 0; i < params.size(); i++)
    {
        int m = pow(2, params[i].first-1) - 1;
        int n = params[i].second;

        vector<vector<ull> > matrix, matrix_masks;
        for(int j = 0; j < m; j++)
        {
            vu t0;
            for(int k = 0; k < n; k++)
            {
                t0.push_back(get_random() % 2);
            }
            matrix.push_back(t0);
        }
        temp = 0;
        temp = measureTime(mpMatrixEncZ22, matrix);

        int leaves = pow(2, params[i].first-1);
        int dn = leaves - 1;
        vu decision_nodes, leaf_nodes;
        for(int j = 0; j < dn; j++)
        {
            decision_nodes.push_back(get_random());
            leaf_nodes.push_back(get_random());
        }
        leaf_nodes.push_back(get_random());
        clockrep tNodeEnc = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            tNodeEnc += measureTime(mpNodeEncZ22, decision_nodes, leaf_nodes);
        }
        cout << "z22," << params[i].first  <<"," << params[i].second << "," << temp <<  "," << tNodeEnc / TOTAL << endl;
        
        
        vu idx;
        for(int j = 0; j < m; j++) idx.push_back(i);
        temp = 0,tNodeEnc = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp += measureTime(mpIndexArray, idx);
            tNodeEnc += measureTime(mpNodeEncZ23, decision_nodes, leaf_nodes);
        }
        
        cout << "z23," << params[i].first  <<"," << params[i].second << "," << temp / TOTAL << "," << tNodeEnc / TOTAL << endl;
        
        temp = measureTime(mpMatrixEncOurs, matrix);
        tNodeEnc = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            tNodeEnc += measureTime(mpNodeEncOurs, decision_nodes, leaf_nodes);
        }
        cout << "Ours," << params[i].first  <<"," << params[i].second << "," << temp << "," << tNodeEnc / TOTAL << endl;
        
        
    }
}






void testNodeEnc(vector<pair<int, int> > params)
{
    clockrep temp = 0;
    for(int i = 0; i < params.size(); i++)
    {
        int leaves = pow(2, params[i].first-1);
        int dn = leaves - 1;

        vu decision_nodes, d_masks;
        for(int j = 0; j < dn; j++)
        {
            decision_nodes.push_back(get_random());
            // d_masks.push_back(get_random());
        }

        vu leaf_nodes, l_masks;
        for(int j = 0; j < dn; j++)
        {
            leaf_nodes.push_back(get_random());
            // l_masks.push_back(get_random());
        }
        temp = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp += measureTime(mpNodeEncZ22, decision_nodes, leaf_nodes);
        }
        cout << "Z22: d=" << params[i].first << " NodeEncTime(ns): " << temp / TOTAL << endl;
        temp = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp += measureTime(mpNodeEncOurs, decision_nodes, leaf_nodes);
        }
        cout << "Mine: d=" << params[i].first << " NodeEncTime(ns): " << temp / TOTAL << endl;
        temp = 0;
        for(int j = 0; j < TOTAL; j++)
        {
            temp += measureTime(mpNodeEncZ23, decision_nodes, leaf_nodes);
        }
        cout << "Mine: d=" << params[i].first << " NodeEncTime(ns): " << temp / TOTAL << endl;
        
    }
    
}


void testMP(vector<pair<int, int> > params)
{
    testMatrixEnc(params);
    // testNodeEnc(params);
}



int main()
{
    freopen("servers_fs.txt", "w", stdout);
    //first: d, second: n;
    vector<pair<int, int> > params;
    params.push_back(make_pair(2, 19));
    params.push_back(make_pair(3, 13));
    params.push_back(make_pair(4, 15));
    params.push_back(make_pair(8, 9));
    params.push_back(make_pair(13, 13));
    params.push_back(make_pair(17, 57));
    cout << "article,d,n,mat_enc_time,node_enc_time" << endl;
    // testClient(params);
    // testMP(params);




}

