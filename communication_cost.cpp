#include "Tree.h"
#include <iostream>
#include <random>
#include <limits>
#include <bitset>
#include <cassert>
#include <utility>
using namespace std;

const double word = 8.0;
const double mb = 1024*1.0;

void testFeatureSelect(vector<pair<int, int> > params)
{
    for(int i = 0; i < params.size(); i++)
    {
        int m = pow(2, params[i].first-1) - 1;
        int n = params[i].second;
        cout << "z22," <<params[i].first << "," << params[i].second << "," << (m*n+n) * word / mb << endl;
        
        cout << "z23," <<params[i].first << "," << params[i].second << ","  << 2*(n+1)*m * word / mb << endl;

        cout << "ours," <<params[i].first << "," << params[i].second << ","  << 0 << endl;
    }
}

void testDecisionNodeEvl(vector<pair<int, int> > params)
{
    for(int i = 0; i < params.size(); i++)
    {
        double tz22 = 0.00006008148193359375 * 1024;
        double tz23 = 0.000030040740966796875 * 1024;
        int m = pow(2, params[i].first-1) - 1;

        cout << "z22," <<params[i].first << "," << params[i].second << "," << tz22 * m <<endl;
        cout << "z23," <<params[i].first << "," << params[i].second << "," << tz23 * m <<endl;
        cout << "ours," <<params[i].first << "," << params[i].second << "," << 0 <<endl;
    }
}

void testInferenceGen(vector<pair<int, int> > params)
{
    for(int i = 0; i < params.size(); i++)
    {
        int leaves = pow(2, params[i].first-1);
        int n = params[i].second;
        cout << "z22," <<params[i].first << "," << params[i].second << "," << leaves*2 * word / mb<<endl;
        cout << "z23," <<params[i].first << "," << params[i].second << "," << (2+params[i].first*leaves*2+1) * word / mb<<endl;
        cout << "ours," <<params[i].first << "," << params[i].second << "," << (2+params[i].first*leaves*2+1) * word / mb<<endl;
    }
    
}

void testCommCost(vector<pair<int, int> > params)
{
    cout<<"===========testFeatureSelect==========="<<endl;
    testFeatureSelect(params);
    cout<<"===========testDecisionNodeEvl==========="<<endl;
    testDecisionNodeEvl(params);
    cout<<"===========testInferenceGen==========="<<endl;
    testInferenceGen(params);
}

int main()
{
    vector<pair<int, int> > params;
    params.push_back(make_pair(2, 19));
    params.push_back(make_pair(3, 13));
    params.push_back(make_pair(4, 15));
    params.push_back(make_pair(8, 9));
    params.push_back(make_pair(13, 13));
    params.push_back(make_pair(17, 57));
    cout << "article,d,n,cost"<<endl;
    testCommCost(params);
}