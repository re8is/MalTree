#include "Utils.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include "Cipher.h"
#include <set>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <assert.h>


#include <openssl/rand.h>
#include <openssl/bn.h>
#include <openssl/ec.h>
#include <openssl/evp.h>
#include <math.h>
#include <iomanip>
#include <vector>

using namespace std;
class Client
{

public:
    const int d;
    const int ser_num;
    const int logR;
    vector<unsigned char *> my_open_str;
    const int n_bytes_agg_bound;
    const BN_ULONG modB;
    const int x0,x1;
    vector<EC_POINT *> encod_vec_LHH;

    // static HomHash hh;
    Cipher cipher;
    // Lagrange lla;
    BIGNUM *b0;
    BIGNUM *b1;

    vector<BN_ULONG> my_vec,A_vec,my_encod_vec0,my_encod_vec1;
    vector<BN_ULONG> PRG_vec0,PRG_vec1;
    vector<BN_ULONG> agg_vec0;
    vector<BN_ULONG> agg_vec1;


public:

    Client(int d,int server_num);
    ~Client();

    int encodegradient(vector<BN_ULONG> &my_vec,vector<BN_ULONG> &A_vec,
                            vector<BN_ULONG> &my_encod_vec1,vector<BN_ULONG> &my_encod_vec2,
                            int x0,int x1);

    int LHHgradient(vector<BN_ULONG> &my_encod_vec0,vector<BN_ULONG> &my_encod_vec1,
                    vector<EC_POINT *> &encod_vec_LHH);
    vector<unsigned char *> PRGgradient(vector<BN_ULONG> &my_encod_vec0,
                            vector<BN_ULONG> &my_encod_vec1,vector<BN_ULONG> &PRG_vec0,
                            vector<BN_ULONG> &PRG_vec1);

    int clverify(unsigned char * agggradient0_buffer, unsigned char * agggradient1_buffer);

    int cltrace(unsigned char * aggencodgradient0_buffer, unsigned char * aggencodgradient1_buffer, EC_POINT *h0, EC_POINT *h1);

    BN_CTX *ctx;

};

void testClient();
vector<bool>  PPRG(vector<bool> src);