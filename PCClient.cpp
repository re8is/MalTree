#include"Client.h"

Client::Client(int d,int ser_num)
:cipher(EVP_aes_256_ofb()),d(d),ser_num(ser_num),my_open_str(ser_num),
n_bytes_agg_bound(ceil((logR + _AGG_BOUND_SLACK_SIZE)/8.0)),
modB(((BN_ULONG)1 << (logR + _AGG_BOUND_SLACK_SIZE)) - 1),
encod_vec_LHH(ser_num),my_vec(d),A_vec(d),my_encod_vec0(d),my_encod_vec1(d),
PRG_vec0(d),PRG_vec1(d),x0(2),x1(3),agg_vec0(d),agg_vec1(d),
logR(24)
{

    BN_ULONG mod = ((BN_ULONG)1 << logR) - 1;
    
    ctx = BN_CTX_new();

    for (int i = 0; i < ser_num; ++i)
		my_open_str[i] = new unsigned char[_ECC_POINT_SIZE];
    b0 = BN_new();
    b1 = BN_new();

    BN_rand_range(b0, cipher.p);
    BN_rand_range(b1, cipher.p);

    for (int i=0;i<ser_num;i++)
        encod_vec_LHH[i] = EC_POINT_new(hh.curve);

    for (int i=0;i<d;i++)
    {
        my_vec[i] = rand()&mod;
        A_vec[i] = rand()&mod;
    }

}
Client::~Client()
{
    for (int i = 0; i < ser_num; ++i)
        delete[] my_open_str[i];  // 释放每个指针

    BN_free(b0);
    BN_free(b1);
    BN_CTX_free(ctx);

    for (auto point : encod_vec_LHH)
        EC_POINT_free(point);

    
}

int Client::encodegradient(vector<BN_ULONG> &my_vec,vector<BN_ULONG> &A_vec,
                           vector<BN_ULONG> &my_encod_vec0,vector<BN_ULONG> &my_encod_vec1,
                           int x0,int x1)
{
    for (int i=0;i<my_vec.size();i++)
    {
        my_encod_vec0[i] = my_vec[i] + A_vec[i]*x0;
        my_encod_vec1[i] = my_vec[i] + A_vec[i]*x1;
    }

/* #ifdef DEBUG
    cout<<my_vec[0]<<"\n";
    cout<<A_vec[0]<<"\n";
    cout<<my_encod_vec0[0]<<"\n";
#endif */

    return 0;
}

HomHash Client::hh(100000);

int Client::LHHgradient(vector<BN_ULONG> &my_encod_vec0,
                        vector<BN_ULONG> &my_encod_vec1,vector<EC_POINT *> &encod_vec_LHH)
{
    hh.Hash(encod_vec_LHH[0],my_encod_vec0);
    hh.Hash(encod_vec_LHH[1],my_encod_vec1);
    
/* #ifdef DEBUG
    EC_POINT *hash = EC_POINT_new(hh.curve);
    EC_POINT *agg_hash = EC_POINT_new(hh.curve);

    vector<BN_ULONG> agg_vec(this->d);
    for(int i=0;i<this->d;i++)
    {
        agg_vec[i] = my_encod_vec0[i] + my_encod_vec1[i];
    }
    hh.Hash(hash,agg_vec);
    hh.Eval(agg_hash,encod_vec_LHH);

    assert( EC_POINT_cmp(hh.curve, agg_hash, hash, ctx) == 0 );

#endif */

    vector<BIGNUM *>hash_bignum(ser_num);
    for(int i=0;i<ser_num;i++)
    {
        hash_bignum[i] = BN_new();
        EC_POINT_point2bn(hh.curve, encod_vec_LHH[i],
		POINT_CONVERSION_COMPRESSED, hash_bignum[i], ctx);

        BN_bn2binpad(hash_bignum[i], my_open_str[i], _ECC_POINT_SIZE);
    }
    //需要发送my_open_str，总共2*_ECC_POINT_SIZE+转发服

/* #ifdef DEBUG
    char* bn_hex = BN_bn2hex(hash_bignum[1]);
    cout << "BIGNUM[" << 1 << "] value: " << bn_hex << endl;
    OPENSSL_free(bn_hex);

    cout << "Binary data: ";
    for (int j = 0; j < _ECC_POINT_SIZE; ++j) {
        cout << hex << setw(2) << setfill('0') 
             << (int)my_open_str[1][j] << " ";
    }
    cout << dec << endl;

    vector<EC_POINT *> bignum_hash(ser_num);
    vector<BIGNUM *> str_bignum(ser_num);
    for(int i=0;i<ser_num;i++)
    {
        str_bignum[i] = BN_new();
        bignum_hash[i] = EC_POINT_new(hh.curve);
        BN_bin2bn(my_open_str[i],_ECC_POINT_SIZE,str_bignum[i]);
        EC_POINT_bn2point(hh.curve, str_bignum[i], bignum_hash[i], ctx);
        assert( EC_POINT_cmp(hh.curve, bignum_hash[i], encod_vec_LHH[i], ctx) == 0 );
    }
#endif */
    for (auto bn : hash_bignum)
        BN_free(bn);

    return 0;

}

vector<unsigned char *> Client::PRGgradient(vector<BN_ULONG> &my_encod_vec0,
                        vector<BN_ULONG> &my_encod_vec1,vector<BN_ULONG> &PRG_vec0,
                        vector<BN_ULONG> &PRG_vec1)
{
    vector<unsigned char *> buffer(ser_num);
    for (auto& ptr : buffer) 
    {
        ptr = new unsigned char[d * n_bytes_agg_bound];  // 按需分配
    }

	int out_len;
	unsigned char *output0 = new unsigned char[_ULTRA_BUFFER_SIZE];
	unsigned char key0[_SYM_KEY_SIZE];//32
	unsigned char iv0[_IV_SIZE];//16
	BN_ULONG mask_entry0;
    unsigned char temp0[_AGREED_KEY_SIZE];
    unsigned char *output1 = new unsigned char[_ULTRA_BUFFER_SIZE];
	unsigned char key1[_SYM_KEY_SIZE];//32
	unsigned char iv1[_IV_SIZE];//16
	BN_ULONG mask_entry1;
    unsigned char temp1[_AGREED_KEY_SIZE];


    BN_bn2binpad(b0, temp0, _AGREED_KEY_SIZE);
	memcpy(key0, temp0, _SYM_KEY_SIZE);
	memcpy(iv0, temp0, _IV_SIZE);
    BN_bn2binpad(b1, temp1, _AGREED_KEY_SIZE);
	memcpy(key1, temp1, _SYM_KEY_SIZE);
	memcpy(iv1, temp1, _IV_SIZE);

    cipher.Encrypt(output0, out_len, _prg_seed_plaintext, _ULTRA_BUFFER_SIZE, key0, iv0);
    cipher.Encrypt(output1, out_len, _prg_seed_plaintext, _ULTRA_BUFFER_SIZE, key1, iv1);

    for(int i=0; i<d; i++)
    {
        mask_entry0 = uchar2ulong(output0 + i*n_bytes_agg_bound, n_bytes_agg_bound);
        PRG_vec0[i] = (my_encod_vec0[i] + (mask_entry0)&modB)&modB;
        /* cout<<my_encod_vec0[0]<<"\n";
        cout<<((mask_entry0)&modB);
        cout<<"\n";
        cout<<PRG_vec0[0]<<"\n"; */

        mask_entry1 = uchar2ulong(output1 + i*n_bytes_agg_bound, n_bytes_agg_bound);
        PRG_vec1[i] = (my_encod_vec1[i] + (mask_entry1)&modB)&modB;
    }

    int buffer_offset = 0;
    for (int i = 0; i < d; ++i)
    {
        ulong2uchar(buffer[0] + buffer_offset, PRG_vec0[i], n_bytes_agg_bound);
        buffer_offset += n_bytes_agg_bound;
    }
    buffer_offset = 0;
    for (int i = 0; i < d; ++i)
    {
        ulong2uchar(buffer[1] + buffer_offset, PRG_vec1[i], n_bytes_agg_bound);
        buffer_offset += n_bytes_agg_bound;
    } 
/* #ifdef DEBUG
    cout << "buffer[0] content (hex): ";
    for (int i = 0; i < d * n_bytes_agg_bound; ++i) 
    {
        cout << hex << setw(2) << setfill('0') 
        << static_cast<int>(buffer[0][i]) << " ";
    }
    cout << dec << endl;


    vector<BN_ULONG>encod_vec(d);
    int offset = 0;
    for(int i=0;i<d;i++)
    {
        encod_vec[i] = uchar2ulong(buffer[0]+ offset, n_bytes_agg_bound);
        offset += n_bytes_agg_bound;
    }

    cout<<encod_vec[3]<<"\n";
    cout<<PRG_vec0[3]<<"\n";


    vector<BN_ULONG>vec0(d);
    BN_bn2binpad(b0, temp0, _AGREED_KEY_SIZE);
	memcpy(key0, temp0, _SYM_KEY_SIZE);
	memcpy(iv0, temp0, _IV_SIZE);
    cipher.Encrypt(output0, out_len, _prg_seed_plaintext, _ULTRA_BUFFER_SIZE, key0, iv0);
    for(int i=0;i<d;i++)
    {
        mask_entry0 = uchar2ulong(output0 + i*n_bytes_agg_bound, n_bytes_agg_bound);
        vec0[i] = (PRG_vec0[i] - (mask_entry0)&modB)&modB;
    }
    for(int i=0;i<d;i++)
    {
        cout<<vec0[i]<<" ....."<<my_encod_vec0[i]<<"\n";
    }
    
    cout<<"------------------------------------------------------------------------------------------------"<<"\n";
#endif

#ifdef DEBUG
    cout << "buffer[0] content (hex): ";
    for (int i = 0; i < 1 * n_bytes_agg_bound; ++i) 
    {
        cout << hex << setw(2) << setfill('0') 
        << static_cast<int>(buffer[0][i]) << " ";
    }
    cout << dec << endl;
#endif */

    return buffer;

    for (auto ptr : buffer)
        delete[] ptr;
    
    delete[] output0;
    delete[] output1;
    

}

int Client::clverify(unsigned char * agggradient0_buffer,unsigned char * agggradient1_buffer)
{

    int offset = 0;
    for(int i=0;i<d;i++)
    {
        agg_vec0[i] = uchar2ulong(agggradient0_buffer + offset, n_bytes_agg_bound);
        offset += n_bytes_agg_bound;
    }
    offset = 0;
    for(int i=0;i<d;i++)
    {
        agg_vec1[i] = uchar2ulong(agggradient1_buffer + offset, n_bytes_agg_bound);
        offset += n_bytes_agg_bound;
    }

    if(agg_vec0 == agg_vec1)
        return 0;
    else
        return 1;
}

int Client::cltrace(unsigned char * aggencodgradient0_buffer, unsigned char * aggencodgradient1_buffer, 
                    EC_POINT *h0, EC_POINT *h1)
{
    vector<BN_ULONG> agg_encodvec0(d);
    vector<BN_ULONG> agg_encodvec1(d);

    int offset = 0;
    for(int i=0;i<d;i++)
    {
        agg_encodvec0[i] = uchar2ulong(aggencodgradient0_buffer + offset, n_bytes_agg_bound);
        offset += n_bytes_agg_bound;
    }
    offset = 0;
    for(int i=0;i<d;i++)
    {
        agg_encodvec1[i] = uchar2ulong(aggencodgradient1_buffer + offset, n_bytes_agg_bound);
        offset += n_bytes_agg_bound;
    }

    EC_POINT *h_0 = EC_POINT_new(hh.curve);
    EC_POINT *h_1 = EC_POINT_new(hh.curve);
    hh.Hash(h_0,agg_encodvec0);
    hh.Hash(h_1,agg_encodvec1);
    if(EC_POINT_cmp(hh.curve, h_0, h0, ctx)==1)
    {
        EC_POINT_free(h_0);
        return 0;
    }
    else if(EC_POINT_cmp(hh.curve, h_1, h1, ctx)==1)
    {
        EC_POINT_free(h_1);
        return 1;
    }
    else
    {
        vector<BN_ULONG> agg_vec(d);
        vector<BIGNUM*> MF0 = lla.Combine(x0,x1);
        BIGNUM *secret1_0 = BN_new();
        BIGNUM *secret2_0 = BN_new();
        BIGNUM *add1_0 = BN_new();
        BIGNUM *add2_0 = BN_new();
        BIGNUM *add_0 = BN_new();

        for(int i=0; i<d; i++)
        {
            BN_set_word(secret1_0,agg_encodvec0[i]);
            BN_set_word(secret2_0,agg_encodvec1[i]);
            BN_mod_mul(add1_0, MF0[0], secret1_0, lla.p, ctx);
            BN_mod_mul(add2_0, MF0[1], secret2_0, lla.p, ctx);
            BN_mod_add(add_0, add2_0, add1_0, cipher.p, lla.ctx);
            agg_vec[i] = BN_get_word(add_0);
        }

        BN_free(secret1_0);
        BN_free(secret2_0);
        BN_free(add1_0);
        BN_free(add2_0);
        BN_free(add_0);
        for (auto bn : MF0) BN_free(bn);
        if(agg_vec != agg_vec0)
            return 0;
        else if(agg_vec != agg_vec1)
            return 1;
        
    }
}



 void testClient(){
// int main(){

    int d = 50;
    int ser_num = 2;

    Client cl(d,ser_num);

    cl.encodegradient(cl.my_vec,cl.A_vec,cl.my_encod_vec0,cl.my_encod_vec1,cl.x0,cl.x1);
    
/*     cl.LHHgradient(cl.my_encod_vec0,cl.my_encod_vec1,cl.encod_vec_LHH);

    cl.PRGgradient(cl.my_encod_vec0,cl.my_encod_vec1,cl.PRG_vec0,cl.PRG_vec1); */

    unsigned char *buffer0 = new unsigned char[d * cl.n_bytes_agg_bound];
    unsigned char *buffer1 = new unsigned char[d * cl.n_bytes_agg_bound];
    int buffer_offset = 0;
    for (int i = 0; i < d; ++i)
    {
        ulong2uchar(buffer0 + buffer_offset, cl.my_encod_vec0[i], cl.n_bytes_agg_bound);
        buffer_offset += cl.n_bytes_agg_bound;
    }
    buffer_offset = 0;
    for (int i = 0; i < d; ++i)
    {
        ulong2uchar(buffer1 + buffer_offset, cl.my_encod_vec1[i], cl.n_bytes_agg_bound);
        buffer_offset += cl.n_bytes_agg_bound;
    }

    int y = cl.clverify(buffer0,buffer1);
    cout<<y<<"\n";

}

//g++ -DDEBUG Client.cpp HomHash.cpp -o Client -std=c++11 -lssl -lcrypto
//g++ Client.cpp HomHash.cpp -o Client -std=c++11 -I/gpfs/home/kezhijie/anaconda3/include -L/gpfs/home/kezhijie/anaconda3/lib -lssl -lcrypto
//export LD_LIBRARY_PATH=/gpfs/home/kezhijie/anaconda3/lib:$LD_LIBRARY_PATH


//需要发送my_open_str，总共2*_ECC_POINT_SIZE + 转发服服务器的2*2*_ECC_POINT_SIZE（与梯度维度无关）4*33
//需要发送buffer 2*d*n_bytes_agg_bound（与梯度维度相关）