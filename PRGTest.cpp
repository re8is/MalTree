#include "PRGTest.h"
#include <vector>
#include <cstdint>
#include <algorithm>

// 将 64 比特 vector<bool> 扩展为 48 字节 unsigned char 数组
void expand(const std::vector<bool>& input, unsigned char output[]) {
    // 1. 检查输入是否为 64 比特
    if (input.size() != 64) {
        throw std::invalid_argument("Input must be 64 bits.");
    }
    for(int i = 0; i < 48; i++)output[i] = 0x00;

    int w = 1, res = 0;
    for (int i = 0; i < 6; ++i) {
        w = 1, res = 0;
        for(int j = 0; j < 8; j++)
        {
            res += w * input[i*8 + j];
            w *= 2;
        }
        output[i] = res;
    }

    for (int i = 6; i < 8; ++i) {
        w = 1, res = 0;
        for(int j = 0; j < 8; j++)
        {
            res += w * input[i*8 + j];
            w *= 2;
        }
        // start from 32
        output[32+i-6] = res;
    }
    for(int i = 6; i < 6+13; i++)
    {
        output[i] = i;
    }
    for(int i = 34; i < 34+7; i++)
    {
        output[i] = i;
    }

}
vector<bool> PPRG(vector<bool> src)
{
    unsigned char *uc = new unsigned char[48];
    expand(src, uc);
    Cipher cipher(EVP_aes_128_ofb());
    
    // cipher.p = BN_new();
	// BN_hex2bn(&(cipher.p), "EB24862D14F93010BA4267549AE549A239FB9805DC25CA9144E789436B0443A836ED");
    // BIGNUM* b0 = BN_new();
    // BN_rand_range(b0, cipher.p);

    unsigned char *output = new unsigned char[_ULTRA_BUFFER_SIZE];
	unsigned char key0[_SYM_KEY_SIZE];//32
	unsigned char iv0[_IV_SIZE];//16
    
	// BN_ULONG mask_entry0;
    // unsigned char temp0[_AGREED_KEY_SIZE];

    // BN_bn2binpad(b0, temp0, _AGREED_KEY_SIZE);
	memcpy(key0, uc, _SYM_KEY_SIZE);
	memcpy(iv0, uc+_SYM_KEY_SIZE, _IV_SIZE);
    int out_len;
    // cout <<"uc:";
    // for(int i = 0; i < 48; i++)cout << (int)uc[i] << " ";cout <<"\nk0:";
    // for(int i = 0; i < 32; i++)cout << (int)key0[i] << " ";cout <<"\niv0:";
    // for(int i = 0; i < 16; i++)cout << (int)iv0[i] << " ";cout <<"\n";
    cipher.Encrypt(output, out_len, _prg_seed_plaintext, _ULTRA_BUFFER_SIZE, key0, iv0);
    // for(int i = 0; i < _ULTRA_BUFFER_SIZE; i++)
    // {
    //     cout <<std::uppercase << std::hex<< (int)output0[i] << " ";
    // }
    // cout << endl;
    int d = 10;
    // need 258, get 33*8 = 264 bit
    int n_bytes_agg_bound = 33;

    // for(int i=0; i<d; i++)
    // {
    //     mask_entry0 = uchar2ulong(output0 + i*n_bytes_agg_bound, n_bytes_agg_bound);
    //     cout << mask_entry0 << endl;
    // } 
    vector<bool> res(264);
    for(int i = 0; i < n_bytes_agg_bound; i++)
    {
        for(int j = 0; j < 8; j++)
        {
            res[i] = ((output[i*8+j] >> i) & 0x01);
        }
    }
    return res;
    
}