#include "Utils.h"
#include <string.h>

// NOTE: both encodePid and decodePid assume pid is an non-negative
//   integer with less than 16 bits
int encodePid(unsigned char *buffer, int pid)
{
	buffer[0] = pid&0xff;
	buffer[1] = (pid >> 8)&0xff;
	return _PID_BYTE_SIZE;
}

int decodePid(int &pid, unsigned char *buffer)
{
	pid = ((buffer[1]&0xff) << 8) + (buffer[0]&0xff);
	return _PID_BYTE_SIZE;
}

void ulong2uchar(unsigned char *buffer, BN_ULONG val, int n_bytes)
{
	for (BN_ULONG i = 0; i < n_bytes; ++i)
		buffer[i] = (val >> (8*i))&0xff;
}



BN_ULONG uchar2ulong(unsigned char *buffer, int n_bytes)
{
	BN_ULONG ret = 0;
	for (BN_ULONG i = 0; i < n_bytes; ++i)
		ret += ((BN_ULONG)(buffer[i]&0xff) << (8*i));
	return ret;
}
int BN_bn2binpad(const BIGNUM *a, unsigned char *to, int tolen) {
    if (tolen <= 0) return 0;

    // 获取实际需要的字节数
    int num_bytes = BN_num_bytes(a);
    
    // 计算需要填充的零字节数
    int pad_len = tolen - num_bytes;

    // 处理截断情况（tolen不足时）
    if (pad_len < 0) {
        // 仅复制高位tolen字节（类似BN_bn2binpad的截断行为）
        memset(to, 0, tolen);
        BN_bn2bin(a, to + (-pad_len));
        return tolen;
    }

    // 处理填充情况（tolen足够时）
    memset(to, 0, pad_len);                      // 高位填充零
    BN_bn2bin(a, to + pad_len);                  // 写入有效数据
    return tolen;
}