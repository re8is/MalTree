#pragma once

#include <openssl/bn.h>

#include "Config.h"

int encodePid(unsigned char *buffer, int pid);

int decodePid(int &pid, unsigned char *buffer);

void ulong2uchar(unsigned char *buffer, BN_ULONG val, int n_bytes);

BN_ULONG uchar2ulong(unsigned char *buffer, int n_bytes);

int BN_bn2binpad(const BIGNUM *a, unsigned char *to, int tolen);

