#pragma once
#define _MAX_BUFFER_SIZE                1024
#define _ULTRA_BUFFER_SIZE              524288                                          // 512KB

#define _ECC_POINT_SIZE                 33                                                      // ECC point compression
#define _SHARE_FIELD_SIZE               (_ECC_POINT_SIZE + 1)
#define _COMMITMENT_SIZE                32                                                      // Folklore hash commitment with SHA256
#define _SYM_KEY_SIZE                   32                                                      // AES-256
#define _IV_SIZE                                16
#define _AGREED_KEY_SIZE                (_SYM_KEY_SIZE + _IV_SIZE)
#define _PID_BYTE_SIZE                  2
#define _SYM_CIPHERTEXT_SIZE    (_SHARE_FIELD_SIZE*4)

//                                                              # of bits
#define _AGG_BOUND_SLACK_SIZE   10
const unsigned char _prg_seed_plaintext[_ULTRA_BUFFER_SIZE] = { 0x00 };