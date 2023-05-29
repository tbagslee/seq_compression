//
// Created by lin on 2021/8/3.
//

#ifndef FINAL_COMPRESS_SPLIT_H
#define FINAL_COMPRESS_SPLIT_H
#include "basic.h"
#include "load_DBG_full.h"
#define K_CODE_SIZE 1
uint64_t res_split(vector<string> &res,map<int,uint64_t> my_count,char *path);
void splitRef(char *p_split_path,uint32_t k_mer,uint8_t thread_num);
map<int,uint64_t> splitNum(char *p_ref_path);
string split_merge(vector<string> &res,uint32_t k_mer,uint64_t count,uint64_t **start_kmer,uint64_t **end_kmer);
int my_strlen(const char *str);
#endif //FINAL_COMPRESS_SPLIT_H
