/*
 * compressing.h
 *
 *  Created on: 2021年6月27日
 *      Author: lin
 */

#ifndef SRC_COMPRESSING_H_
#define SRC_COMPRESSING_H_
#include "BplusTreeBit.h"
#include "exactMatchFMindex.h"
#include "Hash.h"
#include "basic.h"

struct para_multi{
    char **p_ref;
    vector<bool> *is_find;
    uint64_t start;
    uint64_t end;
    uint32_t k_code_size;
    uint32_t k_mer;
    uint64_t **read;
    uint32_t ***arrptr;
};


void compress_main(int argc, char **argv);
char *compressRemove(char *mul_ref, uint64_t mul_ref_length,uint64_t *p_remove_len);
void gene_merge(int argc, char **argv);
void compressing(char *mul_ref, uint64_t mul_ref_length, uint32_t k_mer, uint64_t start, uint64_t end,uint32_t thread_num,\
    char *sort_file_path,char *new_sort_file);
void compress_by_array(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,char *new_sort_file);
void toArray(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,uint32_t thread_num,char *sort_file_path);
void *Parall_compress_ref(void *p);
void saveToFile(char *mul_ref, uint64_t mul_ref_length, uint32_t k_mer,
                uint64_t start, uint64_t end);
#endif /* SRC_COMPRESSING_H_ */
