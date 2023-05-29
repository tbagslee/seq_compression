//
// Created by lin on 2021/7/27.
//

#ifndef FINAL_COMPRESS_UTILS_H
#define FINAL_COMPRESS_UTILS_H
#include "basic.h"
#include "inputRef.h"
#include "load_DBG_full.h"
void utils(int argc, char **argv);
void copy_utils();
void copy();
void sort(int argc, char **argv);
void thread_read_file(int tid, const string& file_path);
void test_detach(const string& file_path);
#endif //FINAL_COMPRESS_UTILS_H
