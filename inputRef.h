/*
 * inputRef.h
 *
 *  Created on: 2021年5月19日
 *      Author: lin
 */

#ifndef INPUTREF_H_
#define INPUTREF_H_
#include "basic.h"
void ReadSeq(char **seq1,uint64_t *seq_length,char* p_ref);
void ReadSeq_v2(char **seq1, uint64_t *seq_length, char *p_ref,uint64_t size);



#endif /* INPUTREF_H_ */
