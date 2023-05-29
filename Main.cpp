/*
 * Main.cpp
 *
 *  Created on: 2021年6月27日
 *      Author: lin
 */
#include "basic.h"
#include "compressing.h"
#include "utils.h"
#include "split.h"
#define KMER_COUNT 10000000 //定义每个子任务的k_mer个数
int main(int argc, char **argv){

	//开始时间
	struct timeval tvs, tve;
	gettimeofday(&tvs, NULL);
	cout << "start..." << endl;

//	//读取命令行参数，确定执行哪个方法
	if(strcmp(argv[1],"-c")==0){
		compress_main(argc, argv);
	}

    if(strcmp(argv[1],"-p")==0){
        utils(argc, argv);
    }

    if(strcmp(argv[1],"-u")==0){
        sort(argc, argv);
    }

    if(strcmp(argv[1],"-s")==0){
        splitRef("compress_ref.fq",9,32);
    }

	//结束时间
	cout << "end..." << endl;
	gettimeofday(&tve, NULL);
	double span = tve.tv_sec - tvs.tv_sec
			+ (tve.tv_usec - tvs.tv_usec) / 1000000.0;
	cout << "seeding time is: " << span << endl;
}



