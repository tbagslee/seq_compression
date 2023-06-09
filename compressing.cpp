/*
 * compressing.cpp
 *
 *  Created on: 2021年6月27日
 *      Author: lin
 */
#include "compressing.h"
#include "inputRef.h"
#include "load_DBG_full.h"
#include "Binary_Search.h"
#define K_CODE_SIZE 1 //k_code的大小
#define MAXSIZE 10000000000

void compress_main(int argc, char **argv) {

	cout << "start compressing" << endl;

	char *p_ref_path;

	//k_mer
	uint32_t k_mer;

	uint64_t start;

	uint64_t end;

	uint32_t thread_num;

	//之前是否已经保存了数组
	uint32_t isPre;

	//要拿来生成数组的sort路径
	char *sort_file_path;

    char *new_sort_file;

	for (int32_t i = 2; i < argc; i++) {

		//读入多个基因组 -R必须放在前面
		if (argv[i][0] == '-' && argv[i][1] == 'R') {
            p_ref_path = argv[i + 1];
		}

        if (argv[i][0] == '-' && argv[i][1] == 'm') {
            sort_file_path = argv[i + 1];
        }

        if (argv[i][0] == '-' && argv[i][1] == 'n') {
            new_sort_file = argv[i + 1];
        }

		if (argv[i][0] == '-' && argv[i][1] == 'k') {
			k_mer = atoi(argv[i + 1]); //读入k_mer长度
		}

		if (argv[i][0] == '-' && argv[i][1] == 's') {
			start = atoi(argv[i + 1]);
		}

		if (argv[i][0] == '-' && argv[i][1] == 'e') {
				end = atoi(argv[i + 1]);
		}

        if (argv[i][0] == '-' && argv[i][1] == 't') {
            thread_num = atoi(argv[i + 1]);
        }

        if (argv[i][0] == '-' && argv[i][1] == 'b') {
            isPre = atoi(argv[i + 1]);
        }
	}

	bool isOpen = false;
	if(isPre==1){
	    isPre = true;
	}


	char *p_ref;
	uint64_t ref_length;

	cout << "读入" << p_ref_path << endl;
	ReadSeq(&p_ref, &ref_length, p_ref_path);


	//去除非AGCT字符串后的字符串
	char *p_remove;
	p_remove = (char*) malloc(sizeof(char) * ref_length);
	uint64_t p_remove_length;
	p_remove = compressRemove(p_ref, ref_length, &p_remove_length);

    compressing(p_remove, p_remove_length, k_mer , start , end,thread_num,sort_file_path,new_sort_file);

	cout << "end compressing" << endl;
}

void gene_merge(int argc, char **argv){
    cout << "start compressing" << endl;

    char **p_ref_path; //ref參考基因 用二维指针 因为可能读入多个基因组
    //基因组的个数和基因组名字长度 暂时定为20
    p_ref_path = (char**) malloc(sizeof(char*) * 100);
    *p_ref_path = (char*) malloc(sizeof(char) * 100);
    uint32_t f_count = 0; //读取的基因组的个数

    //k_mer
    uint32_t k_mer;

    uint64_t start;

    uint64_t end;

    uint32_t thread_num;

    //之前是否已经保存了数组
    uint32_t isPre;

    //要拿来生成数组的sort路径
    char *sort_file_path;

    char *new_sort_file;

    for (int32_t i = 2; i < argc; i++) {

        //读入多个基因组 -R必须放在前面
        if (argv[i][0] == '-' && argv[i][1] == 'R') {
            p_ref_path[f_count++] = argv[i + 1];
            for (int32_t j = 2; j <= argc - 3; j++) {
                if (argv[i + j][0] != '-') {
                    p_ref_path[f_count++] = argv[i + j];
                } else {
                    break;
                }
            }
        }

        if (argv[i][0] == '-' && argv[i][1] == 'm') {
            sort_file_path = argv[i + 1];
        }

        if (argv[i][0] == '-' && argv[i][1] == 'n') {
            new_sort_file = argv[i + 1];
        }

        if (argv[i][0] == '-' && argv[i][1] == 'k') {
            k_mer = atoi(argv[i + 1]); //读入k_mer长度
        }

        if (argv[i][0] == '-' && argv[i][1] == 's') {
            start = atoi(argv[i + 1]);
        }

        if (argv[i][0] == '-' && argv[i][1] == 'e') {
            end = atoi(argv[i + 1]);
        }

        if (argv[i][0] == '-' && argv[i][1] == 't') {
            thread_num = atoi(argv[i + 1]);
        }

    }


    //输出多个基因组合并结果 MAXSIZE是合并后的基因的长度 不同服务器它的长度不一样
    uint64_t maxsize = MAXSIZE;
    char *mul_ref;
    mul_ref = (char*) malloc(sizeof(char) * maxsize);
    uint64_t mul_ref_length;

    //每一条基因和其长度
    char *p_ref;
    uint64_t ref_length;
    for (uint32_t i = 0; i < f_count; i++) {
        ReadSeq(&p_ref, &ref_length, p_ref_path[i]);
        //将多个基因组拼接在一起，最后输出合并后的基因和长度
        strcat(mul_ref, p_ref);
    }
    mul_ref_length = strlen(mul_ref);

    //去除非AGCT字符串后的字符串
    char *p_remove;
    p_remove = (char*) malloc(sizeof(char) * mul_ref_length);
    uint64_t p_remove_length;
    p_remove = compressRemove(mul_ref, mul_ref_length, &p_remove_length);

    compressing(p_remove, p_remove_length, k_mer , start , end,thread_num,sort_file_path,new_sort_file);

    cout << "end compressing" << endl;
}

char* compressRemove(char *mul_ref, uint64_t mul_ref_length,
		uint64_t *p_remove_len) {
	cout << "开始去除非AGCT字符" << endl;
	//去除非ACGT的字符
	char *p_remove;
	p_remove = (char*) malloc(sizeof(char) * mul_ref_length);
	uint64_t p_remove_length = 0;
	for (uint64_t i = 0; i < mul_ref_length; i++) {
		if (mul_ref[i] == 'A' || mul_ref[i] == 'C' || mul_ref[i] == 'G'
				|| mul_ref[i] == 'T') {
			p_remove[p_remove_length++] = mul_ref[i];
		} else {
			continue;
		}
	}
	*p_remove_len = p_remove_length;
	cout << "去除完毕后长度:" << p_remove_length << endl;
	return p_remove;
}

void compressing(char *mul_ref, uint64_t mul_ref_length, uint32_t k_mer, uint64_t start, uint64_t end,uint32_t thread_num,\
    char *sort_file_path,char *new_sort_file) {
    cout << "开始压缩" << endl;
//    saveToFile(mul_ref,mul_ref_length,k_mer,start,end);
//    toArray(mul_ref,mul_ref_length,k_mer,start,end,thread_num,sort_file_path);
    compress_by_array(mul_ref,mul_ref_length,k_mer,start,end,new_sort_file);
    cout << "结束压缩" << endl;
}

void compress_by_array(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,\
    char *new_sort_file){
    start = start * 100000000 +1;
    end = end * 100000000;
    cout << "从第:" << start << "个字符开始" << endl;

    uint32_t k_code_size;
    k_code_size = K_CODE_SIZE;

    //打开之前保存的压缩文件
    char *p_ref_compress;
    uint64_t p_ref_compress_length;
    FILE *fin = fopen("compress_ref.fq","r");
    if(!fin){
        cout << "can't find compress_ref.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    p_ref_compress_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    p_ref_compress = (char *)malloc(sizeof(char)*(p_ref_compress_length+100000000));
    fread(p_ref_compress, sizeof(char), p_ref_compress_length, fin);
    fclose(fin);
    cout << "之前文件保存的compress length：" << p_ref_compress_length << endl;

    if(start == 0){
        start = 1;
    }
    //如果是下一个基因组开始的话，旧加上$区分
    if(start == 1 && p_ref_compress[p_ref_compress_length-1]!='$'){
        p_ref_compress[p_ref_compress_length++] = '$';
    }
    //如果结束位置超出范围，重新定义
    if(end > mul_ref_length - k_mer + 1){
        end = mul_ref_length - k_mer + 1;
    }

    vector<bool> is_find;
    bool temp = false;
    uint64_t vector_size;
    //读二进制文件
    ifstream fin_array("array.bat",ios::binary);
    if(!fin_array){
        cout << "can not find array.bat" << endl;
        for(uint64_t i=0;i<mul_ref_length;i++){
            is_find.push_back(temp);
        }
    }else{
        fin_array.read((char *) &vector_size,sizeof(uint64_t));
        for(uint64_t i=0;i<vector_size;i++){
            fin_array.read((char *) &temp, sizeof(bool));
            is_find.push_back(temp);
        }
    }
    fin_array.close();

    uint64_t unique_len = 0;
    cout << "一共读取了" << is_find.size() << "个数组" << endl;
//    for(uint64_t i=0;i<is_find.size();i++){
//        if(is_find[i] == true){
//            unique_len++;
//        }
//    }
//    cout << "该基因组有" << unique_len << "的kmer已经存储" << endl;

    //插入B+树的参数
    struct bit256KmerPara para_tmp;
    //(1) 初始化B+树
    struct NodeBit **p_tmp;
    p_tmp = bit256initialHashFTable();
    //B+树上的ID，不重复
    uint64_t arrayId = 0;
    //定义放入B+树的k_mer长度和code编码
    uint64_t *k_code;
    //先暂时定义3
    k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);
    set<uint32_t> B_index;  //定义set存储k_mer个数

    //查询二进制文件
    char outputFileName[32] = {0};
    sprintf(outputFileName,"%s%s",new_sort_file,"_S");
    uint64_t total_size = 0;
    uint64_t *read;
//    read = (uint64_t*) malloc(
//            sizeof(uint64_t) * p_ref_compress_length*k_code_size);
    FILE *read_binary_file = fopen(outputFileName, "rb");
    if(!read_binary_file){
        cout << "之前没有排序文件" << endl;
        total_size = 0;
    }else{
        cout << "之前有排序文件" << endl;
        //指针定位到文件尾
        fseek(read_binary_file, 0, 2);
        //计算一共几个64位
        total_size = ftell(read_binary_file)/sizeof(uint64_t);
        //指针定位到文件头
        fseek(read_binary_file, 0, 0);

        read = (uint64_t*) malloc(
                sizeof(uint64_t) * total_size);

        fread(read, sizeof(uint64_t), total_size, read_binary_file);
        fclose(read_binary_file);
        cout << "之前二进制文件保存的uint64个数：" << total_size << endl;
    }
    uint64_t u64_size = 0;

    //产生数组
    uint32_t **arrptr;
    if(total_size != 0){
        arrptr = generate_array(total_size / k_code_size, 0);
    }

    cout << "开始压缩" << endl;
    for (uint64_t i = start - 1; i < end; i++){
        if(i%1000000==0){
            cout << i << endl;
        }
        if(!is_find[i]){
            memset(k_code, 0, sizeof(uint64_t) * k_code_size);
            //1)对每一个k_mer进行赋值计算
            get_para(&para_tmp, k_mer);
            cal_hash_value_directly_256bit(mul_ref + i, k_code, para_tmp);
            para_tmp.kmer64Len += 1;
            uint32_t temp = UINT_MAX;
            if(total_size != 0){
                temp = find_arrindexN(arrptr, read, k_code, k_code_size);
            }
            if (temp == UINT_MAX){
                struct nodeBit c_tmp_hashtable;
                c_tmp_hashtable.hashValue = k_code; //定义B+树结点的hash值
                int64_t x;

                //判断B+树上是否存在k_mer
                x = getHashFTableValue(p_tmp, c_tmp_hashtable.hashValue, para_tmp);

                if (x != -1) {
                    //B+树里面有
                    //如果当前k_mer存在与B+树上
                    //判断上一次循环的k_mer是否是已经存在于B+树，最后一个字符是不是$
                    if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                        continue;
                    } else {
                        //添加一个$符号
                        p_ref_compress[p_ref_compress_length++] = '$';
                    }
                } else {
                    //B+树里面没有
                    //如果当前k_mer不在B+树上
                    c_tmp_hashtable.arrayID = arrayId;
                    bit256insertHashFTable(p_tmp, c_tmp_hashtable, para_tmp);
                    arrayId++;
                    uint32_t tmp = bit256hashFFunction(c_tmp_hashtable.hashValue,
                                                       para_tmp);

                    u64_size++;

                    //如果不是第一个k_mer
                    //判断前一个字符是不是$
                    if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                        //前一个字符是$
                        for (uint32_t t = 0; t < k_mer; t++) {
                            p_ref_compress[p_ref_compress_length++] = mul_ref[i + t];
                        }
                    } else {
                        //前一个字符不是$
                        p_ref_compress[p_ref_compress_length++] = mul_ref[i + k_mer - 1];
                    }
                    B_index.insert(tmp);
                }
            }else{
                if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                    continue;
                } else {
                    //添加一个$符号
                    p_ref_compress[p_ref_compress_length++] = '$';
                }
            }
        }else{
            if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                continue;
            } else {
                //添加一个$符号
                p_ref_compress[p_ref_compress_length++] = '$';
            }
        }
    }

    cout << "start write compress_ref.fq" << endl;
    //写入compress_ref
    FILE *write_ref_file = fopen("compress_ref.fq", "w");
    if(!write_ref_file){
        cout << "can't find compress_ref.fq" << endl;
    }
    fwrite(p_ref_compress, sizeof(char), p_ref_compress_length, write_ref_file);
    fclose(write_ref_file);
    free(p_ref_compress);
    cout << "end write" << endl;

    //写入二进制文件
//	cout << "写入二进制文件" << endl;
    cout << "start search B+ tree" << endl;
    cout << "total_size is:" << total_size << endl;
    if(total_size==0){
        read = (uint64_t *) malloc(sizeof(uint64_t)*(total_size+u64_size*k_code_size));
    }else{
        read = (uint64_t *) realloc(read,sizeof(uint64_t)*(total_size+u64_size*k_code_size));
    }

    set<uint32_t>::iterator it;
    for (it = B_index.begin(); it != B_index.end(); it++) {
        if (p_tmp[*it] != NULL) {
            struct NodeBit *r_tmp = Find_Minimal_Node_bit(p_tmp[*it]);
            while (r_tmp != NULL) {
                for (uint32_t k = 0; k < r_tmp->Node_Size; k++) {
//					cout << *(r_tmp->data[k].hashValue) << endl;
                    for (uint32_t hash_index = 0; hash_index < k_code_size;
                         hash_index++) {
                        read[total_size++] =
                                r_tmp->data[k].hashValue[hash_index];
                    }
                }
                r_tmp = r_tmp->brother;
            }
        }
    }

    cout << "end search" << endl;

    cout << "写入二进制文件" << endl;
    FILE *save_binary_file = fopen(new_sort_file, "wb");
    fwrite(read, sizeof(uint64_t), total_size, save_binary_file);
    fclose(save_binary_file);
    free(read);

//    SortFile_umers_forB(new_sort_file, 1, k_code_size);

    cout << "文件保存的字符个数：" << p_ref_compress_length << endl;
    cout << "二进制文件保存的uint64个数：" << total_size << endl;
    cout << "从第:" << end << "个字符结束" << endl;
}

/**
 * 多线程查询排序文件并生成数组
 */
void toArray(char *mul_ref,uint64_t mul_ref_length,uint32_t k_mer,uint64_t start,uint64_t end,uint32_t thread_num,\
    char *sort_file_path){
    start = start * 100000000 +1;
    end = end * 100000000;
    cout << "从第:" << start << "个字符开始" << endl;


    //多线程
    pthread_t *t;
    t = (pthread_t *)malloc(sizeof(pthread_t)*thread_num);

    uint32_t k_code_size;
    k_code_size = K_CODE_SIZE;

    //查询二进制文件
    FILE *read_binary_file = fopen(sort_file_path, "rb");
    if(!sort_file_path){
        cout << "can't find file" << endl;
    }
    uint64_t total_size = 0;
    //指针定位到文件尾
    fseek(read_binary_file, 0, 2);
    //计算一共几个64位
    total_size = ftell(read_binary_file)/sizeof(uint64_t);
    //指针定位到文件头
    fseek(read_binary_file, 0, 0);
    uint64_t *read;
    read = (uint64_t*) malloc(
            sizeof(uint64_t) * total_size);
    fread(read, sizeof(uint64_t), total_size, read_binary_file);
    fclose(read_binary_file);
    cout << "之前二进制文件保存的uint64个数：" << total_size << endl;

    //产生数组
    uint32_t **arrptr;
    arrptr = generate_array(total_size / k_code_size, 0);


    //如果结束位置大于基因组应该结束的位置，就应该重新定义结束位置
    if(start == 0){
        start = 1;
    }
    if(end > mul_ref_length - k_mer +1){
        end = mul_ref_length - k_mer +1;
    }

    //计算每个线程要处理的kmer数目
    uint64_t task_size = end - start + 1;
    vector<bool> is_find;   //计算之前排序文件的是否存在kmer


    //读二进制文件
    ifstream fin_array("array.bat",ios::binary);
    if(!fin_array){
        cout << "can not find array.bat" << endl;
        is_find = vector<bool>(mul_ref_length - k_mer +1, false);
    }else{
        cout << "读取之前的数组文件" << endl;
        bool temp = false;
        uint64_t vector_size;
        fin_array.read((char *) &vector_size,sizeof(uint64_t));
        for(uint64_t i=0;i<vector_size;i++){
            fin_array.read((char *) &temp, sizeof(bool));
            is_find.push_back(temp);
        }
        fin_array.close();
    }



    //计算多线程任务
    uint64_t r = task_size/thread_num;
    uint64_t s = task_size%thread_num;

    //多线程开始位置
    uint64_t cur = start - 1;

    //构造多线程带的参数
    struct para_multi *para;
    para = (struct para_multi*) malloc(
            sizeof(struct para_multi) * thread_num);

    cout << "多线程开始" << endl;
    for(uint32_t i=0;i<thread_num;i++){
        para[i].start = cur;
        if (i < s) {
            cur = cur + r + 1;
        } else {
            cur = cur + r;
        }
        para[i].end = cur-1;
        para[i].is_find = &is_find;
        para[i].p_ref = &mul_ref;
        para[i].k_code_size = K_CODE_SIZE;
        para[i].k_mer = k_mer;
        para[i].arrptr = &arrptr;
        para[i].read = &read;
        if (pthread_create(t + i, NULL, Parall_compress_ref,
                           (void*) (para + i)) != 0) {
            cout << "error!" << endl;
        }
    }
    for (uint32_t i = 0; i < thread_num; i++) {
        pthread_join(t[i], NULL);
    }
    free(para);
    free(t);

    cout << "多线程执行完毕" << endl;

//    uint64_t unique_len = 0;
//    for(uint64_t i=0;i<is_find.size();i++){
//        if(is_find[i] == true){
//            unique_len++;
//        }
//    }
//    cout << "该基因组有" << unique_len << "的kmer已经存储" << endl;

    //写入二进制文件
    ofstream fout("array.bat",ios::binary);
    if(!fout){
        cout << "can not find array.bat" << endl;
    }
    uint64_t array_size = is_find.size();
    bool temp_w;
    cout << "一共存储" << array_size << "个布尔类型" << endl;
    fout.write((char *) &array_size,sizeof(uint64_t));
    for(uint64_t i=0;i<array_size;i++){
        temp_w = is_find[i];
        fout.write((char *) &temp_w, sizeof(bool));
    }
    fout.close();
    free(read);
    cout << "计算完毕" << endl;
    cout << "在第：" << end << "个字符结束" << endl;
}



/**
 * 用于查找排序文件是否存在kmer
 * @param p
 * @return
 */
void *Parall_compress_ref(void *p){
    struct para_multi *tmp = (struct para_multi *)p;
    struct bit256KmerPara para_tmp;
    //定义放入B+树的k_mer长度和code编码
    uint64_t *k_code;
    k_code = (uint64_t*) malloc(sizeof(uint64_t) * tmp->k_code_size);
    for(uint64_t i = tmp->start;i<=tmp->end;i++){
        memset(k_code, 0, sizeof(uint64_t) * tmp->k_code_size);
        //1)对每一个k_mer进行赋值计算
        get_para(&para_tmp, tmp->k_mer);
        cal_hash_value_directly_256bit(*(tmp->p_ref)+i, k_code, para_tmp);
        uint32_t temp;
        temp = find_arrindexN(*(tmp->arrptr), *(tmp->read), k_code, tmp->k_code_size);
        if(temp != UINT_MAX){
            (*(tmp->is_find))[i] = true;
        }
    }
}



void saveToFile(char *mul_ref, uint64_t mul_ref_length, uint32_t k_mer,
                uint64_t start, uint64_t end) {
    cout << "start save first file" << endl;
    start = start * 100000000 +1;
    end = end * 100000000;
    cout << "从第:" << start << "个字符开始" << endl;

    if(start == 0){
        start = 1;
    }
    if(end > mul_ref_length - k_mer +1){
        end = mul_ref_length - k_mer +1;
    }

    uint32_t k_code_size;
    k_code_size = K_CODE_SIZE;

    //	cout << "未压缩前的字符串：" << mul_ref << endl;
    //每个子串压缩后的p_ref
    char *p_ref_compress;
    p_ref_compress = (char*) malloc(sizeof(char) * (end - start)*2);
    uint64_t p_ref_compress_length = 0;

    //插入B+树的参数
    struct bit256KmerPara para_tmp;
    //(1) 初始化B+树
    struct NodeBit **p_tmp;
    p_tmp = bit256initialHashFTable();

    //B+树上的ID，不重复
    uint64_t arrayId = 0;
    //定义放入B+树的k_mer长度和code编码
    uint64_t *k_code;
    //定义为3,因为k_mer长度不超过96
    k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);

    set<uint32_t> B_index;  //定义set存储插入了哪一颗B+树

    //第一个子任务存多少个k_mer
    for (uint64_t i = start - 1; i < end; i++) {
        if(i%1000000==0){
            cout << i << endl;
        }
        //1)对每一个k_mer进行赋值计算
        memset(k_code, 0, sizeof(uint64_t) * k_code_size);
        get_para(&para_tmp, k_mer);
        cal_hash_value_directly_256bit(mul_ref + i, k_code, para_tmp);
//		for (uint32_t i = 0; i < K_CODE_SIZE; i++) {
//			cout << k_code[i] << endl;
//		}
        para_tmp.kmer64Len += 1;

        struct nodeBit c_tmp_hashtable;
        c_tmp_hashtable.hashValue = k_code; //定义B+树结点的hash值
        int64_t x;


        //判断B+树上是否存在k_mer
        x = getHashFTableValue(p_tmp, c_tmp_hashtable.hashValue, para_tmp);

        if (x != -1) {
            //如果当前k_mer存在与B+树上
            //判断上一次循环的k_mer是否是已经存在于B+树，最后一个字符是不是$
            if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                continue;
            } else {
                //添加一个$符号
                p_ref_compress[p_ref_compress_length++] = '$';
            }
        } else {
            //如果当前k_mer不在B+树上
            c_tmp_hashtable.arrayID = arrayId;
            bit256insertHashFTable(p_tmp, c_tmp_hashtable, para_tmp);
            arrayId++;
            uint32_t tmp = bit256hashFFunction(c_tmp_hashtable.hashValue,
                                               para_tmp);
            //set容器为空说明还没插入
            if (B_index.empty()) {
                //如果第一个k_mer
                for (uint32_t t = 0; t < k_mer; t++) {
                    p_ref_compress[p_ref_compress_length++] = mul_ref[i + t];
                }
            } else {
                //如果不是第一个k_mer
                //判断前一个字符是不是$
                if (p_ref_compress[p_ref_compress_length - 1] == '$') {
                    //前一个字符是$
                    for (uint32_t t = 0; t < k_mer; t++) {
                        p_ref_compress[p_ref_compress_length++] = mul_ref[i + t];
                    }
                } else {
                    //前一个字符不是$
                    p_ref_compress[p_ref_compress_length++] = mul_ref[i + k_mer - 1];
                }
            }
            B_index.insert(tmp);
        }
    }

    //存第一个子串
    FILE *save_compress_file = fopen("compress_ref.fq", "w");
    fwrite(p_ref_compress, sizeof(char), p_ref_compress_length, save_compress_file);
    fclose(save_compress_file);
    free(p_ref_compress);
    //存二进制
    uint64_t *kmer_code;
    kmer_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size * (end - start -k_mer)*2);
    uint64_t binary_size = 0;
    set<uint32_t>::iterator it;
    for (it = B_index.begin(); it != B_index.end(); it++) {
        if (p_tmp[*it] != NULL) {
            struct NodeBit *r_tmp = Find_Minimal_Node_bit(p_tmp[*it]);
            while (r_tmp != NULL) {
                for (uint32_t k = 0; k < r_tmp->Node_Size; k++) {
//					cout << *(r_tmp->data[k].hashValue) << endl;
                    for (uint32_t hash_index = 0; hash_index < k_code_size;
                         hash_index++) {
                        kmer_code[binary_size++] =
                                r_tmp->data[k].hashValue[hash_index];
                    }
                }
                r_tmp = r_tmp->brother;
            }
        }
    }

    //写二进制文件
    FILE *save_binary_file = fopen("kmer_binary.bat", "wb");
    fwrite(kmer_code, sizeof(uint64_t), binary_size, save_binary_file);
    fclose(save_binary_file);

    SortFile_umers_forB("kmer_binary.bat", 1, k_code_size);

    cout << "第一个文件保存了" << binary_size / k_code_size << "个k_mer" << endl;
    cout << "压缩后的长度：" << p_ref_compress_length << endl;
    cout << "end save first file" << endl;
}