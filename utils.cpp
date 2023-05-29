//cd ..
// Created by lin on 2021/7/27.
//
#include "utils.h"
#define MAXSIZE 10000000000
#define COPY_SIZE 2005647323
#define K_CODE_SIZE 1

void utils(int argc, char **argv){

    cout << "start merge" << endl;

    char **p_ref_path; //ref參考基因 用二维指针 因为可能读入多个基因组
    //基因组的个数和基因组名字长度 暂时定为20
    p_ref_path = (char**) malloc(sizeof(char*) * 100);
    *p_ref_path = (char*) malloc(sizeof(char) * 100);
    uint32_t f_count = 0; //读取的基因组的个数

    for (int32_t i = 2; i < argc; i++){
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

    //写入compress_ref
    FILE *write_ref_file = fopen("compress_ref.fq", "w");
    if(!write_ref_file){
        cout << "can't find compress_ref.fq" << endl;
    }
    fwrite(mul_ref, sizeof(char), mul_ref_length, write_ref_file);
    fclose(write_ref_file);
    cout << "write length:" << mul_ref_length << endl;

    cout << "check out write" << endl;
    char *check_ref;
    uint64_t check_length;
    FILE *fin = fopen("compress_ref.fq","r");
    if(!fin){
        cout << "can't find compress_ref.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    check_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    check_ref = (char *)malloc(sizeof(char)*check_length);
    fread(check_ref, sizeof(char), check_length, fin);
    fclose(fin);
    cout << "read length:" << check_length << endl;

    cout << "end merge" << endl;
}

void copy(){
    cout << "utils for copy correct char" << endl;
    char *p_ref;
    uint_fast64_t ref_length;
    char *p_ref_path = "compress_ref.fq";
    ReadSeq(&p_ref, &ref_length, p_ref_path);
    uint64_t copy_size = COPY_SIZE;
    char *copy;
    copy = (char *)malloc(sizeof(char)*ref_length);
    for(uint64_t i=0;i<copy_size;i++){
        copy[i]=p_ref[i];
    }
    uint64_t copy_size2;
    copy_size2 = copy_size + 502145707;
    for(uint64_t i=copy_size2-1;i<ref_length;i++){
        copy[copy_size++] = p_ref[i];
    }

    char *p_ref_compress;
    uint64_t p_ref_compress_length;
    //写入compress_ref
    FILE *write_ref_file = fopen("compress_ref_v3.fq", "w");
    if(!write_ref_file){
        cout << "can't find compress_ref.fq" << endl;
    }
    fwrite(copy, sizeof(char), copy_size, write_ref_file);
    fclose(write_ref_file);
    cout << "end copy" << endl;

    cout << "test my copy" << endl;
    FILE *fin = fopen("compress_ref_v3.fq","r");
    if(!fin){
        cout << "can't find compress_ref_v3.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    p_ref_compress_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    p_ref_compress = (char *)malloc(sizeof(char)*p_ref_compress_length);
    fread(p_ref_compress, sizeof(char), p_ref_compress_length, fin);
    fclose(fin);
    cout << "文件保存的compress length：" << p_ref_compress_length << endl;
    cout << "last char:" << p_ref_compress[p_ref_compress_length-1] << endl;
}

void copy_utils(){
    cout << "utils for copy correct char" << endl;
    char *p_ref;
    uint_fast64_t ref_length;
    char *p_ref_path = "compress_ref_v2.fq";
    ReadSeq(&p_ref, &ref_length, p_ref_path);
    uint64_t copy_size = COPY_SIZE;
    char *copy;
    copy = (char *)malloc(sizeof(char)*copy_size);
    for(uint64_t i=0;i<copy_size;i++){
        copy[i]=p_ref[i];
    }

    //打开之前保存的压缩文件
    char *p_ref_compress;
    uint64_t p_ref_compress_length;
    //写入compress_ref
    FILE *write_ref_file = fopen("compress_ref_v2.fq", "w");
    if(!write_ref_file){
        cout << "can't find compress_ref.fq" << endl;
    }
    fwrite(copy, sizeof(char), copy_size, write_ref_file);
    fclose(write_ref_file);
    cout << "end copy" << endl;

    cout << "test my copy" << endl;
    FILE *fin = fopen("compress_ref_v2.fq","r");
    if(!fin){
        cout << "can't find compress_ref_v2.fq" << endl;
    }
    fseek(fin, 0, 2);
    //计算序列总长度
    p_ref_compress_length = ftell(fin);
    //指针定位到文件头
    fseek(fin, 0, 0);
    p_ref_compress = (char *)malloc(sizeof(char)*p_ref_compress_length);
    fread(p_ref_compress, sizeof(char), p_ref_compress_length, fin);
    fclose(fin);
    cout << "文件保存的compress length：" << p_ref_compress_length << endl;
    cout << "last char:" << p_ref_compress[p_ref_compress_length-1] << endl;
}



void sort(int argc, char **argv){
    //要拿来生成数组的sort路径
    char *sort_file_path;

    for (int32_t i = 2; i < argc; i++){
        if (argv[i][0] == '-' && argv[i][1] == 'm') {
            sort_file_path = argv[i + 1];
        }
    }
    SortFile_umers_forB(sort_file_path, 32, K_CODE_SIZE);
}


//读文件的函数
void thread_read_file(int tid, const string& file_path)
{
    //tid是指线程序号
    ifstream file(file_path.c_str(), ios::in);	//ios::in表示以输入的方式打开一个			文件，并从里面读取数据
    if (!file.good()) {		//file.good() 是判断文件是否正常，这是文件流的一个函数，当文件出错时进行此循环
        stringstream ss;	//<sstream>
        ss << "Thread " << tid << " failed to open file: " << file_path << "\n";
        cout << ss.str();
        return;
    }

    int pos;
    if (tid == 0) pos = 0;
    else pos = tid * 10;

    file.seekg(pos, ios::beg);
    /*
    * seekg 对输入文件定位（seek—寻找，g—get），有两个参数：
    * 第一个：表示偏移量，可正可负，正表示向后，负表示向前
    * 第二个：偏移的基地址，可以是：
    * ios::beg 输入流的开始
    * ios::cur 输入流的当前位置
    * ios::end 输入流的结束
    */

    string line;
    getline(file, line);	//读文件file的内容存入line中
    stringstream ss;
    ss << "Thread " << tid << ", pos=" << pos << ": " << line << endl;
    cout << ss.str();
}

void test_detach(const string& file_path)
{
    for (int i = 0; i < 10; ++i) {

        //std::thread  th = thread (thread_read_file, i, file_path);—和下面一样的
        std::thread  th(thread_read_file, i, file_path);//生成多线程，std::thread构造函数 在 #include<thread> 头文件中声明
        //i和file_path是thread_read_file的参数，i是指线程号
        th.detach();
    }
}

//void test_join(const string& file_path)
//{
//	vector<std::thread> vec_threads;	//线程型数组vec_threads
//	for (int i = 0; i < 10; ++i) {
//		std::thread  th(thread_read_file, i, file_path);
//
//		//添加一个新元素到结束的容器。该元件是构成在就地，即没有复制或移动操作进行。
//		//std::move函数可以以非常简单的方式将左值引用转换为右值引用
//		vec_threads.emplace_back(std::move(th));  // push_back() is also OK
//	}
//
//	auto it = vec_threads.begin();
//	for (; it != vec_threads.end(); ++it) {
//		(*it).join();	//相当于将动态数组里的变量（线程）全设置成join
//	}
//}


//int main()
//{
//    string file_path = "./1.txt";
//
//    test_detach(file_path);
//
//    //表示当前线程休眠一段时间，休眠期间不与其他线程竞争CPU，根据线程需求，等待若干时间。
//    std::this_thread::sleep_for(std::chrono::seconds(1));  // wait for detached threads done，std::chrono处理时间的类
//
//    //test_join(file_path);
//    //test_detach(file_path);
//    return 0;
//}



