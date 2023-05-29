//
// Created by lin on 2021/8/3.
//

#include "split.h"
void splitRef(char *p_split_path,uint32_t k_mer,uint8_t thread_num){
    cout << "splitRef start" << endl;

    map<int,uint64_t> my_count;
    cout << "计算$数目" << endl;
    my_count = splitNum(p_split_path);
    cout << "得到$的数目:" << endl;
    vector<string> res;
    vector<string> res_v;
    //数组的数目
    cout << "得到数组的数目" << endl;
    uint64_t count = res_split(res,my_count,p_split_path);
    cout << "数组的数目是：" << count << endl;
    uint64_t *start_kmer;
    uint64_t *end_kmer;
    string s_res = split_merge(res,k_mer,count,&start_kmer,&end_kmer);
    cout << "合并的字符串长度" << s_res.length() << endl;
    ofstream w_file("res_str.fq");
    w_file << s_res;
    w_file.close();
    cout << "splitRef end" << endl;
}

bool cmp_by_value(const pair<uint64_t,uint64_t>& lhs, const pair<uint64_t,uint64_t>& rhs) {
    return lhs.second < rhs.second;
}

string split_merge(vector<string> &res,uint32_t k_mer,uint64_t count,uint64_t **start_kmer,uint64_t **end_kmer){
    cout << "开始合并成一个字符串" << endl;
    string res_str;
    uint64_t k_code_size = K_CODE_SIZE;
    struct bit256KmerPara para_tmp;
    unordered_set<uint64_t> set_index;
    map<uint64_t,uint64_t> map_index;
    vector<string> res_v;

    *start_kmer = (uint64_t *) malloc(sizeof(uint64_t)*count*(k_code_size+1));
    *end_kmer = (uint64_t *) malloc(sizeof(uint64_t)*count*(k_code_size+1));

    uint64_t *start_k_code;
    //定义为3,因为k_mer长度不超过96
    start_k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);

    uint64_t *end_k_code;
    //定义为3,因为k_mer长度不超过96
    end_k_code = (uint64_t*) malloc(sizeof(uint64_t) * k_code_size);

    /**
     * 开始kmer的code和结束kmer的code以及在数组的位置
     */
    for(uint64_t i=0;i < count;i++){
        memset(start_k_code,0,sizeof(uint64_t) * k_code_size);
        memset(end_k_code,0,sizeof(uint64_t) * k_code_size);
        get_para(&para_tmp, k_mer);
        //分別求开始和结束的k_mer的哈希值
        cal_hash_value_directly_256bit_const(res[i].c_str(), start_k_code, para_tmp);
        cal_hash_value_directly_256bit_const(res[i].c_str()+ res[i].length()-k_mer, end_k_code, para_tmp);
        for(uint32_t k=0;k<k_code_size;k++){
            (*start_kmer)[i*(k_code_size+1)+k] = start_k_code[k];
            (*end_kmer)[i*(k_code_size+1)+k] = end_k_code[k];
        }
        (*start_kmer)[i*(k_code_size+1)+k_code_size] = i;
        (*end_kmer)[i*(k_code_size+1)+k_code_size] = i;
    }

    cout << "开始排序" << endl;
    uint64_t * sort_start_kmer = SortFile_umers_for_split(count*(k_code_size+1),*start_kmer,1,k_code_size+1);
    uint64_t * sort_end_kmer = SortFile_umers_for_split(count*(k_code_size+1),*end_kmer,1,k_code_size+1);
    cout << "排序完毕" << endl;

    uint64_t i=0;
    uint64_t j=0;
    while(i<count&& j<count){
        if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==2
            &&sort_end_kmer[i*(k_code_size+1)+k_code_size] != sort_start_kmer[j*(k_code_size+1)+k_code_size]){

//            string start_kmer = res[sort_start_kmer[j*(k_code_size+1)+k_code_size]];
//            string end_kmer = res[sort_end_kmer[i*(k_code_size+1)+k_code_size]];
//            string merge = end_kmer + start_kmer.substr(k_mer);
//            res[sort_end_kmer[i*(k_code_size+1)+k_code_size]] = merge;
            map_index.insert(pair<uint64_t,uint64_t>(sort_end_kmer[i*(k_code_size+1)+k_code_size],sort_start_kmer[j*(k_code_size+1)+k_code_size]));
            set_index.insert(sort_end_kmer[i*(k_code_size+1)+k_code_size]);
            set_index.insert(sort_start_kmer[j*(k_code_size+1)+k_code_size]);
            i++;
            j++;
        } else if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==1){
            j++;
        } else if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+j*(k_code_size+1),k_code_size)==0){
            i++;
        }else{
            if(cmp256BitKmer(sort_end_kmer+i*(k_code_size+1),sort_start_kmer+(j+1)*(k_code_size+1),k_code_size)==2){
                j++;
            }else{
                i++;
            }
        }
    }

    cout << "共有" << count - set_index.size() << "个数组无需拼接" << endl;

//    map<uint64_t,uint64_t>::iterator it_map;
//    for(it_map = map_index.begin();it_map != map_index.end();it_map++){
//        string end_kmer = res[(*it_map).first];
//        string start_kmer = res[(*it_map).second];
//        string merge = end_kmer + start_kmer.substr(k_mer);
//        res[(*it_map).first] = merge;
//    }
    cout << "组合无需拼接的数组" << endl;
    for(uint64_t i = 0;i<count;i++){
        if(!set_index.count(i)){
            res_str += res[i];
            res_str += '$';
        }
    }
    cout << "组合完毕" << endl;

    cout << "开始拼接数组" << endl;
    cout << "需要拼接的数组:" << map_index.size() << endl;

    vector<pair<uint64_t,uint64_t>> key_sort;
    vector<pair<uint64_t,uint64_t>> value_sort;
    map<uint64_t,uint64_t>::iterator it_map_v1;
    map<uint64_t,uint64_t>::iterator it_map_v2;
    set<uint64_t> have_set;
    for(auto temp:map_index){
        key_sort.push_back(temp);
        value_sort.push_back(temp);
    }
    sort(value_sort.begin(),value_sort.end(), cmp_by_value);
    uint64_t test_count = 0;
    for(uint64_t i=0;i<value_sort.size();i++){
        if(!have_set.count(value_sort[i].first)){
            string end_kmer = res[value_sort[i].first];
            string start_kmer = res[value_sort[i].second];
            string merge = end_kmer + start_kmer.substr(k_mer);
            res[value_sort[i].first] = merge;
            have_set.insert(value_sort[i].second);
        }
        for(uint64_t j=0;j<key_sort.size();j++){
            if(value_sort[i].second < key_sort[j].first){
                break;
            }
            if(key_sort[j].first == UINT64_MAX){
                continue;
            }
            if(value_sort[i].second == key_sort[j].first){
                string end_kmer = res[value_sort[i].first];
                string start_kmer = res[value_sort[j].second];
                string merge = end_kmer + start_kmer.substr(k_mer);
                res[value_sort[i].first] = merge;
                have_set.insert(key_sort[j].second);
                value_sort[i].second = key_sort[j].second;
                key_sort[j].first = UINT64_MAX;
                key_sort[j].second = UINT64_MAX;
            }
        }
//        cout << "count =" << ++test_count << endl;
    }

    for(uint64_t i=0;i<value_sort.size();i++){
        if(!have_set.count(value_sort[i].first)){
            res_str += res[value_sort[i].first];
            res_str += '$';
        }
    }

//
//    uint64_t test_count = 0;
//    for(it_map_v1 = map_index.begin();it_map_v1 != map_index.end();it_map_v1++){
//        if((*it_map_v1).second == UINT64_MAX){
//            continue;
//        }
//        string end_kmer = res[(*it_map_v1).first];
//        string start_kmer = res[(*it_map_v1).second];
//        string merge = end_kmer + start_kmer.substr(k_mer);
//        res[(*it_map_v1).first] = merge;
//        for(it_map_v2 = it_map_v1;it_map_v2 != map_index.end();it_map_v2++){
//            if((*it_map_v1).second == (*it_map_v2).first){
//                string end_kmer = res[(*it_map_v1).first];
//                string start_kmer = res[(*it_map_v2).second];
//                string merge = end_kmer + start_kmer.substr(k_mer);
//                (*it_map_v1).second = (*it_map_v2).second;
//                res[(*it_map_v1).first] = merge;
////                map_index.erase(it_map_v2);
//                (*it_map_v2).second = UINT64_MAX;
//
//            }
////            cout << "数组大小" << map_index.size() << endl;
//        }
//        cout << "count =" << ++count << endl;
//    }
    cout << "拼接完毕" << endl;

//    cout << "组合拼接的数组" << endl;
//    for(it_map_v1 = map_index.begin();it_map_v1!=map_index.end();it_map_v1++){
//        if((*it_map_v1).second != UINT64_MAX){
//            res_str += res[(*it_map_v1).first];
//            res_str += '$';
//        }
//    }


    cout << "组合完毕" << endl;

    return res_str;
}

uint64_t res_split(vector<string> &res,map<int,uint64_t> my_count,char *path){
    string temp;
    ifstream readFile(path);
    getline(readFile,temp);
    uint64_t start=0;
    uint64_t end=0;
    uint64_t count = 0;
    for(uint64_t i=0;i<my_count.size()-1;i++){
        if(my_count.at(i)!=my_count.at(i+1)){
            if(temp[my_count.at(i)]=='$'){
                start=my_count.at(i)+1;
            }else{
                start=my_count.at(i);
            }
            if(temp[my_count.at(i+1)]=='$'){
                end=my_count.at(i+1)-1;
            }else{
                end=my_count.at(i+1);
            }
            if(start!=end){
                res.push_back(temp.substr(start,end-start+1).c_str());
                count++;
            }
        }
    }
    cout << "共保存了" << res.size() << "个数组" << endl;
    return count;
}


map<int,uint64_t> splitNum(char *p_ref_path){

    //每个子串压缩后的p_ref
    FILE *read_ref_file = fopen(p_ref_path, "r");
    uint64_t total_size = 0;
    //指针定位到文件尾
    fseek(read_ref_file, 0, 2);
    //计算序列总长度
    total_size = ftell(read_ref_file);
    //指针定位到文件头
    fseek(read_ref_file, 0, 0);

    char *p_ref_compress;
    p_ref_compress = (char*) malloc(sizeof(char) * total_size);
    fread(p_ref_compress, sizeof(char), total_size, read_ref_file);
    fclose(read_ref_file);

    map<int,uint64_t> my_count;

    uint64_t count = 0;
    my_count.insert(pair<int,uint64_t>(0,0));

    for(uint64_t i = 0; i < total_size; i++){
        if (p_ref_compress[i] == '\n' || p_ref_compress[i] == '\0') {
            break;
        }
        if(p_ref_compress[i] == '$'){
            my_count.insert(pair<int,uint64_t>(++count,i));
        }
    }

    my_count.insert(pair<int,uint64_t>(++count,total_size-1));
    cout << "拼接前序列长度是：" << total_size << endl;
    return my_count;
}

//递归方法（不使用中间变量）
int my_strlen(const char *str)
{
    if (*str == '\0')
    {
        return 0;
    }
    else
    {
        return 1 + my_strlen(str + 1);
    }
}