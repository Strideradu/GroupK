#include <iostream>
#include <vector>

using namespace std;

int get_suffix_idx(int sa_val, int left, int right, vector<int> &seq_len_array){
    //seq_len_array: an array storing the sequence lengths
    
    int mid = (left+right)/2;
    if(mid==left){
        if(sa_val>=seq_len_array[left]) return right;
        else return left;
    }
    if(sa_val<seq_len_array[mid]){
        right = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
    else{
        left = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
}

int main(){
    int arr[] = {2, 5, 8, 15};
    vector<int> array(arr, arr+4);
    cout<<array.size()<<endl;
    int sa_val = 6;
    int a = get_suffix_idx(sa_val, 0, array.size()-1, array);
    cout<<a<<endl;

    return 0;
}
