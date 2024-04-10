#include <vector>
using namespace std;

class Solution {
public:
    int removeDuplicates(vector<int>& nums)
    {
        vector<int> new_vec;
        int head =0;
        new_vec.push_back(nums[0]);
        for(int cnt=1;cnt<nums.size();cnt++)
        {
            if(nums[cnt]!=new_vec[head])
            {
                new_vec.push_back(nums[cnt]);
                head++;
            }

        }
        for(int cnt=0;cnt<nums.size();cnt++)
        {
            if(cnt<new_vec.size())
                nums[cnt] = new_vec[cnt];
            else
                nums[cnt] =0;
        }
        return new_vec.size();
    }
};


int main(void)
{
    Solution solution;
    vector<int> nums;
    int arr[] = {0,0,1,1,1,2,2,3,3,4};
    for(int i=0;i<(sizeof(arr)/sizeof(int));i++)
    {
        nums.push_back(arr[i]);
    }
    solution.removeDuplicates(nums);
    return 0;
}
