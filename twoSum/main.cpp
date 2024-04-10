#include <iostream>
#include <unordered_map>
#include <vector>
#define ALG  0//0 ==FAST | 1==MEMORY
using namespace std;

#if ALG
class Solution
{
public:
    vector<int> twoSum(vector<int>& nums, int target)
    {
        vector<int> output_vec;
        bool find_flag = false;
        for(int i=0;(i<nums.size() && !find_flag);i++)
        {
            int remain = target - nums.at(i);
            for(int j=0;(j<nums.size() && !find_flag);j++)
            {
                if(i != j)
                {
                    if(remain - nums.at(j) == 0)
                    {
                        output_vec.push_back(i);
                        output_vec.push_back(j);
                        find_flag = true;
                    }
                }
            }
        }
        return output_vec;
    }
};
#else

class Solution {
public:
    vector<int> twoSum(vector<int>& nums, int target)
    {
        unordered_map<int, int> num_indices;
        vector<int> output_vec;

        for (int i = 0; i < nums.size(); ++i)
        {
            int remain = target - nums[i];
            // Check if the complement (target - current element) exists in the hash table
            if (num_indices.count(remain))
            {
                output_vec.push_back(num_indices[remain]);
                output_vec.push_back(i);
                return output_vec; // Early exit once a pair is found
            }
            // Add the current element and its index to the hash table
            num_indices[nums[i]] = i;
        }
        // No pair found
        return output_vec;
    }
};

#endif

int main(int argc, char *argv[])
{

    Solution solution;
    vector<int> input_vec= {2,7,11,14};
    vector<int> output_vect= solution.twoSum(input_vec,9);
    for(int i=0;i<output_vect.size();i++)
        cout<<output_vect[i]<<endl;

    return 0;
}


