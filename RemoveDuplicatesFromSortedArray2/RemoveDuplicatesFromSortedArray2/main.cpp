#include <vector>
using namespace std;

class Solution {
public:
    int removeDuplicates(vector<int>& nums)
    {
        int current_element = nums[0];
        int cnt_element =0;
        int end_index = nums.size();
        int current_index =1;
        int cnt_element_num =1;
        for(int cnt=1;current_index<end_index;cnt++)
        {
            if(nums[current_index]==current_element)
            {
                cnt_element++;
                current_index++;
                if(cnt_element>1)
                {
                    for(int i=current_index-1;i<end_index-1;i++)
                    {
                        nums[i] = nums[i+1];
                    }
                    end_index--;
                    current_index--;
                }
            }
            else
            {
                current_element = nums[current_index];
                cnt_element =0;
                current_index++;
            }
        }
        return end_index;
    }
};


int main(void)
{
    Solution solution;
    vector<int> nums;
    int arr[] = {1,1,1,2,2,3,4,4,4,4,5,5,5,5,5,6,7,7,7,7,8};//{0,0,1,1,1,1,2,3,3};
/*                 1,1,2,2,2,3
                 1,1,2,2,2,3
                 1,1,2,2,3,3
*/
    for(int i=0;i<(sizeof(arr)/sizeof(int));i++)
    {
        nums.push_back(arr[i]);
    }
    solution.removeDuplicates(nums);
    return 0;
}
