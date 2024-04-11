#include "qdebug.h"
#include "qlogging.h"
#include <vector>
#include <cmath> // for log10 function
using namespace std;
class Solution {
public:
    bool isPalindrome(int x)
    {
        int input_num = x;
        if(input_num<0)
            return false;
        else
        {
            if(input_num==0)
                return true;
            int num_digits = floor(log10(input_num)) + 1; // floor rounds down to nearest integer
            vector<int> digits;
            while(input_num!=0)
            {
                digits.push_back(input_num%10);
                input_num = input_num/10;
            }

            for (int cnt = 0; cnt < (num_digits/2); cnt++)
            {
                if(digits[cnt]!=digits[num_digits-1-cnt])
                    return false;
            }
            return true;
        }
    }
};



int main(void)
{
    Solution solution;
    bool val = solution.isPalindrome(1234321);
    qDebug()<<"Value = "<<val;
    return 0;
}
