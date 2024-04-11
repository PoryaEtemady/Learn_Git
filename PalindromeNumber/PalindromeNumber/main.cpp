#include "qdebug.h"
#include "qlogging.h"
#include <vector>
#include <cmath> // for log10 function

#define SOL1    0//first try
#define SOL2    1//enhanced solution
#define SOL3    2//FASTEST WAY

#define SOLUTION  SOL3

using namespace std;
#if SOLUTION==SOL1
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
#elif SOLUTION==SOL2
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
            long reverse=0;
            while(input_num!=0)
            {
                reverse = reverse*10+input_num%10;
                input_num = input_num/10;
            }

            if(reverse == x)
                return true;
            else
                return false;
        }
    }
};
#elif SOLUTION==SOL3
class Solution
{
public:
    bool isPalindrome(int x)
    {
        string s = to_string(x);
        for (long long i = 0; i < s.length(); i++)
        {
            if (s[i] != s[s.length() - i - 1])
                return false;
        }
        return true;
    }
};
#endif


int main(void)
{
    Solution solution;
    bool val = solution.isPalindrome(123456789);
    qDebug()<<"Value = "<<val;
    return 0;
}
