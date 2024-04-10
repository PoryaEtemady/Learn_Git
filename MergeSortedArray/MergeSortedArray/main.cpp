#include <stdint.h>
#define  M  3U
#define  N  3U


void merge(int* nums1, int nums1Size, int m, int* nums2, int nums2Size, int n)
{
    int temp=0;
    for(int i = m;i<m+n;i++)
        nums1[i] = nums2[i-m];

    for(int i=1;i<nums1Size;i++)
    {
        for(int j=1;j<nums1Size;j++)
        {
            if(nums1[j-1]>nums1[j])
            {
                temp = nums1[j-1];
                nums1[j-1] = nums1[j];
                nums1[j] = temp;
            }
        }
    }
}


int main(void)
{
    int nums1[M+N] = {4,5,6,0,0,0};
    int nums2[N] = {1,2,3};
    merge(&nums1[0],M+N,M,&nums2[0],N,N);


    return 0;
}
