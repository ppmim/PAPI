/*-----------------------------------------------------------------------------------------
Function SELECT() 
from book: Numerical Recipes in C, p. 341, (by W.H. Press); Second Edition

SELECT() finds the k_th smallest element of an array[1...n] by rearranging the array
in a way that arr[k] holds the k_th smallest element, with all smaller elements moved to
arr[1...k-1] (in arbitrary order) and all larger elements moved to arr[k+1...n] (in 
arbitrary order).
This is the fastest general method for selection, O(n) (rather than O(n*logn) as 
for quicksort). The method rearranges the array and uses partitioning.

float select(unsigned long k, unsigned long n, float arr[])
returns the k_th smallest element out of an array of n elements  

MODIFIED VERSION (RF:21.8.02):
The code was changed to comply with C-array-indices: arr[0], arr[1],...,arr[n-1]
Now, k is really the k_th smallest array element: k=1 being the smallest,..., k=n the 
largest.
n is now the number of array elements, as it is supposed to. 
-----------------------------------------------------------------------------------------*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;		/* Macro definition for function SELECT()*/


float select(unsigned long k, unsigned long n, float arr[]);		/* SELECT() prototype */


float select(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i, ir, j, l, mid;
  float a, temp;

  l=0;			/* modified from: l=1 */
  ir = n-1;		/* modified from: ir=n */

  for(;;)
    {
      if(ir <= l+1)
	{
	  if(ir == l+1 && arr[ir] < arr[l])
	    {
	      SWAP(arr[l],arr[ir])
	    }
	  return arr[k-1];		/* modified from: arr[k] */
	}else{
	  mid=(l+ir) >> 1;
	  SWAP(arr[mid],arr[l+1])
	    if(arr[l] > arr[ir])
	      {
		SWAP(arr[l],arr[ir])
	      }
	  if(arr[l+1] > arr[ir])
	    {
	      SWAP(arr[l+1],arr[ir])
	    }
	  if(arr[l] > arr[l+1])
	    {
	      SWAP(arr[l],arr[l+1])
	    }
	  i=l+1;
	  j=ir;
	  a=arr[l+1];
	  for(;;)
	    {
	      do i++; while(arr[i] < a);
	      do j--; while(arr[j] > a);
	      if(j < i) break;
	      SWAP(arr[i],arr[j])
	    }
	  arr[l+1]=arr[j];
	  arr[j]=a;
	  if(j >= k-1) ir=j-1;		/* modified from: if(j>0k) */
	  if(j <= k-1) l=i;		/* modified from: if(j<=k) */
	}
    } 
}
