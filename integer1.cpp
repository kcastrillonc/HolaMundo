#include <iostream>

void multiplos3(int nmin, int nmax);
int main(void)
{
  const int m = 1;
  const int n = 100;

  multiplos3(m, n);

  std::cout << "\n";
  return 0;
  
}

void multiplos3(int nmin, int nmax)
{
  for(int ii=nmin;ii<=nmax;ii=ii+1)
  {
    if(ii%3==0){
      std::cout << ii << " ";
    }
   
  }
}