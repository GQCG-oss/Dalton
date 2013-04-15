


template<class T>
int pf(T x)
{
  return x;
}

struct fundat
{
  const char *name;
  int (*f)(int);
};

template<int N>
struct funtem
{
  static fundat dat;
};

template<> fundat funtem<0>::dat = {"knas",pf<int>};
