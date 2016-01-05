#ifndef TOOLS_H_
#define TOOLS_H_

#include <cmath>
#include <ctime>

// Some functions in this header are inlined to avoid "multiple declaration"
// linker errors resulting from including this file in multiple C files.

// draws a normally distributed variable using the Marsaglia polar method
template<typename T>
static T randn(T mu = 0, T sigma = 1)
{
  static bool have_deviate = false;
  static T stored_deviate;
  if (have_deviate)
  {
    have_deviate = false;
    return stored_deviate * sigma + mu;
  }
  else
  {
    T v1, v2, s;
    do
    {
      v1 = 2 * (std::rand() / (T)RAND_MAX) - 1;
      v2 = 2 * (std::rand() / (T)RAND_MAX) - 1;
      s = v1 * v1 + v2 * v2;
    } while (s >= 1 || s == 0);
    const T polar = std::sqrt(-2 * std::log(s) / s);
    stored_deviate = v1 * polar;
    have_deviate = true;
    return v2 * polar * sigma + mu;
  }
}

// returns the number of seconds elapsed since start
template<typename T>
inline T time_elapsed(const timespec& start)
{
  timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return static_cast<T>(now.tv_sec - start.tv_sec + 1e-9 * (now.tv_nsec - start.tv_nsec));
}

// prints the current datetime in a nice format to a char array
inline void timeNow(char * now)
{
  time_t rawtime;
  std::time(&rawtime);
  struct tm * timeinfo = localtime(&rawtime);
  strftime(now, 20, "%Y-%m-%d %H:%M:%S", timeinfo);
}

// a branchless, type-safe, accurate, fast implementation of the signum function
template <typename T>
inline int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

#endif /* TOOLS_H_ */
