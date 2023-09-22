#include <iostream>
#include <numeric>
#include <array>
#include <chrono>
#include <immintrin.h>
//-- FOR AVX
constexpr auto AVX_FLOAT_COUNT = 8u;
//-- FOR AVX512
//constexpr auto AVX_FLOAT_COUNT = 16u;

using namespace std;

//#########################################################
void filter1(int Nx, int Nh, float* x,float* c,float* y) {
  //#pragma omp parallel for
  for (int i = 0; i < Nx; i++) {
    y[i] = 0.f;
    for (int j = 0; j < Nh; j++) {
      y[i] += x[i + j] * c[j];
    }
  }
  return;
}
//#########################################################
void filter2(int Nx, int Nh, float* x,float* c,float* y) {

  for (int i = 0; i < Nx; i++) {
    y[i] = 0.f;

    for (int j = 0; j < Nh; j += 4) {
      y[i] += x[i + j] * c[j] + 
              x[i + j + 1] * c[j + 1] +
              x[i + j + 2] * c[j + 2] + 
              x[i + j + 3] * c[j + 3];
    }
  }
  return;
}
//#########################################################
void filter3(int Nx, int Nh, float* x,float* c,float* y) {
  
  // A fixed-size array to move the data from registers into
  std::array<float, AVX_FLOAT_COUNT> outStore;

  for (int i = 0; i < Nx; i++) {
    auto outChunk = _mm256_setzero_ps();

    for (int j = 0u; j < Nh; j += AVX_FLOAT_COUNT) {
      auto xChunk = _mm256_loadu_ps(x + i + j);
      auto cChunk = _mm256_loadu_ps(c + j);
      auto temp = _mm256_mul_ps(xChunk, cChunk);
      outChunk = _mm256_add_ps(outChunk, temp);
    }
    _mm256_storeu_ps(outStore.data(), outChunk);

    y[i] = std::accumulate(outStore.begin(), outStore.end(), 0.f);
  }

  return;
}
//#########################################################
void filter4(int Nx, int Nh, float* x,float* c,float* y) {

  for (int i = 0; i < Nx; i +=4) {
    y[i] = 0.f;
    y[i + 1] = 0.f;
    y[i + 2] = 0.f;
    y[i + 3] = 0.f;

    for (int j = 0; j < Nh; j++) {
      y[i] += x[i + j] * c[j];
      y[i + 1] += x[i + j + 1] * c[j];
      y[i + 2] += x[i + j + 2] * c[j];
      y[i + 3] += x[i + j + 3] * c[j];
    }
  }
  return;
}
//#########################################################
void filter5(int Nx, int Nh, float* x,float* c,float* y) {

  #pragma omp parallel for
  for (int i = 0; i < Nx; i += AVX_FLOAT_COUNT) {
    auto yChunk = _mm256_setzero_ps();

    for (auto j = 0u; j < Nh; ++j) {
      auto xChunk = _mm256_loadu_ps(x + i + j);

      //-- this loads just one c but copies to all registers
      auto cChunk = _mm256_set1_ps(c[j]);

      // Element-wise multiplication
      auto temp = _mm256_mul_ps(xChunk, cChunk);

      // Add to the accumulators
      yChunk = _mm256_add_ps(yChunk, temp);
    }

    // Store 8 computed values in the result vector
    _mm256_storeu_ps(y + i, yChunk);
  }
}
//#########################################################
//#########################################################
//#########################################################
int main() {
    int Nx = 65536*4;
    int Nh = 16*4;
    float *x = new float[Nx];
    float *y = new float[Nx];
    float *c = new float[Nh];

    //-- initialization
    for (int i=0;i<Nx;i++) x[i] = i;
    for (int i=0;i<Nh;i++) {
        c[i] = 0;
        if (i==Nh-1) c[i] = 1.0;
    }
    //-- initialization

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (int k=0;k<100;k++) {
        filter1(Nx,Nh,&x[0],&c[0],&y[0]);
    }

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/100 << "[µs]" << std::endl;



  //for (int i=0;i<Nx;i++) printf("%f, ",y[i]);
  //printf("\n");

    delete [] x;
    delete [] y;
    delete [] c;
    
}
//#########################################################
