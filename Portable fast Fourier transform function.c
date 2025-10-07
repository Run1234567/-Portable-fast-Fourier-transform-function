#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* 圆周率常量定义 */
#endif

/* 复数结构体定义 */
struct ComplexNumber {
    float real_part;      /* 复数的实部 */
    float imag_part;      /* 复数的虚部 */
};

/* 
 * 函数功能：复数乘法运算
 * 参数说明：
 *   - first_complex: 第一个复数
 *   - second_complex: 第二个复数
 * 返回值：两个复数的乘积
 */
struct ComplexNumber complex_multiply(struct ComplexNumber first_complex, 
                                     struct ComplexNumber second_complex) {
    struct ComplexNumber result_complex;
    
    /* 复数乘法公式：(a+bi) * (c+di) = (ac - bd) + (ad + bc)i */
    result_complex.real_part = first_complex.real_part * second_complex.real_part - 
                              first_complex.imag_part * second_complex.imag_part;
    result_complex.imag_part = first_complex.real_part * second_complex.imag_part + 
                              first_complex.imag_part * second_complex.real_part;
    
    return result_complex;
}

/* 
 * 函数功能：快速傅立叶变换 (FFT)
 * 参数说明：
 *   - complex_signal: 输入信号数组（复数形式）
 *   - fft_size: FFT变换的点数，必须是2的幂次方
 * 说明：该函数使用原地计算，会直接修改输入数组
 */
void fast_fourier_transform(struct ComplexNumber *complex_signal, int fft_size) {
    int loop_index_i, loop_index_j, loop_index_k, butterfly_stage, total_stages;
    int butterfly_size, half_butterfly_size, pair_index;
    struct ComplexNumber twiddle_factor, temp_complex;
    
    /* 检查FFT点数是否为2的幂次方 */
    total_stages = 0;
    int temp_size = fft_size;
    while (temp_size > 1) {
        if (temp_size % 2 != 0) {
            fprintf(stderr, "错误：FFT点数必须是2的幂次方。\n");
            return;
        }
        temp_size /= 2;
        total_stages++;
    }

    /* 第一步：位反转重新排列输入数据 */
    loop_index_j = 0;
    for (loop_index_i = 0; loop_index_i < fft_size - 1; loop_index_i++) {
        /* 交换位置i和j的数据 */
        if (loop_index_i < loop_index_j) {
            temp_complex = complex_signal[loop_index_i];
            complex_signal[loop_index_i] = complex_signal[loop_index_j];
            complex_signal[loop_index_j] = temp_complex;
        }
        
        /* 计算下一个位反转索引 */
        loop_index_k = fft_size >> 1;  /* 相当于除以2 */
        while (loop_index_k <= loop_index_j) {
            loop_index_j -= loop_index_k;
            loop_index_k >>= 1;        /* 相当于除以2 */
        }
        loop_index_j += loop_index_k;
    }

    /* 第二步：蝶形运算 */
    for (butterfly_stage = 1; butterfly_stage <= total_stages; butterfly_stage++) {
        butterfly_size = 1 << butterfly_stage;           /* 计算当前级的蝶形大小，2^butterfly_stage */
        half_butterfly_size = butterfly_size >> 1;       /* 蝶形大小的一半 */
        
        /* 计算每个蝶形组的旋转因子 */
        for (loop_index_j = 0; loop_index_j < half_butterfly_size; loop_index_j++) {
            /* 计算旋转因子（twiddle factor）W_N^k = e^(-j*2πk/N) */
            double angle_radian = -2.0 * M_PI * loop_index_j / butterfly_size;
            twiddle_factor.real_part = cos(angle_radian);
            twiddle_factor.imag_part = sin(angle_radian);
            
            /* 对每个蝶形进行计算 */
            for (loop_index_i = loop_index_j; loop_index_i < fft_size; loop_index_i += butterfly_size) {
                pair_index = loop_index_i + half_butterfly_size;
                
                /* 蝶形运算：下半部分 = 上半部分 - 旋转因子 * 下半部分 */
                temp_complex = complex_multiply(complex_signal[pair_index], twiddle_factor);
                
                complex_signal[pair_index].real_part = complex_signal[loop_index_i].real_part - temp_complex.real_part;
                complex_signal[pair_index].imag_part = complex_signal[loop_index_i].imag_part - temp_complex.imag_part;
                
                /* 蝶形运算：上半部分 = 上半部分 + 旋转因子 * 下半部分 */
                complex_signal[loop_index_i].real_part += temp_complex.real_part;
                complex_signal[loop_index_i].imag_part += temp_complex.imag_part;
            }
        }
    }
}

/* 主函数 */
int main() {
    const int FFT_POINTS = 64;  /* FFT点数，必须是2的幂次方 */
    struct ComplexNumber input_signal[FFT_POINTS];  /* 输入信号数组 */
    int index_i;  /* 循环索引变量 */
    
    /* 生成测试信号：10Hz正弦波 + 0.5 * 20Hz正弦波 */
    for (index_i = 0; index_i < FFT_POINTS; index_i++) {
        /* 生成包含10Hz和20Hz正弦波的测试信号 */
        input_signal[index_i].real_part = sin(2 * M_PI * 8 * index_i / FFT_POINTS) + 
                                        0.5 * sin(2 * M_PI * 20 * index_i / FFT_POINTS);
        input_signal[index_i].imag_part = 0.0;  /* 虚部设为0，表示实信号 */
    }
    
    /* 执行FFT变换 */
    fast_fourier_transform(input_signal, FFT_POINTS);
    
    /* 计算并打印幅度谱 */
    printf("频率分量幅度:\n");
    for (index_i = 0; index_i < FFT_POINTS / 2; index_i++) {  
        /* 只显示前一半频率bin，因为实信号的频谱是对称的 */
        float magnitude_value = sqrt(input_signal[index_i].real_part * input_signal[index_i].real_part + 
                                   input_signal[index_i].imag_part * input_signal[index_i].imag_part);
        printf("频率 bin %d: %.3f\n", index_i, magnitude_value);
    }
    
    return 0;  /* 程序正常结束 */
}
