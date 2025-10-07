#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* Բ���ʳ������� */
#endif

/* �����ṹ�嶨�� */
struct ComplexNumber {
    float real_part;      /* ������ʵ�� */
    float imag_part;      /* �������鲿 */
};

/* 
 * �������ܣ������˷�����
 * ����˵����
 *   - first_complex: ��һ������
 *   - second_complex: �ڶ�������
 * ����ֵ�����������ĳ˻�
 */
struct ComplexNumber complex_multiply(struct ComplexNumber first_complex, 
                                     struct ComplexNumber second_complex) {
    struct ComplexNumber result_complex;
    
    /* �����˷���ʽ��(a+bi) * (c+di) = (ac - bd) + (ad + bc)i */
    result_complex.real_part = first_complex.real_part * second_complex.real_part - 
                              first_complex.imag_part * second_complex.imag_part;
    result_complex.imag_part = first_complex.real_part * second_complex.imag_part + 
                              first_complex.imag_part * second_complex.real_part;
    
    return result_complex;
}

/* 
 * �������ܣ����ٸ���Ҷ�任 (FFT)
 * ����˵����
 *   - complex_signal: �����ź����飨������ʽ��
 *   - fft_size: FFT�任�ĵ�����������2���ݴη�
 * ˵�����ú���ʹ��ԭ�ؼ��㣬��ֱ���޸���������
 */
void fast_fourier_transform(struct ComplexNumber *complex_signal, int fft_size) {
    int loop_index_i, loop_index_j, loop_index_k, butterfly_stage, total_stages;
    int butterfly_size, half_butterfly_size, pair_index;
    struct ComplexNumber twiddle_factor, temp_complex;
    
    /* ���FFT�����Ƿ�Ϊ2���ݴη� */
    total_stages = 0;
    int temp_size = fft_size;
    while (temp_size > 1) {
        if (temp_size % 2 != 0) {
            fprintf(stderr, "����FFT����������2���ݴη���\n");
            return;
        }
        temp_size /= 2;
        total_stages++;
    }

    /* ��һ����λ��ת���������������� */
    loop_index_j = 0;
    for (loop_index_i = 0; loop_index_i < fft_size - 1; loop_index_i++) {
        /* ����λ��i��j������ */
        if (loop_index_i < loop_index_j) {
            temp_complex = complex_signal[loop_index_i];
            complex_signal[loop_index_i] = complex_signal[loop_index_j];
            complex_signal[loop_index_j] = temp_complex;
        }
        
        /* ������һ��λ��ת���� */
        loop_index_k = fft_size >> 1;  /* �൱�ڳ���2 */
        while (loop_index_k <= loop_index_j) {
            loop_index_j -= loop_index_k;
            loop_index_k >>= 1;        /* �൱�ڳ���2 */
        }
        loop_index_j += loop_index_k;
    }

    /* �ڶ������������� */
    for (butterfly_stage = 1; butterfly_stage <= total_stages; butterfly_stage++) {
        butterfly_size = 1 << butterfly_stage;           /* ���㵱ǰ���ĵ��δ�С��2^butterfly_stage */
        half_butterfly_size = butterfly_size >> 1;       /* ���δ�С��һ�� */
        
        /* ����ÿ�����������ת���� */
        for (loop_index_j = 0; loop_index_j < half_butterfly_size; loop_index_j++) {
            /* ������ת���ӣ�twiddle factor��W_N^k = e^(-j*2��k/N) */
            double angle_radian = -2.0 * M_PI * loop_index_j / butterfly_size;
            twiddle_factor.real_part = cos(angle_radian);
            twiddle_factor.imag_part = sin(angle_radian);
            
            /* ��ÿ�����ν��м��� */
            for (loop_index_i = loop_index_j; loop_index_i < fft_size; loop_index_i += butterfly_size) {
                pair_index = loop_index_i + half_butterfly_size;
                
                /* �������㣺�°벿�� = �ϰ벿�� - ��ת���� * �°벿�� */
                temp_complex = complex_multiply(complex_signal[pair_index], twiddle_factor);
                
                complex_signal[pair_index].real_part = complex_signal[loop_index_i].real_part - temp_complex.real_part;
                complex_signal[pair_index].imag_part = complex_signal[loop_index_i].imag_part - temp_complex.imag_part;
                
                /* �������㣺�ϰ벿�� = �ϰ벿�� + ��ת���� * �°벿�� */
                complex_signal[loop_index_i].real_part += temp_complex.real_part;
                complex_signal[loop_index_i].imag_part += temp_complex.imag_part;
            }
        }
    }
}

/* ������ */
int main() {
    const int FFT_POINTS = 64;  /* FFT������������2���ݴη� */
    struct ComplexNumber input_signal[FFT_POINTS];  /* �����ź����� */
    int index_i;  /* ѭ���������� */
    
    /* ���ɲ����źţ�10Hz���Ҳ� + 0.5 * 20Hz���Ҳ� */
    for (index_i = 0; index_i < FFT_POINTS; index_i++) {
        /* ���ɰ���10Hz��20Hz���Ҳ��Ĳ����ź� */
        input_signal[index_i].real_part = sin(2 * M_PI * 8 * index_i / FFT_POINTS) + 
                                        0.5 * sin(2 * M_PI * 20 * index_i / FFT_POINTS);
        input_signal[index_i].imag_part = 0.0;  /* �鲿��Ϊ0����ʾʵ�ź� */
    }
    
    /* ִ��FFT�任 */
    fast_fourier_transform(input_signal, FFT_POINTS);
    
    /* ���㲢��ӡ������ */
    printf("Ƶ�ʷ�������:\n");
    for (index_i = 0; index_i < FFT_POINTS / 2; index_i++) {  
        /* ֻ��ʾǰһ��Ƶ��bin����Ϊʵ�źŵ�Ƶ���ǶԳƵ� */
        float magnitude_value = sqrt(input_signal[index_i].real_part * input_signal[index_i].real_part + 
                                   input_signal[index_i].imag_part * input_signal[index_i].imag_part);
        printf("Ƶ�� bin %d: %.3f\n", index_i, magnitude_value);
    }
    
    return 0;  /* ������������ */
}
