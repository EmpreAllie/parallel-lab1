// hello world
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

int get_variant(int A, int B);
float generate_uniform(unsigned int* seed, float min, float max);
float* generate_m1(int A, int N, unsigned int* seed_ptr);
float* generate_m2(int A, int N, unsigned int* seed_ptr);
void print_array(float* arr, int size, FILE* file);
void map_m1(float* m1, int size);
void map_m2(float* m2, int size);
void merge(float* m1, float* m2, int size_m2);
void comb_sort(float* m2, int size);
float reduce(float* m2, int size, float target);
float find_min_above_zero(float* m2, int size);


int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Запуск: %s <N>\n", argv[0]);
        return 1;
    }

    FILE* file = fopen("output", "w");

    // информация о пользователе
    int l_name = 10;
    int l_surname = 7;
    int l_patr = 9;

    int A = l_name * l_patr * l_surname;

    int B1 = 7, B2 = 8, B3 = 6, B4 = 6;
    int var1 = get_variant(A, B1);
    int var2 = get_variant(A, B2);
    int var3 = get_variant(A, B3);
    int var4 = get_variant(A, B4);
    printf("\tИнформация о пользователе:\n");
    printf("A = %d", A);
    printf("\nВарианты в таблицах:\n\t1: %d\n\t2: %d\n\t3: %d\n\t4: %d\n", var1, var2, var3, var4);
    int N = atoi(argv[1]);

    int i;
    long delta_ms;
    struct timeval T1, T2;
    gettimeofday(&T1, NULL);

    unsigned int seed = time(NULL);
    for (i = 0; i < 100; i++) {

    // 1. Generate
        float* m1 = generate_m1(A, N, &seed);
        float* m2 = generate_m2(A, N, &seed);

    // 2. Map
        map_m1(m1, N);
        map_m2(m2, N/2);

    // 3. Merge
        merge(m1, m2, N/2);

    // 4. Sort
        comb_sort(m2, N/2);

    // 5. Reduce
        float target = find_min_above_zero(m2, N/2);
        float x = reduce(m2, N/2, target);

        fprintf(file, "\tИтерация № %d\nМассив m1:\n", i + 1);
        print_array(m1, N, file);
        fprintf(file, "Массив m2:\n");
     print_array(m2, N/2, file);
        fprintf(file, "Итоговое число X: %.2f\n\n\n", x);

        free(m1);
        free(m2);
    }

    gettimeofday(&T2, NULL);
    delta_ms = (T2.tv_sec - T1.tv_sec) * 1000 + (T2.tv_usec - T1.tv_usec) / 1000;
    printf("\nN = %d\nМиллисекунд прошло: %ld\n", N, delta_ms);

    fclose(file);

    return 0;
}

int get_variant(int A, int B) {
    return 1 + ((A % 47) % B);
}


float generate_uniform(unsigned int* seed, float min, float max) {
    return min + (max - min) * (rand_r(seed) / (RAND_MAX + 1.0));
}



float* generate_m1(int A, int N, unsigned int* seed_ptr) {
    float* m1 = (float*) malloc(N * sizeof(float));

    for (int i = 0; i < N; i++) {
        m1[i] = generate_uniform(seed_ptr, 1.0, (float) A);
    }

    return m1;
}

float* generate_m2(int A, int N, unsigned int* seed_ptr) {
    float* m2 = (float*) malloc(N/2 * sizeof(float));

    for (int i = 0; i < N/2; i++) {
     m2[i] = generate_uniform(seed_ptr, (float) A, 10.0 * A);
    }

    return m2;
}


void print_array(float* arr, int size, FILE* file) {
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.2f ", arr[i]);
    }
    fprintf(file, "\n\n");
}

void map_m1(float* m1, int size) {
    for (int i = 0; i < size; i++) {
        m1[i] = (float) cbrt(m1[i] / M_E);
    }
}

void map_m2(float* m2, int size) {
    float* m2_copy = (float*) malloc(size * sizeof(float));
    memcpy(m2_copy, m2, size * sizeof(float));

    for (int i = 1; i < size; i++) {
        float prev = m2_copy[i - 1];
        float sum = m2[i] + prev;
        float tan_val = tan(sum);

        if (fabs(tan_val) < 1e-10) {
            m2[i] = 0.0f;
        }
        else {
            m2[i] = fabs(1.0 / tan_val);
        }

//        m2[i] = fabs(1.0 / tan(m2[i] + m2_copy[i - 1]));
    }

    free(m2_copy);
}


void merge(float* m1, float* m2, int size_m2) {
    for (int i = 0; i < size_m2; i++) {
        m2[i] = m1[i] / m2[i];
    }
}


void comb_sort(float* m2, int size) {
    int gap = size;
    float k = 1.3;
    bool sorted = false;

    while(!sorted) {
        gap = (int) (gap / k);

        // последний проход
        if (gap <= 1) {
            gap = 1;
            sorted = true;
        }

        for (int i = 0; i + gap < size; i++) {
            if (m2[i] > m2[i + gap]) {
                // swap
                float temp = m2[i];
                m2[i] = m2[i + gap];
                m2[i + gap] = temp;
                sorted = false;
            }
        }
    }
}

float find_min_above_zero(float* m2, int size) {
    float target = FLT_MAX;
    for (int i = 0; i < size; i++) {
        if (m2[i] < target && m2[i] > 0) {
            target = m2[i];
        }
    }
    return target;
}

float reduce(float* m2, int size, float target) {
    float sum = 0;

    for (int i = 0; i < size; i++) {
        int integ = (int) (m2[i] / target);
        if (integ % 2 == 0) {
            sum += sin(m2[i]);
        }
    }
    return sum;
}

