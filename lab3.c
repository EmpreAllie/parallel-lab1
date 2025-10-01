#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <mkl/mkl.h>

#ifdef _OPENMP
#include <omp.h>
#endif


int get_variant(int A, int B);
double generate_uniform(unsigned int* seed, double min, double max);
double* generate_m1(int A, int N, unsigned int* seed_ptr);
double* generate_m2(int A, int N, unsigned int* seed_ptr);
void print_array(double* arr, int size, FILE* file);
void map_m1(double* m1, int size);
void map_m2(double* m2, int size);
void merge(double* m1, double* m2, int size_m2);
void comb_sort(double* m2, int size);
double reduce(double* m2, int size, double target);
double find_min_above_zero(double* m2, int size);


int main(int argc, char* argv[]) {

    if (argc < 3) {
        printf("Запуск: %s <Размер_массива> <Число_потоков>\n", argv[0]);
        return 1;
    }

    FILE* file = fopen("output", "w");

    // информация о пользователе
    int l_name = 10;
    int l_surname = 7;
    int l_patr = 9;

    int A = l_name * l_patr * l_surname;
/*
    int B1 = 7, B2 = 8, B3 = 6, B4 = 6;
    int var1 = get_variant(A, B1);
    int var2 = get_variant(A, B2);
    int var3 = get_variant(A, B3);
    int var4 = get_variant(A, B4);
//    printf("\tИнформация о пользователе:\n");
//    printf("A = %d", A);
//    printf("\nВарианты в таблицах:\n\t1: %d\n\t2: %d\n\t3: %d\n\t4: %d\n", var1, var2, var3, var4);
*/
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    mkl_set_num_threads(M);

    int i;
    long delta_ms;
    struct timeval T1, T2;
    gettimeofday(&T1, NULL);

    unsigned int seed = time(NULL);
    for (i = 0; i < 100; i++) {

    // 1. Generate
        double* m1 = generate_m1(A, N, &seed);
        double* m2 = generate_m2(A, N, &seed);

    // 2. Map
        map_m1(m1, N);
        map_m2(m2, N/2);

    // 3. Merge
        merge(m1, m2, N/2);

    // 4. Sort
        comb_sort(m2, N/2);

    // 5. Reduce
        double target = find_min_above_zero(m2, N/2);
        double x = reduce(m2, N/2, target);

        fprintf(file, "\tИтерация № %d\nМассив m1:\n", i + 1);
        print_array(m1, N, file);
        fprintf(file, "Массив m2:\n");
        print_array(m2, N/2, file);
        fprintf(file, "Итоговое число X: %.2lf\n\n\n", x);

        free(m1);
        free(m2);
    }

    gettimeofday(&T2, NULL);
    delta_ms = (T2.tv_sec - T1.tv_sec) * 1000 + (T2.tv_usec - T1.tv_usec) / 1000;
   // printf("\nN = %d\nВремя выполнения программы: %ld\n", N, delta_ms);
    printf("%ld\n", delta_ms);

    fclose(file);

    return 0;
}

int get_variant(int A, int B) {
    return 1 + ((A % 47) % B);
}


double generate_uniform(unsigned int* seed, double min, double max) {
    return min + (max - min) * (rand_r(seed) / (RAND_MAX + 1.0));
}



double* generate_m1(int A, int N, unsigned int* seed_ptr) {
    double* m1 = (double*) malloc(N * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(A, N, m1, seed_ptr)
#endif

    for (int i = 0; i < N; i++) {
        unsigned int seed_priv;

    #ifdef _OPENMP
        seed_priv = *seed_ptr + (unsigned int) (i + omp_get_thread_num());
    #else
        seed_priv = *seed_ptr + (unsigned int) i;
    #endif

        m1[i] = generate_uniform(&seed_priv, 1.0, (double) A);
    }

    return m1;
}

double* generate_m2(int A, int N, unsigned int* seed_ptr) {
    double* m2 = (double*) malloc((N/2) * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(A, N, m2, seed_ptr)
#endif

    for (int i = 0; i < N/2; i++) {
        unsigned int seed_priv;

    #ifdef _OPENMP
        seed_priv = *seed_ptr + (unsigned int)(i + omp_get_thread_num());
    #else
        seed_priv = *seed_ptr + (unsigned int) i;
    #endif

        m2[i] = generate_uniform(&seed_priv, (double) A, 10.0 * A);
    }

    return m2;
}


void print_array(double* arr, int size, FILE* file) {
    for (int i = 0; i < size; i++) {
        fprintf(file, "%.2lf ", arr[i]);
    }
    fprintf(file, "\n\n");
}

void map_m1(double* m1, int size) {
    double e = M_E;
    double* e_arr = (double*) malloc(size * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(e_arr, e, size)
#endif

    for (int i = 0; i < size; i++) {
        e_arr[i] = e;
    }

    double* temp = (double*) malloc(size * sizeof(double));

    vdDiv(size, m1, e_arr, temp);
    vdCbrt(size, temp, m1);

    free(temp);
}

void map_m2(double* m2, int size) {
    double* m2_copy = (double*) malloc(size * sizeof(double));
    memcpy(m2_copy, m2, size * sizeof(double));

    // actual size equals (size-1) because the first element remains unchanged
    double* args = (double*) malloc((size - 1) * sizeof(double));
    double* tans = (double*) malloc((size - 1) * sizeof(double));

    // filling the array of tan() function arguments that will later be used in vdTan( . . . );
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(args, m2, m2_copy, size)
#endif
    for (int i = 1; i < size; i++) {
        args[i - 1] = m2[i] + m2_copy[i - 1];
    }

    // vectorized tangent function
    vdTan(size - 1, args, tans);

    vdInv(size - 1, tans, m2);
    vdAbs(size - 1, m2, m2);

    // we're leaving this loop as it is, because parallelization is impossible here
    for (int i = size - 1; i < 1; i--) {
        m2[i] = m2[i - 1];
    }
    m2[0] = 0.0;

    free(m2_copy);
    free(args);
    free(tans);
}


void merge(double* m1, double* m2, int size_m2) {
    vdDiv(size_m2, m1, m2, m2); // m2[i] = m1[i] / m2[i]
}


void comb_sort(double* m2, int size) {
    int gap = size;
    double k = 1.3;
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
                double temp = m2[i];
                m2[i] = m2[i + gap];
                m2[i + gap] = temp;
                sorted = false;
            }
        }
    }
}

double find_min_above_zero(double* m2, int size) {
    double target = DBL_MAX;
    for (int i = 0; i < size; i++) {
        if (m2[i] < target && m2[i] > 0) {
            target = m2[i];
        }
    }
    return target;
}

double reduce(double* m2, int size, double target) {
    double sum = 0;

    for (int i = 0; i < size; i++) {
        int integ = (int) (m2[i] / target);
        if (integ % 2 == 0) {
            sum += sin(m2[i]);
        }
    }
    return sum;
}

