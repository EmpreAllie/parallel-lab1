#!/bin/bash

# Генерация размеров массивов: 400 + 19960*i до <=200000
start=400
step=19960
max=200000

sizes=()
i=0
while true; do
    size=$((start + step * i))
    if [ $size -gt $max ]; then
        break
    fi
    sizes+=($size)
    ((i++))
done

# Исполняемые файлы
executables=("lab1-seq" "lab1-par-1" "lab1-par-2" "lab1-par-4" "lab1-par-8")

# Количество повторов
repeats=10

# Итоговый файл
output_file="results_avg.txt"
> "$output_file"  # очищаем файл перед записью

for size in "${sizes[@]}"; do
    echo "Размер массивов: $size" | tee -a "$output_file"
    
    declare -A avg_times  # массив для хранения средних значений каждой версии

    for exe in "${executables[@]}"; do
        echo "Программа: $exe" | tee -a "$output_file"

        # временный файл для хранения всех запусков текущей версии
        temp_file="temp_times.txt"
        > "$temp_file"

        # повторяем запуск программы
        for i in $(seq 1 $repeats); do
            ./$exe $size >> "$temp_file" 2>/dev/null
        done

        # усредняем все значения delta_ms, как выводит программа
        avg=$(awk '{sum+=$1; count++} END{if(count>0) print sum/count; else print "NA"}' "$temp_file")
        echo "Среднее время: $avg ms" | tee -a "$output_file"

        # сохраняем для вычисления ускорения
        avg_times[$exe]=$avg
    done

    # вычисляем ускорение параллельных версий относительно seq
    t_seq=${avg_times["lab1-seq"]}
    for exe in "${executables[@]:1}"; do
        t_par=${avg_times[$exe]}
        if [[ "$t_seq" != "NA" && "$t_par" != "NA" && "$t_par" != 0 ]]; then
            speedup=$(echo "scale=2; $t_seq / $t_par" | bc)
            echo "Ускорение $exe: $speedup" | tee -a "$output_file"
        fi
    done

    echo "-----------------------------" | tee -a "$output_file"
done
