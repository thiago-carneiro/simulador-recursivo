//
//  main.c
//  Recursive Simulator
//
//  Created by Thiago Pacheco on 02/01/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
//  #include <png.h>
//  #include <glib.h>

#define LOSS_PER_DEPTH (1 - 5e-4)
#define MIN_SIGNAL 1e-6

typedef struct
{
    float depth;
    float speed;
    float rho;
    float impedance;
    float time;
} layer;

typedef struct
{
    unsigned int n_layers;
    layer *layers;
    float total_depth;
    float *reflectivity_downwards;
    float *transmissibility_downwards;
    float *reflectivity_upwards;
    float *transmissibility_upwards;
    float *times;
} geological_model_1d;

typedef struct
{
    geological_model_1d *columns;
    float depth;
    unsigned int nx;
} geological_model_2d;

void init_layer(layer *l, float depth, float speed, float rho)
{
    l->depth = depth;
    l->speed = speed;
    l->rho = rho;
    l->impedance = speed * rho;
    l->time = depth / speed;
}

void init_geological_model_1d(geological_model_1d *model,
                              unsigned int number_layers, layer *layers)
{
    model->n_layers = number_layers;
    model->layers = layers;
    model->reflectivity_downwards = (float *)malloc((number_layers) * sizeof(float));
    model->transmissibility_downwards = (float *)malloc((number_layers) * sizeof(float));
    model->reflectivity_upwards = (float *)malloc((number_layers) * sizeof(float));
    model->transmissibility_upwards = (float *)malloc((number_layers) * sizeof(float));
    model->times = (float *)malloc((number_layers) * sizeof(float));
    model->total_depth = 0;
    float time = 0;
    for (unsigned int i = 0; i < number_layers - 1; i++)
    {
        float z1 = layers[i].impedance;
        float z2 = layers[i + 1].impedance;
        time += layers[i].time;
        model->reflectivity_downwards[i] = (z2 - z1) / (z1 + z2);
        model->transmissibility_downwards[i] = 2 * z1 / (z1 + z2);
        model->reflectivity_upwards[i] = (z1 - z2) / (z1 + z2);
        model->transmissibility_upwards[i] = 2 * z2 / (z1 + z2);
        model->times[i] = time;
        model->total_depth += layers[i].depth;
    }
    model->reflectivity_downwards[number_layers - 1] = 0;
    model->transmissibility_downwards[number_layers - 1] = 1;
    model->reflectivity_upwards[number_layers - 1] = 0;
    model->transmissibility_upwards[number_layers - 1] = 1;
    model->times[number_layers - 1] = time + layers[number_layers - 1].time;
    model->total_depth += layers[number_layers - 1].depth;
}

void free_geological_model_1d(geological_model_1d *model)
{
    free(model->reflectivity_downwards);
    free(model->transmissibility_downwards);
    free(model->reflectivity_upwards);
    free(model->transmissibility_upwards);
    free(model->times);
}

void init_geological_model_2d(geological_model_2d *model, unsigned int nx,
                              geological_model_1d *columns)
{
    model->columns = columns;
    model->nx = nx;
    model->depth = 0;
    for (unsigned int i = 0; i < nx; i++)
    {
        model->depth =
            (model->depth > columns[i].total_depth) ? model->depth : columns[i].total_depth;
    }
}

void free_geological_model_2d(geological_model_2d *model)
{
    for (unsigned int x = 0; x < model->nx; x++)
    {
        free_geological_model_1d(model->columns + x);
    }
}

geological_model_2d new_senoidal_model(unsigned int nx, unsigned int nz,
                                       float min_depth, float max_depth,
                                       float min_speed, float max_speed,
                                       float min_rho, float max_rho)
{
    float *speed = (float *)malloc(nx * nz * sizeof(float));
    float *rho = (float *)malloc(nx * nz * sizeof(float));
    for (unsigned int x = 0; x < nx; x++)
    {
        for (unsigned int z = 0; z < nz; z++)
        {
            speed[x * nz + z] = min_speed;
            rho[x * nz + z] = min_rho;
        }
    }
    float current_z = nz;
    do
    {
        float depth = ((rand() / (float)RAND_MAX) * (max_depth - min_depth)) + min_depth;
        float amplitude = (rand() / (float)RAND_MAX) * max_depth;
        float k = M_PI / ((rand() / (float)RAND_MAX) * 50);
        float phase = 2 * M_PI * (rand() / (float)RAND_MAX);
        float layer_speed = ((rand() / (float)RAND_MAX) * (max_speed - min_speed)) + min_speed;
        float layer_rho = ((rand() / (float)RAND_MAX) * (max_rho - min_rho)) + min_rho;
        current_z -= depth;
        for (unsigned int x = 0; x < nx; x++)
        {
            int layer_bottom =
                (int)(current_z + amplitude * sin(k * x + phase));
            layer_bottom = (layer_bottom < 0) ? 0 : layer_bottom;
            layer_bottom = ((unsigned int)layer_bottom > nz) ? nz : layer_bottom;
            for (unsigned int z = 0; z < (unsigned int)layer_bottom; z++)
            {
                speed[x * nz + z] = layer_speed;
                rho[x * nz + z] = layer_rho;
            }
        }
    } while (current_z > 0);
    geological_model_1d *columns =
        (geological_model_1d *)malloc(nx * sizeof(geological_model_1d));
    for (unsigned int x = 0; x < nx; x++)
    {
        unsigned int z = 0;
        unsigned int n_layers = 0;
        while (z < nz)
        {
            float current_speed = speed[x * nz + z];
            float current_rho = rho[x * nz + z];
            unsigned int current_depth = 0;
            while ((z + current_depth < nz) &&
                   (current_speed == speed[x * nz + z + current_depth]) &&
                   (current_rho == rho[x * nz + z + current_depth]))
            {
                current_depth++;
            }
            n_layers++;
            z += current_depth; // +1?
        }
        z = 0;
        unsigned int i = 0;
        layer *layers = (layer *)malloc(n_layers * sizeof(layer));
        while (z < nz)
        {
            float current_speed = speed[x * nz + z];
            float current_rho = rho[x * nz + z];
            unsigned int current_depth = 0;
            while ((z + current_depth < nz) &&
                   (current_speed == speed[x * nz + z + current_depth]) &&
                   (current_rho == rho[x * nz + z + current_depth]))
            {
                current_depth++;
            };
            init_layer(layers + i, current_depth, current_speed, current_rho);
            i++;
            z += current_depth; // +1?
        }
        init_geological_model_1d(columns + x, n_layers, layers);
    }
    geological_model_2d model;
    init_geological_model_2d(&model, nx, columns);
    return model;
}

typedef struct record
{
    float time;
    float amplitude;
    struct record *next;
} record;

typedef struct
{
    unsigned int number_records;
    record *first;
    record *last;
} registry;

registry *join_registries(registry *a, registry *b)
{
    if (!a)
        return b;
    if (!b)
        return a;
    a->number_records += b->number_records;
    b->last->next = a->first;
    a->first = b->first;
    return a;
}

typedef enum
{
    UPWARD,
    DOWNWARD
} direction;

registry *pulse(geological_model_1d *model,
                unsigned int i, float base_time, float signal, direction way,
                float max_time, unsigned int remove_multiples)
{
    registry *downward = NULL;
    registry *upward = NULL;
    registry *ret = NULL;
    signal *= pow(LOSS_PER_DEPTH, model->layers[i].depth);
    float time_i = base_time + model->layers[i].time;
    float return_time;
    if (way == DOWNWARD)
    {
        return_time = time_i + model->times[i];
    }
    else
    {
        return_time = base_time + model->times[i];
    }
    if ((return_time < max_time) && (signal > MIN_SIGNAL))
    {
        if ((way == DOWNWARD) && (i < model->n_layers))
        {
            downward = pulse(model, i + 1, time_i,
                             signal * model->transmissibility_downwards[i],
                             DOWNWARD, max_time, remove_multiples);
            upward = pulse(model, i, time_i, signal * model->reflectivity_downwards[i],
                           UPWARD, max_time, remove_multiples);
        }
        if (way == UPWARD)
        {
            if (i == 0)
            {
                ret = (registry *)malloc(sizeof(registry));
                ret->first = (record *)malloc(sizeof(record));
                ret->last = ret->first;
                ret->number_records = 1;
                ret->first->time = time_i;
                ret->first->amplitude = signal;
                ret->first->next = NULL;
                if (!remove_multiples)
                {
                    downward = pulse(model, i, time_i, -signal, DOWNWARD, max_time, remove_multiples);
                }
            }
            else
            {
                upward = pulse(model, i - 1, time_i, signal * model->transmissibility_upwards[i - 1],
                               UPWARD, max_time, remove_multiples);
                if (!remove_multiples)
                {
                    downward = pulse(model, i, time_i, signal * model->reflectivity_upwards[i - 1],
                                     DOWNWARD, max_time, remove_multiples);
                }
            }
        }
    }
    ret = join_registries(ret, downward);
    if ((downward) && (downward != ret))
        free(downward);
    ret = join_registries(ret, upward);
    if ((upward) && (upward != ret))
        free(upward);
    return ret;
}

float *measure_signals(geological_model_1d *model,
                       float max_time, unsigned int remove_multiples, unsigned int nz)
{
    registry *m = pulse(model, 0, 0., 1., DOWNWARD, max_time, remove_multiples);
    float *signals = (float *)malloc(nz * sizeof(float));
    for (unsigned int i = 0; i < nz; i++)
        signals[i] = 0.;
    float bins_per_time = nz / max_time;
    record *signal = NULL;
    if (m)
        signal = m->first;
    while (signal)
    {
        unsigned int bin = floor(signal->time * bins_per_time);
        signals[bin] += signal->amplitude;
        signal = signal->next;
    }
    free(m);
    return signals;
}

float *acquire_seismic_image(geological_model_2d *model,
                             float max_time, unsigned int remove_multiples, unsigned int nz)
{
    float *seismic_image = (float *)malloc(model->nx * nz * (sizeof(float)));
    for (unsigned int x = 0; x < model->nx; x++)
    {
        float *column_measurement = measure_signals(model->columns + x, max_time,
                                                    remove_multiples, nz);
        memcpy(seismic_image + (x * nz), column_measurement, nz);
        free(column_measurement);
    }
    return seismic_image;
}

int main(int argc, const char *argv[])
{
    srand((unsigned int)time(NULL));
    geological_model_2d model = new_senoidal_model(512, 10000, 500., 1000., 1500., 4500., 1., 2.5);
    float *seismic_image = acquire_seismic_image(&model, 3., 0, 256);
    FILE *file = fopen("seismic_image.bin", "wb");
    if (file == NULL)
    {
        printf("Error opening file!");
        return 1;
    }
    size_t count = fwrite(seismic_image, sizeof(float), 512 * 256, file);
    if (count < 512 * 256)
    {
        printf("Error writing to file!");
        return 1;
    }
    return 0;
}
