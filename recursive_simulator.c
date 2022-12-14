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
#include <zlib.h>
#include <mpi.h>
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
    free(model->layers);
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
    free(model->columns);
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
        float k = M_PI / ((rand() / (float)RAND_MAX) * 900 + 100);
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
    free(speed);
    free(rho);
    return model;
}

geological_model_2d new_parallel_model(unsigned int nx, unsigned int n_layers,
                                       float *depths, float *speeds, float *rhos)
{
    if (!n_layers)
    {
        n_layers = 4;
    }
    float depths_v[] = {1000., 1000., 1000., 1000.};
    if (!depths)
    {

        depths = depths_v;
    }
    float speeds_v[] = {1500., 2800., 4500., 3200.};
    if (!speeds)
    {
        speeds = speeds_v;
    }
    float rhos_v[] = {1., 2., 2.16, 2.6};
    if (!rhos)
    {
        rhos = rhos_v;
    }
    geological_model_1d *columns = (geological_model_1d *)malloc(nx * sizeof(geological_model_1d));
    for (unsigned int x = 0; x < nx; x++)
    {
        layer *layers = (layer *)malloc(n_layers * sizeof(layer));
        for (unsigned int l = 0; l < n_layers; l++)
            init_layer(layers + l, depths[l], speeds[l], rhos[l]);
        init_geological_model_1d(columns + x, n_layers, layers);
    }
    geological_model_2d model;
    init_geological_model_2d(&model, nx, columns);
    return model;
}

geological_model_2d new_dome_model(unsigned int nx, float *speeds, float *rhos)
{
    float speeds_v[] = {1500., 2800., 4500., 3200.};
    if (!speeds)
    {
        speeds = speeds_v;
    }
    float rhos_v[] = {1., 2., 2.16, 2.6};
    if (!rhos)
    {
        rhos = rhos_v;
    }
    geological_model_1d *columns = (geological_model_1d *)malloc(nx * sizeof(geological_model_1d));
    float ray = nx / 6.;
    for (unsigned int x = 0; x < nx; x++)
    {
        float cathetus_x = fabs(x - nx / 2.);
        float cathetus_z;
        unsigned int include_dome;
        if (cathetus_x < ray)
        {
            cathetus_z = 7. * sqrt(ray * ray - cathetus_x * cathetus_x);
            include_dome = 1;
        }
        else
        {
            cathetus_z = 0;
            include_dome = 0;
        }
        layer *layers = (layer *)malloc((3 + include_dome) * sizeof(layer));
        init_layer(layers, 1000, speeds[0], rhos[0]);
        init_layer(layers + 1, 1000 - cathetus_z, speeds[1], rhos[1]);
        if (include_dome)
        {
            init_layer(layers + 2, cathetus_z, speeds[2], rhos[2]);
        }
        init_layer(layers + 2 + include_dome, 1000, speeds[3], rhos[3]);
        init_geological_model_1d(columns + x, 3 + include_dome, layers);
    }
    geological_model_2d model;
    init_geological_model_2d(&model, nx, columns);
    return model;
}

geological_model_2d new_presalt_model(unsigned int nx, float *speeds, float *rhos)
{
    float speeds_v[] = {1500., 2800., 4500., 5000., 4500., 2000., 3200.};
    if (!speeds)
    {
        speeds = speeds_v;
    }
    float rhos_v[] = {1., 2., 2.16, 2.5, 2.16, 1.5, 2.6};
    if (!rhos)
    {
        rhos = rhos_v;
    }
    geological_model_1d *columns = (geological_model_1d *)malloc(nx * sizeof(geological_model_1d));
    float ray = nx / 6.;
    for (unsigned int x = 0; x < nx; x++)
    {
        float cathetus_x = fabs(x - nx / 2.);
        float cathetus_z;
        unsigned int include_dome;
        if (cathetus_x < ray)
        {
            cathetus_z = sqrt(ray * ray - cathetus_x * cathetus_x);
            include_dome = 1;
        }
        else
        {
            cathetus_z = 0;
            include_dome = 0;
        }
        layer *layers = (layer *)malloc((3 + 4 * include_dome) * sizeof(layer));
        init_layer(layers, 1000, speeds[0], rhos[0]);
        init_layer(layers + 1, 1000 - 7 * cathetus_z, speeds[1], rhos[1]);
        if (include_dome)
        {
            init_layer(layers + 2, 2 * cathetus_z, speeds[2], rhos[2]);
            init_layer(layers + 3, 2 * cathetus_z, speeds[3], rhos[3]);
            init_layer(layers + 4, 2 * cathetus_z, speeds[4], rhos[4]);
            init_layer(layers + 5, cathetus_z, speeds[5], rhos[5]);
        }
        init_layer(layers + 2 + 4 * include_dome, 1000, speeds[6], rhos[6]);
        init_geological_model_1d(columns + x, 3 + 4 * include_dome, layers);
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
    record *first;
    record *last;
} registry;

registry *join_registries(registry *a, registry *b)
{
    if (!a)
        return b;
    if (!b)
        return a;
    b->last->next = a->first;
    b->last = a->last;
    free(a);
    return b;
}

void free_registry(registry *a)
{
    if (!a)
        return;
    record *rec = a->first;
    record *next;
    while ((next = rec->next))
    {
        free(rec);
        rec = next;
    }
    free(rec);
    return;
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
    if ((return_time < max_time) && (fabs(signal) > MIN_SIGNAL))
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
    ret = join_registries(ret, upward);
    return ret;
}

struct parameters
{
    geological_model_1d *model;
    float max_time;
    unsigned int remove_multiples;
    unsigned int nz;
};

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
        unsigned int bin = (signal->time * bins_per_time);
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
        memcpy(seismic_image + (x * nz), column_measurement, nz * sizeof(float));
        free(column_measurement);
    }
    return seismic_image;
}

float *rasterize_geology_spd(geological_model_2d *model, unsigned int nz)
{
    float *image = (float *)malloc(model->nx * nz * sizeof(float));
    float z_ratio = nz / model->depth;
    for (unsigned int x = 0; x < model->nx; x++)
    {
        float current_z = 0.;
        for (unsigned int l = 0; l < model->columns[x].n_layers; l++)
        {
            float height = model->columns[x].layers[l].depth * z_ratio;
            for (unsigned int z = current_z; z < current_z + height; z++)
            {
                image[x * nz + z] = model->columns[x].layers[l].speed;
            }
            current_z += height;
        }
    }
    return image;
}

void batch(size_t number_simulations, unsigned int nx, unsigned int nz, char *filename)
{
    int number_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &number_procs);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t batchsize = number_simulations / number_procs;

    size_t buffer_size = 2 * batchsize * nx * nz * sizeof(float);

    uLong source_size = 2 * number_simulations * nx * nz * sizeof(float);
    if (source_size != 2 * number_simulations * nx * nz * sizeof(float))
    {
        printf("Error: buffer too large for zlib compress2.");
        exit(EXIT_FAILURE);
    }

    printf("Getting ready to generate %s\n", filename);
    geological_model_2d *models = (geological_model_2d *)malloc(batchsize * sizeof(geological_model_2d));
    if (!models)
    {
        printf("Failed to allocate memory for the geological models!\n");
        exit(EXIT_FAILURE);
    }
    float *seismic_images = (float *)malloc(buffer_size);
    if (!seismic_images)
    {
        printf("Failed to allocate memory for the seismic images!\n");
        exit(EXIT_FAILURE);
    }
    float *gathered_data = NULL;
    Bytef *compression_buffer = NULL;
    if (!rank)
    {
        gathered_data = (float *)malloc(source_size);
        if (!gathered_data)
        {
            printf("Failed to allocate memory for the results!\n");
            exit(EXIT_FAILURE);
        }
        compression_buffer = (Bytef *)malloc(source_size);
        if (!compression_buffer)
        {
            printf("Failed to allocate memory for the zlib compression!\n");
            exit(EXIT_FAILURE);
        }
    }

    struct timespec start, loop_time, last_time;
    clock_gettime(CLOCK_REALTIME, &start);
    clock_gettime(CLOCK_REALTIME, &last_time);

    for (size_t i = 0; i < batchsize; i++)
    {
        printf("%s :: Rank %d :: ", filename, rank);
        printf("Simulation %zu of %zu. ", i + 1, batchsize);
        models[i] = new_senoidal_model(nx, 10000, 500., 1000., 1500., 4500., 1., 2.5);
        float *seismic_image_complete = acquire_seismic_image(models + i, 3., 0, nz);
        memcpy(seismic_images + 2 * i * nx * nz, seismic_image_complete, nx * nz * sizeof(float));
        float *seismic_image_clean = acquire_seismic_image(models + i, 3., 1, nz);
        memcpy(seismic_images + (2 * i + 1) * nx * nz, seismic_image_clean, nx * nz * sizeof(float));
        clock_gettime(CLOCK_REALTIME, &loop_time);
        float total_elapsed = loop_time.tv_sec - start.tv_sec + (loop_time.tv_nsec - start.tv_nsec) * 1e-9;
        float instant_elapsed = loop_time.tv_sec - last_time.tv_sec + (loop_time.tv_nsec - last_time.tv_nsec) * 1e-9;
        float total_remaining = (batchsize - i) * total_elapsed / (i + 1);
        printf("Loop time: %.2fs. Running time: %.2fh. Remaining time: %.2fh.\n", instant_elapsed, total_elapsed / 3600., total_remaining / 3600.);
        free_geological_model_2d(models + i);
        free(seismic_image_complete);
        free(seismic_image_clean);
        clock_gettime(CLOCK_REALTIME, &last_time);
    }
    MPI_Gather(seismic_images, buffer_size, MPI_FLOAT, gathered_data, buffer_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    uLongf compressed_size = source_size;
    FILE *file;

    if (!rank)
    {
        printf("\nFinished simulations for %s. Now compressing.\n", filename);
        if (Z_OK != compress2((Bytef *)compression_buffer, &compressed_size, (Bytef *)gathered_data, source_size, Z_BEST_COMPRESSION))
        {
            printf("Error during compression.\n");
            exit(EXIT_FAILURE);
        }
        printf("Compression sucessful.\n");

        file = fopen(filename, "wb");
        if (file == NULL)
        {
            printf("Error opening file!\n");
            exit(EXIT_FAILURE);
        }
        fwrite(compression_buffer, 1, compressed_size, file);
        if (ferror(file))
        {
            printf("Error writing to file!\n");
            exit(EXIT_FAILURE);
        }
        fclose(file);

        printf("\tUncompressed data: %lu\n", source_size);
        printf("\tCompressed size: %lu\n", compressed_size);
        printf("\tCompression ratio: %.2f\n\n", compressed_size / (float)source_size);
        free(compression_buffer);
        free(gathered_data);
    }
    free(models);
    free(seismic_images);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    srand((unsigned int)time(NULL));

    batch(100000, 512, 256, "train.zlib");
    batch(30000, 512, 256, "test.zlib");
    batch(30000, 512, 256, "val.zlib");
    MPI_Finalize();
    return 0;
}