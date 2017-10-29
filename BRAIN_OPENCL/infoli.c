/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: info@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all functions.
 * The main function allocates the necessary memory, initializes the system
 * state and runs the model calculations.
 *
 */
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

#include "infoli.h"
#include "init.h"
#include "ioFile.h"
#include <pthread.h>
#define THREADED 0

const char *getErrorString(cl_int error);

typedef unsigned long long timestamp_t;

timestamp_t tWriteFile, tWriteFileStart, tWriteFileEnd;
FILE *pOutFile;
char temp[100]; // warning: this buffer may overflow

static timestamp_t get_timestamp()
{
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

struct write_to_file_param
{
    int i;
    cl_mod_prec *cellStatePtr;
};

void *write_to_file(void *args)
{
    // write output to file
    struct write_to_file_param *data = (struct write_to_file_param *)args;
    int i = data->i;
    cl_mod_prec *cellStatePtr = data->cellStatePtr;

    if (EXTRA_TIMING)
    {
        tWriteFileStart = get_timestamp();
    }
    int index = ((i % 2) ^ 1) * IO_NETWORK_SIZE * STATE_SIZE + AXON_V;
    int j, k;
    for (j = 0; j < IO_NETWORK_DIM1; j++)
    {
        for (k = 0; k < IO_NETWORK_DIM2; k++)
        {
            writeOutputDouble(
                temp,
                cellStatePtr[index + (k * IO_NETWORK_DIM1 + j) * STATE_SIZE],
                pOutFile);
        }
    }

    writeOutput(temp, ("\n"), pOutFile);
    if (EXTRA_TIMING)
    {
        tWriteFileEnd = get_timestamp();
        tWriteFile += (tWriteFileEnd - tWriteFileStart);
    }

    pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
    pthread_t write_thread;
    struct write_to_file_param args;
    //DEBUG_PRINT(("Threaded: %d\n", THREADED));
    char *outFileName = "InferiorOlive_Output.txt";
    cl_uint i, j, k, p, q;
    int simSteps = 0;
    int simTime = 0;
    int inputFromFile = 0;
    int initSteps;
    cl_mod_prec *cellStatePtr;
    cl_mod_prec *cellVDendPtr;
    cl_mod_prec *cellVAxonPtr;
    cl_mod_prec iApp;
    int seedvar;
    timestamp_t t0, t1, usecs, tNeighbourStart, tNeighbourEnd, tComputeStart,
        tComputeEnd, tReadStart, tReadEnd, tWriteStateStart, tWriteStateEnd,
        tWriteCompStart, tWriteCompEnd, tInitStart, tInitEnd, tLoopStart, tLoopEnd;
    timestamp_t tNeighbour, tCompute, tUpdate, tRead, tInit, tLoop;
    tNeighbour = tCompute = tUpdate = tRead = tWriteFile = tInit = tLoop = 0;

    cl_event writeDone, neighbourDone, computeDone, readDone;
    cl_int status;

    t0 = get_timestamp();
    if (EXTRA_TIMING)
    {
        tInitStart = get_timestamp();
    }

    simTime = SIMTIME; // in miliseconds
    simSteps = ceil(simTime / DELTA);

    DEBUG_PRINT(("Inferior Olive Model (%d x %d cell mesh)\n", IO_NETWORK_DIM1, IO_NETWORK_DIM2));

    //Open output file
    pOutFile = fopen(outFileName, "w");
    if (pOutFile == NULL)
    {
        printf("Error: Couldn't create %s\n", outFileName);
        exit(EXIT_FAILURE);
    }
    writeOutput(temp, ("#simSteps Time(ms) Input(Iapp) Output(V_axon)\n"), pOutFile);

    //Malloc for the array of cellStates and cellCompParams
    mallocCells(&cellStatePtr, &cellVDendPtr, &cellVAxonPtr);

    //Write initial state values
    InitState(cellStatePtr, cellVDendPtr);

    //Initialize g_CaL
    //init_g_CaL(cellStatePtr);

    //-----------------------------------------------------
    // STEP 1: Discover and initialize the platforms
    //-----------------------------------------------------
    cl_uint numPlatforms = 0;
    cl_platform_id *platforms = NULL;

    //Retrieve number of platforms
    status = clGetPlatformIDs(0, NULL, &numPlatforms);

    if (status != CL_SUCCESS)
    {
        printf("Error in retrieving platform IDs\n");
        exit(-1);
    }

    //Allocate space for each platform
    platforms = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));

    // Fill in platforms with clGetPlatformIDs()
    status = clGetPlatformIDs(numPlatforms, platforms, NULL);

    if (status != CL_SUCCESS)
    {
        printf("Error in filling platform IDs\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 2: Discover and initialize the devices
    //-----------------------------------------------------
    cl_uint numDevices = 0;
    cl_device_id *devices = NULL;

    // Use clGetDeviceIDs() to retrieve the number of
    // devices present
    status = clGetDeviceIDs(
        platforms[0],
        CL_DEVICE_TYPE_ALL,
        0,
        NULL,
        &numDevices);
    if (status != CL_SUCCESS)
    {
        printf("Error in retrieving number of devices present\n");
        exit(-1);
    }

    // Allocate enough space for each device
    devices = (cl_device_id *)malloc(numDevices * sizeof(cl_device_id));

    // Fill in devices with clGetDeviceIDs()
    status |= clGetDeviceIDs(
        platforms[0],
        CL_DEVICE_TYPE_ALL,
        numDevices,
        devices,
        NULL);

    // select the device which will be used
    int device_id = 0;

    if ((status != CL_SUCCESS) || (device_id >= numDevices))
    {
        printf("error in initializing devices\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 3: Create a context
    //-----------------------------------------------------
    cl_context context = NULL;

    // Create a context using clCreateContext() and
    // associate it with the devices
    context = clCreateContext(
        NULL,
        1,
        &devices[device_id],
        NULL,
        NULL,
        &status);
    if (status != CL_SUCCESS)
    {
        printf("error in creating context\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 4: Create a command queue
    //-----------------------------------------------------
    cl_command_queue cmdQueue;

    // Create a command queue using clCreateCommandQueue(),
    // and associate it with the device you want to execute
    // on
    cmdQueue = clCreateCommandQueue(
        context,
        devices[device_id],
        0,
        &status);

    if (status != CL_SUCCESS)
    {
        printf("error in creating command queue\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 5: Create device buffers
    //-----------------------------------------------------
    cl_mem bufferCellState, bufferCellVDend, bufferCellVAxon;

    bufferCellState = clCreateBuffer(
        context,
        CL_MEM_READ_WRITE,
        IO_NETWORK_DIM1 * IO_NETWORK_DIM2 * PARAM_SIZE * sizeof(cl_mod_prec),
        NULL,
        &status);

    if (status != CL_SUCCESS)
    {
        printf("error in step 5, creating buffer for bufferCellState\n");
        exit(-1);
    }

    bufferCellVDend = clCreateBuffer(
        context,
        CL_MEM_READ_WRITE,
        IO_NETWORK_SIZE * sizeof(cl_mod_prec),
        NULL,
        &status);

    if (status != CL_SUCCESS)
    {
        printf("error in step 5, creating buffer for bufferCellVDend\n");
        exit(-1);
    }

    bufferCellVAxon = clCreateBuffer(
        context,
        CL_MEM_WRITE_ONLY,
        IO_NETWORK_SIZE * sizeof(cl_mod_prec),
        NULL,
        &status);

    if (status != CL_SUCCESS)
    {
        printf("error in step 5, creating buffer for bufferCellVAxon\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 6: Write host data to device buffers
    //-----------------------------------------------------
    status = clEnqueueWriteBuffer(
        cmdQueue,
        bufferCellState,
        CL_FALSE,
        0,
        IO_NETWORK_SIZE * PARAM_SIZE * sizeof(cl_mod_prec),
        cellStatePtr,
        0,
        NULL,
        &writeDone);

    status |= clEnqueueWriteBuffer(
        cmdQueue,
        bufferCellVDend,
        CL_FALSE,
        0,
        IO_NETWORK_SIZE * sizeof(cl_mod_prec),
        cellVDendPtr,
        0,
        NULL,
        &writeDone);

    if (status != CL_SUCCESS)
    {
        printf("error in step 6, writing data\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 7: Create and compile the program
    //-----------------------------------------------------
    /**
    ** Please check if your kernel files can be opened, 
    ** especially when running on the server
    **/
    char *computeFileName, *neighbourFileName;
    neighbourFileName = "neighbour_kernel.cl";
    computeFileName = "compute_kernel.cl";

    FILE *computeFile, *neighbourFile;
    neighbourFile = fopen(neighbourFileName, "r");
    if (neighbourFile == NULL)
    {
        printf("cannot open neighbour file\n");
        printf("current path: %s\n", neighbourFileName);
        exit(EXIT_FAILURE);
    }

    computeFile = fopen(computeFileName, "r");
    if (computeFile == NULL)
    {
        printf("cannot open compute file\n");
        printf("current path: %s\n", computeFileName);
        exit(EXIT_FAILURE);
    }

    // Get neighbour and compute filesize
    fseek(neighbourFile, 0, SEEK_END);
    size_t neighbourSize = ftell(neighbourFile);
    rewind(neighbourFile);

    fseek(computeFile, 0, SEEK_END);
    size_t computeSize = ftell(computeFile);
    rewind(computeFile);

    // read kernel sources into buffers
    char *neighbourBuffer, *computeBuffer;
    neighbourBuffer = (char *)malloc(neighbourSize + 1);
    if (neighbourBuffer == NULL)
    {
        printf("NeighbourBuffer malloc failed\n");
        exit(-1);
    }
    neighbourBuffer[neighbourSize] = '\0';
    fread(neighbourBuffer, sizeof(char), neighbourSize, neighbourFile);
    fclose(neighbourFile);

    computeBuffer = (char *)malloc(computeSize + 1);
    if (computeBuffer == NULL)
    {
        printf("NeighbourBuffer malloc failed\n");
        exit(-1);
    }
    computeBuffer[computeSize] = '\0';
    fread(computeBuffer, sizeof(char), computeSize, computeFile);
    fclose(computeFile);

    // Create programs
    cl_program neighbourProgram = clCreateProgramWithSource(
        context,
        1,
        (const char **)&neighbourBuffer,
        &neighbourSize,
        &status);
    free(neighbourBuffer);

    cl_program computeProgram = clCreateProgramWithSource(
        context,
        1,
        (const char **)&computeBuffer,
        &computeSize,
        &status);
    free(computeBuffer);

    // Compile program
    const char options[] = "-cl-std=CL1.2";
    status = clBuildProgram(
        neighbourProgram,
        1,
        &devices[device_id],
        options,
        NULL,
        NULL);
    if (status != CL_SUCCESS)
    {
        printf("error in step 7, neighbourProgram\n");
        exit(-1);
    }
    status = clBuildProgram(
        computeProgram,
        1,
        &devices[device_id],
        options,
        NULL,
        NULL);
    if (status != CL_SUCCESS)
    {
        printf("error in step 7, computeProgram, error code %d\n", status);
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 8: Create the kernel
    //----------------------------------------------------
    cl_kernel neighbourKernel = NULL;
    cl_kernel computeKernel = NULL;
    // Use clCreateKernel() to create a kernel
    neighbourKernel = clCreateKernel(neighbourProgram, "neighbour_kernel", &status);
    if (status != CL_SUCCESS)
    {
        printf("error in step 8, neighbour\n");
        exit(-1);
    }
    computeKernel = clCreateKernel(computeProgram, "compute_kernel", &status);
    if (status != CL_SUCCESS)
    {
        printf("error in step 8, compute\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 9: Set the kernel arguments
    //-----------------------------------------------------
    status = clSetKernelArg(
        neighbourKernel,
        0,
        sizeof(cl_mem),
        &bufferCellState);
    if (status != CL_SUCCESS)
    {
        printf("error in step 9.1\n");
        exit(-1);
    }
    status |= clSetKernelArg(
        neighbourKernel,
        1,
        sizeof(cl_mem),
        &bufferCellVDend);
    if (status != CL_SUCCESS)
    {
        printf("error in step 9.2\n");
        exit(-1);
    }

    status |= clSetKernelArg(
        computeKernel,
        0,
        sizeof(cl_mem),
        &bufferCellState);
    if (status != CL_SUCCESS)
    {
        printf("error in step 9.3\n");
        exit(-1);
    }
    status |= clSetKernelArg(
        computeKernel,
        1,
        sizeof(cl_mem),
        &bufferCellVDend);
    status |= clSetKernelArg(
        computeKernel,
        2,
        sizeof(cl_mem),
        &bufferCellVAxon);
    if (status != CL_SUCCESS)
    {
        printf("error in step 9.4\n");
        exit(-1);
    }

    //-----------------------------------------------------
    // STEP 10: Configure the work-item structure
    //-----------------------------------------------------

    size_t globalWorkSize[2];
    globalWorkSize[0] = IO_NETWORK_DIM1;
    globalWorkSize[1] = IO_NETWORK_DIM2;

    size_t local_work_size1 = max(1, min(IO_NETWORK_DIM1/8, 16));
    size_t local_work_size2 = max(1, min(IO_NETWORK_DIM2/8, 16));
    size_t localWorkSize[] = {local_work_size1, local_work_size2};

    if (EXTRA_TIMING)
    {
        tInitEnd = get_timestamp();
    }
    if (EXTRA_TIMING)
    {
        tLoopStart = get_timestamp();
    }

    for (i = 0; i < simSteps; i++)
    {
        //Compute one sim step for all cells
        if (i > 20000 - 1 && i < 20500 - 1)
        {
            iApp = 6;
        } // start @ 1 because skipping initial values
        else
        {
            iApp = 0;
        }

        // TODO this needs to optimized to similar in CUDA version
        // Set iApp and simulation step parameters
        /* status = clSetKernelArg(
            neighbourKernel,
            2,
            sizeof(cl_uint),
            &i); */

        status = clSetKernelArg(
            computeKernel,
            3,
            sizeof(cl_mod_prec),
            &iApp);
        /* status |= clSetKernelArg(
            computeKernel,
            3,
            sizeof(cl_uint),
            &i); */
        if (status != CL_SUCCESS)
        {
            printf("error in step 11.0\n");
            exit(-1);
        }

        if (EXTRA_TIMING)
        {
            tNeighbourStart = get_timestamp();
        }

        //-----------------------------------------------------
        // STEP 11.1: Run neighbour kernel
        //-----------------------------------------------------
        status = clEnqueueNDRangeKernel(
            cmdQueue,
            neighbourKernel,
            2,
            NULL,
            globalWorkSize,
            localWorkSize,
            1,
            &writeDone, // Wait for initial write
            &neighbourDone);
        if (status != CL_SUCCESS)
        {
            printf("error in step 11.1 neighbour kernel\n");
            exit(-1);
        }

        if (EXTRA_TIMING)
        {
            status = clWaitForEvents(1, &neighbourDone);
            tNeighbourEnd = get_timestamp();
            tNeighbour += (tNeighbourEnd - tNeighbourStart);
            tComputeStart = get_timestamp();
        }

        //-----------------------------------------------------
        // STEP 11.2: Run compute kernel
        //-----------------------------------------------------
        status = clEnqueueNDRangeKernel(
            cmdQueue,
            computeKernel,
            2,
            NULL,
            globalWorkSize,
            localWorkSize,
            1,
            &neighbourDone, // Wait for neighbour done
            &computeDone);
        if (status != CL_SUCCESS)
        {
            printf("error in step 11.2 compute kernel\n");
            exit(-1);
        }

        if (EXTRA_TIMING)
        {
            status = clWaitForEvents(1, &computeDone);
            tComputeEnd = get_timestamp();
            tCompute += (tComputeEnd - tComputeStart);
            tReadStart = get_timestamp();
        }

        if (status != CL_SUCCESS)
        {
            printf("error in loop, compute\n");
            exit(EXIT_FAILURE);
        }

        // transfer data from device to CPU
        if (WRITE_OUTPUT)
        {
            //-----------------------------------------------------
            // STEP 11.3: Read output data from device
            //-----------------------------------------------------
            clEnqueueReadBuffer(
                cmdQueue,
                bufferCellVAxon,
                CL_TRUE,
                0,
                IO_NETWORK_SIZE * sizeof(cl_mod_prec),
                cellVAxonPtr,
                1,
                &computeDone,
                NULL);
            if (status != CL_SUCCESS)
            {
                printf("error in reading data\n");
                exit(-1);
            }
        }

        if (EXTRA_TIMING)
        {
            tReadEnd = get_timestamp();
            tRead += (tReadEnd - tReadStart);
        }

        // write output to file
        if (WRITE_OUTPUT)
        {
            if (THREADED)
            {
                if (EXTRA_TIMING)
                {
                    tWriteFileStart = get_timestamp();
                }
                // Wait for previous data to be written to file
                if (i)
                    pthread_join(write_thread, NULL);

                // prepare data
                args.i = i;
                args.cellStatePtr = cellStatePtr; // Doesn't change

                // Could also be moved to pthread.
                sprintf(temp, "%d %.2f %.1f ", i + 1, i * 0.05,
                        iApp); // start @ 1 because skipping initial values
                fputs(temp, pOutFile);
                // Create thread
                pthread_create(&write_thread, NULL, write_to_file,
                               (void *)&args);
            }
            else
            {
                sprintf(temp, "%d %.2f %.1f ", i + 1, i * 0.05,
                        iApp); // start @ 1 because skipping initial values
                fputs(temp, pOutFile);
                if (EXTRA_TIMING)
                {
                    tWriteFileStart = get_timestamp();
                }
                for (j = 0; j < IO_NETWORK_SIZE; j++)
                {
                    writeOutputDouble(
                        temp, cellVAxonPtr[j],
                        pOutFile);
                }
            }
            writeOutput(temp, ("\n"), pOutFile);
            if (EXTRA_TIMING)
            {
                tWriteFileEnd = get_timestamp();
                tWriteFile += (tWriteFileEnd - tWriteFileStart);
            }
        }
    }

    // Wait for final data to be written to file
    if (WRITE_OUTPUT && THREADED)
        pthread_join(write_thread, NULL);

    if (EXTRA_TIMING)
    {
        tLoopEnd = get_timestamp();
    }

    t1 = get_timestamp();
    usecs = (t1 - t0); // / 1000000;
    DEBUG_PRINT(("%d ms of brain time in %d simulation steps\n", simTime, simSteps));
    DEBUG_PRINT((" %lld usecs real time \n", usecs));

    if (EXTRA_TIMING)
    {
        tInit = (tInitEnd - tInitStart);
        tLoop = (tLoopEnd - tLoopStart);

        //DEBUG_PRINT(("\n"));
        DEBUG_PRINT(("----------------------------------\n"));
        DEBUG_PRINT(("tInit: \t\t %lld \n", tInit));
        DEBUG_PRINT(("tLoop: \t\t %lld \n", tLoop));
        DEBUG_PRINT(("\ttNeighbour: \t %lld \n", tNeighbour));
        DEBUG_PRINT(("\ttCompute: \t %lld \n", tCompute));
        DEBUG_PRINT(("\ttRead: \t\t %lld \n", tRead));
        DEBUG_PRINT(("\ttWriteFile: \t %lld \n", tWriteFile));
        DEBUG_PRINT(("\t----------- + \n"));
        DEBUG_PRINT(("\ttSumLoop: \t %lld \n", (tWriteFile + tCompute + tNeighbour + tRead)));
        DEBUG_PRINT(("----------------------------------\n"));
        DEBUG_PRINT(("tSum: \t %lld \n", (tInit + tLoop)));
    }

    //-----------------------------------------------------
    // STEP 12: Release OpenCL resources
    //-----------------------------------------------------
    clReleaseKernel(neighbourKernel);
    clReleaseKernel(computeKernel);
    clReleaseProgram(neighbourProgram);
    clReleaseProgram(computeProgram);
    clReleaseCommandQueue(cmdQueue);
    clReleaseMemObject(bufferCellState);
    clReleaseMemObject(bufferCellVDend);
    clReleaseContext(context);

    //Free up memory and close files
    free(cellStatePtr);
    free(cellVDendPtr);
    fclose(pOutFile);

    return EXIT_SUCCESS;
}

const char *getErrorString(cl_int error)
{
    switch (error)
    {
    // run-time and JIT compiler errors
    case 0:
        return "CL_SUCCESS";
    case -1:
        return "CL_DEVICE_NOT_FOUND";
    case -2:
        return "CL_DEVICE_NOT_AVAILABLE";
    case -3:
        return "CL_COMPILER_NOT_AVAILABLE";
    case -4:
        return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5:
        return "CL_OUT_OF_RESOURCES";
    case -6:
        return "CL_OUT_OF_HOST_MEMORY";
    case -7:
        return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8:
        return "CL_MEM_COPY_OVERLAP";
    case -9:
        return "CL_IMAGE_FORMAT_MISMATCH";
    case -10:
        return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11:
        return "CL_BUILD_PROGRAM_FAILURE";
    case -12:
        return "CL_MAP_FAILURE";
    case -13:
        return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14:
        return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15:
        return "CL_COMPILE_PROGRAM_FAILURE";
    case -16:
        return "CL_LINKER_NOT_AVAILABLE";
    case -17:
        return "CL_LINK_PROGRAM_FAILURE";
    case -18:
        return "CL_DEVICE_PARTITION_FAILED";
    case -19:
        return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30:
        return "CL_INVALID_VALUE";
    case -31:
        return "CL_INVALID_DEVICE_TYPE";
    case -32:
        return "CL_INVALID_PLATFORM";
    case -33:
        return "CL_INVALID_DEVICE";
    case -34:
        return "CL_INVALID_CONTEXT";
    case -35:
        return "CL_INVALID_QUEUE_PROPERTIES";
    case -36:
        return "CL_INVALID_COMMAND_QUEUE";
    case -37:
        return "CL_INVALID_HOST_PTR";
    case -38:
        return "CL_INVALID_MEM_OBJECT";
    case -39:
        return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40:
        return "CL_INVALID_IMAGE_SIZE";
    case -41:
        return "CL_INVALID_SAMPLER";
    case -42:
        return "CL_INVALID_BINARY";
    case -43:
        return "CL_INVALID_BUILD_OPTIONS";
    case -44:
        return "CL_INVALID_PROGRAM";
    case -45:
        return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46:
        return "CL_INVALID_KERNEL_NAME";
    case -47:
        return "CL_INVALID_KERNEL_DEFINITION";
    case -48:
        return "CL_INVALID_KERNEL";
    case -49:
        return "CL_INVALID_ARG_INDEX";
    case -50:
        return "CL_INVALID_ARG_VALUE";
    case -51:
        return "CL_INVALID_ARG_SIZE";
    case -52:
        return "CL_INVALID_KERNEL_ARGS";
    case -53:
        return "CL_INVALID_WORK_DIMENSION";
    case -54:
        return "CL_INVALID_WORK_GROUP_SIZE";
    case -55:
        return "CL_INVALID_WORK_ITEM_SIZE";
    case -56:
        return "CL_INVALID_GLOBAL_OFFSET";
    case -57:
        return "CL_INVALID_EVENT_WAIT_LIST";
    case -58:
        return "CL_INVALID_EVENT";
    case -59:
        return "CL_INVALID_OPERATION";
    case -60:
        return "CL_INVALID_GL_OBJECT";
    case -61:
        return "CL_INVALID_BUFFER_SIZE";
    case -62:
        return "CL_INVALID_MIP_LEVEL";
    case -63:
        return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64:
        return "CL_INVALID_PROPERTY";
    case -65:
        return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66:
        return "CL_INVALID_COMPILER_OPTIONS";
    case -67:
        return "CL_INVALID_LINKER_OPTIONS";
    case -68:
        return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000:
        return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001:
        return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002:
        return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003:
        return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004:
        return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005:
        return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default:
        return "Unknown OpenCL error";
    }
}
