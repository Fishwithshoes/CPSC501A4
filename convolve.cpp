#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// int16_t
using namespace std;

char* parseCompleteFile(char *currName);
int extractIntAtIndex(char *currFile, int index);
float* intToFloat(int16_t* inArray, int size);
int16_t* floatToInt (float* inArray, int size);
float * normalFloat (float* inArray, int size);
void convolve(float x[], int N, float h[], int M, float y[], int P);
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Test tone duration in seconds  */
#define DURATION          2.0				

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

int main(int argc, char *argv[]) {
    
    if (argc != 4) {
        cout << "Please provide the correct arguements in the following order:\n"
        << "inputfile IRfile outputfile\n";
    }
    
    char *inName = argv[1];
    char *inputFileComplete;
    int16_t *inData;
    float *x;
    int N;
    
    char *IRName = argv[2];
    char *IRFileComplete;
    int16_t *IRData;
    float *h;
    int M;
    
    char *outName = argv[3];
    float *y;
    int16_t *outData;
    int P;
   
    inputFileComplete = parseCompleteFile(inName);
    int inChunkSize = extractIntAtIndex(inputFileComplete, 4);
    cout << "in ChunkSize: " << inChunkSize << endl;
    int inSubChunk1Size = extractIntAtIndex(inputFileComplete, 16);
    cout << "in subChunk1Size: " << inSubChunk1Size << endl;
    int inSubChunk1SizeCorrected = inSubChunk1Size - 16;
    cout << "in subChunk1SizeCorrected: " << inSubChunk1SizeCorrected << endl;
    memcpy(&N, inputFileComplete + (40 + inSubChunk1SizeCorrected), 4);
    cout << "N: " << N <<endl;
    inData = new int16_t[N];
    memcpy(inData, inputFileComplete + (44 + inSubChunk1SizeCorrected), N); //no vals over INT16_MAX, less -INT16_MAX
    x = new float[N];
    x = intToFloat(inData, N); // no vals over 1.0, less -1.0
    
    IRFileComplete = parseCompleteFile(IRName);
    int IRChunkSize = extractIntAtIndex(IRFileComplete, 4);
    cout << "IR ChunkSize: " << IRChunkSize << endl;
    int IRSubChunk1Size = extractIntAtIndex(IRFileComplete, 16);
    cout << "IR subChunk1Size: " << IRSubChunk1Size << endl;
    int IRSubChunk1SizeCorrected = IRSubChunk1Size - 16;
    cout << "IR subChunk1SizeCorrected: " << IRSubChunk1SizeCorrected << endl;
    memcpy(&M, IRFileComplete + (40 + IRSubChunk1SizeCorrected), 4);
        cout << "M: " << M <<endl;
    IRData = new int16_t[M];
    memcpy(IRData, IRFileComplete + (44 + IRSubChunk1SizeCorrected), M); //no vals over INT16_MAX, less -INT16_MAX
    h = new float[M];
    h = intToFloat(IRData, M); //no vals over 1.0, less -1.0
    
    P = N + M - 1;
    y = new float[P];
    
    convolve(x, N, h, M, y, P);
    y = normalFloat(y, P);
    
    outData = new int16_t[M];
    
    outData = floatToInt(y, P);
    //outData = floatToInt(h, M);
    //outData = floatToInt(x, N);
    for (int i = 0; i < P; i++) {
        cout << outData[i] << " ";
    }
    
    FILE *outputFileStream = fopen(outName, "wb");
    if (outputFileStream == NULL) {
        fprintf(stderr, "File %s cannot be opened for writing\n", outName);
    }
    
    writeWaveFileHeader(1,P,44100.0, outputFileStream);
    for (int i = 0; i < P; i++) {
        fwriteShortLSB(outData[i], outputFileStream);
    }
    
    cout << "P: " << P <<endl;
    
    fclose(outputFileStream);
    
    return 0;
}

void convolve(float x[], int N, float h[], int M, float y[], int P) {
    cout << "Convolving... " << endl;
    int n, m;
    for (int i = 0; i < P; i++) {
        y[i] = 0.0;
    }
    for (n = 0; n < N; n++) {
        for (m = 0; m < M; m++) {
            y[n+m] += x[n] * h[m];
        }
    }
}

char* parseCompleteFile(char *currName) {
    streampos size;
    char* memBlock;
    
    ifstream file (currName, ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size = file.tellg();
        memBlock = new char [size];
        file.seekg(0, ios::beg);
        file.read(memBlock,size);
        file.close();
        
        /*for (int i = 0; i < 50; i++) {
            cout << memBlock[i] << " ";
        }*/
    cout << endl;
    return memBlock;
    }
    else {
    cout << "File read error" << endl;
    }
}

int extractIntAtIndex(char *currFile, int index) {
    int returnVal;
    memcpy(&returnVal, currFile+index, sizeof(returnVal));
    return returnVal;
}

float * intToFloat (int16_t* inArray, int size) {
    float *outArray = new float[size];
    for (int i = 0; i < size; i++) {
        float f = (float) inArray[i]/(INT16_MAX + 1); //currMax;
        outArray[i] = f;
    }
    return outArray;
}

int16_t* floatToInt (float* inArray, int size) {
    int16_t *outArray = new int16_t[size];
    for (int i = 0; i < size; i++) {
        int j = floor(inArray[i]*INT16_MAX);
        outArray[i] = j;
        if (outArray[i] > INT16_MAX)
            cout << outArray[i] << endl;
    }
    return outArray;
}

float * normalFloat (float* inArray, int size) {
    int maxValue = 0;
    for (int j =0; j < size; j++){
        if (abs(inArray[j]) > maxValue)
            maxValue = inArray[j];
    }
        
    float *outArray = new float[size];
    for (int i = 0; i < size; i++) {
        float f = inArray[i]/maxValue; //currMax;
        outArray[i] = f;
    }
    return outArray;
}


//This code courtesy of Leonard Manzara, CPSC 501, F16
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples;// * BYTES_PER_SAMPLE NOT NEEDED
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}