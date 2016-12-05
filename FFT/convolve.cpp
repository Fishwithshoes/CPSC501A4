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
double* intToDouble(int16_t* inArray, int size);
int16_t* doubleToInt (double* inArray, int size);
double * normalDouble (double* inArray, int size);
void convolve(double x[], double h[], int L, double y[], int P);
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

int formatN(int N, int M);
void four1(double data[], int nn, int isign);
void complexMulti(double X[], double H[], double Y[], int L);

#define PI                3.14159265358979

#define SIZE       8
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

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
    double *raw_x;
    double *x;
    int N;
    
    char *IRName = argv[2];
    char *IRFileComplete;
    int16_t *IRData;
    double *raw_h;
    double *h;
    int M;
    
    char *outName = argv[3];
    double *y;
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
    raw_x = new double[N];
    raw_x = intToDouble(inData, N); // no vals over 1.0, less -1.0
    
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
    raw_h = new double[M];
    raw_h = intToDouble(IRData, M); //no vals over 1.0, less -1.0

    
    int L = formatN(N,M);
    cout << "L: " << L << endl;
    
    x = new double[L*2];
    h = new double[L*2];

    for (int i; i < L*2; i++) {
        x[i] = 0.0;
        h[i] = 0.0;
    }

    memcpy(x, raw_x, N);
    memcpy(h, raw_h, M);
     
    P = N + M - 1;
    y = new double[L*2];
    
    convolve(x, h, L, y, P);
    
    outData = new short[L];
    
    outData = doubleToInt(y, L);
    
    FILE *outputFileStream = fopen(outName, "wb");
    if (outputFileStream == NULL) {
        fprintf(stderr, "File %s cannot be opened for writing\n", outName);
    }
    
    writeWaveFileHeader(1,L,44100.0, outputFileStream);
    for (int i = 0; i < L; i++) {
        fwriteShortLSB(outData[i], outputFileStream);
    }
    
    //;cout << "P: " << P <<endl;
    
    fclose(outputFileStream);
    
    return 0;
}

void convolve(double x[], double h[], int L, double y[], int P) {
    four1(x-1, L, 1); //now X
    four1(h-1, L, 1); //now H
    for (int k = 0, i =0; k < L; k++, i+=2){    //scale output
        x[i] /= (double) L;
        x[i+1] /= (double) L;
        h[i] /= (double) L;
        h[i+1] /= (double) L;
        
        x[k] = sqrt(x[i] * x[i] + x[i+1] * x[i+1]);     //calc amplitude
        h[k] = sqrt(h[i] * h[i] + h[i+1] * h[i+1]);
        }
                        for (int z; z < L; z++) {
            cout << x[z] << " ";
    }
    
    for (int k = 1, j = L-1; k < L/2; k++, j--) {           //combine pos/neg amplitudes
        x[k] = x[k] + x[j];
        y[k] = y[k] + y[j];
    }
    
    complexMulti(x,h,y,L);
    four1(y-1, L, -1);
    

}

 void complexMulti(double X[], double H[], double Y[], int L) {
    for (int k = 0, i =0; k < L; k++, i+=2){
        Y[i] = ((X[i] * H[i]) - (X[i+1] * H[i+1]));
        Y[i+1] = ((X[i+1] * H[i]) + (X[i] * H[i+1]));
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

double * intToDouble (int16_t* inArray, int size) {
    double *outArray = new double[size];
    for (int i = 0; i < size; i++) {
        double f = (double) inArray[i]/(INT16_MAX + 1); //currMax;
        outArray[i] = f;
    }
    return outArray;
}

int16_t* doubleToInt (double* inArray, int size) {
    int16_t *outArray = new int16_t[size];
    for (int i = 0; i < size; i++) {
        int16_t j = ceil(inArray[i]*INT16_MAX);
        outArray[i] = j;
    }
    return outArray;
}

double * normalDouble (double* inArray, int size) {
    double maxValue = 0.0;
    for (int j =0; j < size; j++){
        if (abs(inArray[j]) > maxValue)
            maxValue = abs(inArray[j]);
    }
        
    double *outArray = new double[size];
    for (int i = 0; i < size; i++) {
        double f = inArray[i]/maxValue; //currMax;
        outArray[i] = f;
    }
    return outArray;
}
int formatN(int N, int M) {
    int val, newN;
    if (N > M) {
        val = N;
    }
    else {
        val = M;
    }
    
    for (int i = 1; i < 20; i++) {
        newN = pow(2, i);
        if (val <= newN) {
            return newN;
        }
    }
}

//All following code courtesy of Leonard Manzara, CPSC 501, F16
//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

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