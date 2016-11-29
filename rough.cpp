#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
using namespace std;

unsigned int* parseCompleteFile(char *currName);
unsigned int extractIntAtIndex(char *currFile, int index);

float* intToFloat (unsigned int* inArray, unsigned int size);
unsigned int* floatToInt (float* inArray, unsigned int size);


int main(int argc, char *argv[]) {
    
    if (argc != 3) {
        cout << "Please provide the correct arguements in the following order:\n"
        << "inputfile IRfile outputfile\n";
    }
    
    char *inName = argv[1];
    unsigned int *inputFileComplete;
    //unsigned int *x;
    float *x;
    int N;
    
    char *IRName = argv[2];
    char *IRFileComplete;
    //unsigned int *h;
    int *h;
    unsigned int M;
    
    char *outName = argv [3];
    char *outputFileComplete;
    //unsigned int *y;
    int *y;
    unsigned int P;
   
    inputFileComplete = parseCompleteFile(inName);
    /*unsigned int chunckSize = extractIntAtIndex(inputFileComplete, 4);
    unsigned int inSubChunk1Size = extractIntAtIndex(inputFileComplete, 16);
    unsigned int inSubChunk1SizeCorrected = inSubChunk1Size - 16;*/
    int chunkSize = inputFileComplete[4];
    cout << "chunkSize: " << chunkSize << endl;
    int inSubChunk1Size = inputFileComplete[16];
    int inSubChunk1SizeCorrected = inSubChunk1Size - 16;
    N = inputFileComplete[40 + inSubChunk1SizeCorrected]; //correct?
    unsigned int * inData = &inputFileComplete[44 + inSubChunk1SizeCorrected];
    cout << "N = 1069514 for DryDrums: " << N << endl; //seems to be working for Sitar.wav

    cout << "in spare if this not 16: " << inSubChunk1Size << endl;
    float *q = intToFloat(inData, N);
    for (int i = 0; i < N; i++) {
        cout << "test" << endl;
        cout << q[i] << " ";
    }

    
    /*IRFileComplete = parseCompleteFile(IRName);
    unsigned int IRSubChunk1Size = extractIntAtIndex(IRFileComplete, 16);
    unsigned int IRSubChunk1SizeCorrected = IRSubChunk1Size - 16;
    M = extractIntAtIndex(IRFileComplete, 40 + IRSubChunk1SizeCorrected); //correct?
    char * IRData = &IRFileComplete[44 + IRSubChunk1SizeCorrected];

    cout << "N = 1069514 for DryDrums: " << N << endl; //seems to be working for Sitar.wav
    cout << "chunkSize = 1069552 for DryDrums: " << chunckSize << endl;
    cout << "in spare if this not 16: " << inSubChunk1Size << endl;
    cout << "IR spare if this not 16: " << IRSubChunk1Size << endl;
    cout << "M = 1069514 for DryDrums: " << M << endl; //seems to be working for Sitar.wav
    x = inData;
    //memcpy(&x, inData, N);
    //memcpy(&h, IRData, M);
    float *xFloat = intToFloat(x, N);*/
    
    
    return 0;
}

float* intToFloat (unsigned int* inArray, unsigned int size) {
    float *outArray = new float[size];
    for (int i = 0; i < size; i++) {
        float f = inArray[i]/32786;
        cout << f << '\0' << endl;
        outArray[i] = f;
    }
    return outArray;
}

unsigned int* floatToInt (float* inArray, unsigned int size) {
    unsigned int *outArray = new unsigned int[size];
    for (int i = 0; i < size; i++) {
        int j = floor(inArray[i]*32787);
        outArray[i] = j;
    }
    return outArray;
}

unsigned int* parseCompleteFile(char *currName) {
    streampos size;
     unsigned int *memBlock;
    //char* memBlock;
    
    ifstream file (currName, ios::in|ios::binary|ios::ate);
    if (file.is_open()) {
        size = file.tellg();
        memBlock = new unsigned int [size];
        file.seekg(0, ios::beg);
        file.read((char*)memBlock,size);
        file.close();
        
        for (int i = 0; i < 50; i++) {
            cout << memBlock[i] << " ";
        }
    cout << endl;
    return memBlock;
    }
    else {
    cout << "File read error" << endl;
    }
}

unsigned int extractIntAtIndex(char *currFile, int index) {
    unsigned int returnVal;
    memcpy(&returnVal, currFile+index, sizeof(returnVal));
    return returnVal;
}