#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
using namespace std;

char* parseCompleteFile(char *currName);
unsigned int extractIntAtIndex(char *currFile, unsigned int index);
float* intToFloat (int* inArray, int size);

int main(int argc, char *argv[]) {
    
    if (argc != 3) {
        cout << "Please provide the correct arguements in the following order:\n"
        << "inputfile IRfile outputfile\n";
    }
    
    char *inName = argv[1];
    char *inputFileComplete ;
    int *inData;
    float *x;
    int N;
   
    inputFileComplete = parseCompleteFile(inName);
    unsigned int chunkSize;
    unsigned int chunkSize = extractIntAtIndex(inputFileComplete, 4);
    unsigned int inSubChunk1Size = extractIntAtIndex(inputFileComplete, 16);
    unsigned int inSubChunk1SizeCorrected = inSubChunk1Size - 16;
    memcpy(&N, inputFileComplete + (40 + inSubChunk1SizeCorrected), 4);
    cout << "N: " << endl;
    inData = new int[N];
    memcpy(inData, inputFileComplete + (44 + inSubChunk1SizeCorrected), N);
    x = new float[N];
    x = intToFloat(inData, N);
    for (int i = 0; i < N; i++) {
        cout << x[i] << " ";
    }

    
    return 0;
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

unsigned int extractIntAtIndex(char *currFile, unsigned int index) {
    unsigned int returnVal;
    memcpy(&returnVal, currFile+index, sizeof(returnVal));
    return returnVal;
}

float * intToFloat (int* inArray, int size) {
    int currMax;
    for (int j = 0; j < size; j++) {
        if (inArray[j] > currMax)
            currMax = inArray[j];
    }
    
    float *outArray = new float[size];
    for (int i = 0; i < size; i++) {
        float f = (float) inArray[i]/1073741824;
        outArray[i] = f;
    }
    return outArray;
}