#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

char* parseCompleteFile(char *currName);
int extractIntAtIndex(char *currFile, int index);

int main(int argc, char *argv[]) {
    
    if (argc != 3) {
        cout << "Please provide the correct arguements in the following order:\n"
        << "file1, file2";
    }
    
    char *f1Name = argv[1];
    char *f1Complete;
    int16_t *f1Data;
    int N;
    
    char *f2Name = argv[2];
    char *f2Complete;
    int16_t *f2Data;
    int M;

    int P;
   
    f1Complete = parseCompleteFile(f1Name);
    int f1ChunkSize = extractIntAtIndex(f1Complete, 4);
    int f1SubChunk1Size = extractIntAtIndex(f1Complete, 16);
    int f1SubChunk1SizeCorrected = f1SubChunk1Size - 16;
    memcpy(&N, f1Complete + (40 + f1SubChunk1SizeCorrected), 4); //file1 size
    f1Data = new int16_t[N];
    memcpy(f1Data, f1Complete + (44 + f1SubChunk1SizeCorrected), N); //no vals over INT16_MAX, less -INT16_MAX
    
    f2Complete = parseCompleteFile(f2Name);
    int f2ChunkSize = extractIntAtIndex(f2Complete, 4);
    int f2SubChunk1Size = extractIntAtIndex(f2Complete, 16);
    int f2SubChunk1SizeCorrected = f2SubChunk1Size - 16;
    memcpy(&M, f2Complete + (40 + f2SubChunk1SizeCorrected), 4); //file2 size
    f2Data = new int16_t[N];
    memcpy(f2Data, f2Complete + (44 + f2SubChunk1SizeCorrected), M); //no vals over INT16_MAX, less -INT16_MAX
    
    if (N == M) {
        cout << "data section size equal: " << N << endl;
        
        for (int i = 0; i < N; i++) {
            if(f1Data[i] != f2Data[i]) {
                cout << "data not equal at sample: " << i << endl;
                break;
            }
        }
    }
    else {
        cout << "data sections not equal with file1: " << N << "    and file2: " << M;
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
