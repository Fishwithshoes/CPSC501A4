#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

char* parseCompleteFile(char *currName);
unsigned int extractIntAtIndex(char *currFile, int index);


int main(int argc, char *argv[]) {
    
    if (argc != 3) {
        cout << "Please provide the correct arguements in the following order:\n"
        << "inputfile IRfile outputfile\n";
    }
    
    char *inName = argv[1];
    char *inputFileComplete;
    unsigned int *x;
    unsigned int N;
    
    /*char *IRName = argv[2];
    char *IRFileComplete;
    char *h;
    long M;
    
    char *outName = argv [3];
    char *outputFileComplete;
    char *y;
    long P;*/
   
    inputFileComplete = parseCompleteFile(inName);
    unsigned int chunckSize = extractIntAtIndex(inputFileComplete, 4);
    cout << "chunkSize = 1069552 for DryDrums: " << chunckSize << endl;
    unsigned int subChunk1Size = extractIntAtIndex(inputFileComplete, 16);
    unsigned int subChunk1SizeCorrected = subChunk1Size - 16;
    cout << "spare if this not 16: " << subChunk1Size << endl;
    N = extractIntAtIndex(inputFileComplete, 40 + subChunk1SizeCorrected); //correct?
    char * data = &inputFileComplete[44 + subChunk1SizeCorrected];
    cout << "N = 1069514 for DryDrums: " << N << endl; //seems to be working for Sitar.wav
    memcpy(&x, data, N);
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
        
        /*for (int i = 0; i < size; i++) {
            cout << (int) memBlock[i] << " ";
        }*/
    cout << endl;
    return memBlock;
    }
    else {
    cout << "File read error" << endl;
    }
}

unsigned int extractIntAtIndex(char *currFile, int index) {
    unsigned int returnVal;
    char * temp = &currFile[index];
    memcpy(&returnVal, temp, sizeof(returnVal));
    return returnVal;
}