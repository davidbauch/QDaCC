#include <math.h>
#include <iostream>

using namespace std;

#define PI 3.14159265
int n;

int main( int argc, char **argv ) {
    cout << "Enter the size: ";
    cin >> n;
    double inputData[n][n];
    cout << "Enter the 2D elements ";
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            cin >> inputData[i][j];

    double realOut[n][n];
    double imagOut[n][n];
    double amplitudeOut[n][n];

    int height = n;
    int width = n;

    // Two outer loops iterate on output data.
    for ( int yWave = 0; yWave < height; yWave++ ) {
        for ( int xWave = 0; xWave < width; xWave++ ) {
            // Two inner loops iterate on input data.
            for ( int ySpace = 0; ySpace < height; ySpace++ ) {
                for ( int xSpace = 0; xSpace < width; xSpace++ ) {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += ( inputData[ySpace][xSpace] * cos( 2 * PI * ( ( 1.0 * xWave * xSpace / width ) + ( 1.0 * yWave * ySpace / height ) ) ) ) / sqrt( width * height );
                    imagOut[yWave][xWave] -= ( inputData[ySpace][xSpace] * sin( 2 * PI * ( ( 1.0 * xWave * xSpace / width ) + ( 1.0 * yWave * ySpace / height ) ) ) ) / sqrt( width * height );
                    amplitudeOut[yWave][xWave] = sqrt( realOut[yWave][xWave] * realOut[yWave][xWave] + imagOut[yWave][xWave] * imagOut[yWave][xWave] );
                }
                cout << realOut[yWave][xWave] << " + " << imagOut[yWave][xWave] << " i (" << amplitudeOut[yWave][xWave] << ")\n";
            }
        }
    }
}