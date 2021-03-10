/*****************************************************************************************************
* @file: audio_classifier.h 
* @brief: check wether the current audio contain music
* @notice: sample rate of the input audio must be 32000, all the params fit to 1 sec duration input
* @author: liushenghua
* @email: 1203096438@qq.com
* @time: 2021-03-10
******************************************************************************************************/
#pragma once
#ifndef AUDIO_CLASSIFIER_H
#define AUDIO_CLASSIFIER_H
#include <cmath>
#include <cstdint>

typedef unsigned char uchar;     /* uint8_t */
typedef unsigned short ushort;   /* uint16_t */
typedef unsigned int uint;       /* uint32_t */

struct complex {
   float x, y;
   complex(float xx=0.0,float yy=0.0):x(xx),y(yy){}
   ~complex(){}
   complex operator + (complex a) { return complex(x + a.x, y + a.y); }
   complex operator - (complex a) { return complex(x - a.x, y - a.y); }
   complex operator * (complex a) { return complex(x * a.x - y * a.y, x * a.y + y * a.x); }
};


const float PI = std::acos(-1);

struct algoParam {
   /* pre-emphasis factor */
   float emphasisFactor = 0.97; 
   /* params for framing */
   ushort sampleRate = 32000; /* frequence */
   ushort windowLen = 1024;   /* 32 ms */
   ushort stepLen = 512;      /* 16 ms */
   /* centers for low energy frame rate algorithm */
   float lowEnergyCenters[4][4]= {
      /* {speak,           music,               mix,                  noise} */
      {0.393327961287313,  0.5054023977272730,  0.504118400000000000, 0.492891335299903000},    /* 1.00 */
      {0.185053230317164,  0.0266393568181818,  0.021431416585365900, 0.001813431907571290},    /* 0.75 */
      {0.124097763432836,  0.0223547045454545,  0.000319872195121951, 0.000636717059980334},    /* 0.50 */
      {0.053385740000000,	0.0143442600000000,	0.000239904000000000, 0.000620598000000000}     /* 0.25 */
   };
   float lowEnergyFactors[4] = { 1.00, 0.75, 0.50, 0.25 };
   /* centes for sub-band energy ratio algotithm */
   float subBandCenters[4][3]={ 
      /* {(20,64),   (64,1200),  (1200,4000)} */
      {0.01556305,	0.7921705,	0.117853},   /* speak */
      {0.02127115,	0.7563860,	0.140953},   /* music */
      {0.01849620,	0.6927620,	0.173328},   /* mix */
      {0.01242895,	0.6725795,	0.195659}    /* noise */
   };
   ushort fftLen = 16000;  /* spectrum length for fft, freqnence / 2 */
   ushort nfft = 63998;    /* nfft / 2 + 1 = fftLen * 2 = sampleRate */
};

class AudioClassifier {
private:
   ushort mWindowLen;
   float* mHamming;
   ushort mStepCount;
   bool* mReferRet;
   ushort mReferNum;
   bool mReferRetInited;
   ushort mReferCnt;
   bool mRealRet;
private:
   /*
   * @brief: pre-emphasis
   * @return: result after pre-emphasis
   * @param datas: original audio signal
   * @param len: length of original audio signal
   * @param factor: factor for pre-emphasis
   */
   float* emphasis(uchar* audio, ushort len, float factor);
   /*
   * @brief: add window, and divede audio to frames
   * @return: audio frams
   * @param audio: audio signal
   * @param audioLen: length of window
   * @param stepLen: lenth of window step
   */
   float** frameSeg(float* audio, ushort audioLen, ushort windowLen, ushort stepLen);
   /*
   * @brief: hamming window
   * @return: pointer for hamming window
   * @param len: window length
   * @param amp: amplitude
   */
   float* hammingWindow(ushort len, float amp);
   /*
   * @brief: create 2d array
   * @return: 2d array
   * @param h: height of 2d array
   * @param w: width of 2d array
   * @param psize: size of data type
   */
   void** array2dCreate(uint h, uint w, size_t psize);
   /*
   * @brief: delete ad array buffer
   * @param buff: buffer for 2d array
   * @param h: height of 2d array
   */
   void array2dDelete(void** buff, uint h);
   /*
   * @brief: hadman multiply
   * @result: array1
   * @param array1, array2: arrays for hadman multiply
   * @param len: length of arrays
   */
   void hadmanMulti(float* array1, float* array2, uint len);
   /*
   * @brief: check audio type by low energy frame rate
   * @return: retue-contain music; false-does not contain music
   * @param frames: audio frames
   * @param frmLen: length of frame
   * @param frmNum: numbers of frames
   * @param factors: factors for low energy frame rate
   * @param centers: centers for low energy frame rate
   */
   bool musicCheckByLowEnergyFrameRate(float** frames, ushort frmLen, uint frmNum, float factors[4], float centers[4][4]);
   /*
   * @brief: calculate the average low energy frame rate of all frames
   * @return: average low energy frame rate
   * @param frames: audio frames
   * @param frmLen: length of frame
   * @param frmNum: numbers of frames
   * @param factor: low energy frame rate factor
   */
   float lowEnergyFrameRate(float** frames, ushort frmLen, uint frmNum, float factor);
   /*
   * @brief: calculate the distance of 2 arrays
   * @return: distance of 2 arrays
   * @param array1, array2: arrays for calculate
   * @param len: length of arrays
   */
   float distance(float* array1, float* array2, uint len);
   /*
   * @brief: check audio type by fft
   * @return: true-contain music; false-does not contain music
   * @param audio: audio signal
   * @param len: length of audio signal
   * @param nfft: nfft / 2 + 1 = sample rate
   * @param centers: centers for sub-band energy ratio algotithm
   */
   bool musicCheckByFft(float* audio, uint len, ushort nfft, float** centers);
   /*
   * @brief: calculate sub-band energy rate
   * @return: sub-suband energy rate
   * @param spectrum: spectrum (fft result)
   */
   float* subBandEnergyRate(float* spectrum);
   /*
   * @brief: fast fourier transform
   * @return: spectrum of signal
   * @param frame: audio frame
   * @param len: length of audio signal
   * @param nfft: nfft / 2 + 1 = sample rate
   */
   float* fft(float* frame, uint len, ushort nfft);
   /*
   * @brief: calculate the nearest center
   * @return: pointer of nearest center
   * @param data: data to check
   * @param centers: centers for check
   * @param centerNum: numbers of center
   * @param centerLen: length of center
   */
   float* nearestCenter(float* data, float** centers, ushort centerNum, ushort centerLen);
   /*
   * @brief: reference the previous result
   * @return: result after reference the previous result
   * @param currentRet: current result
   * @param referNum: number to reference
   */
   bool preRetRefer(bool currentRet, ushort referNum);
public:
   AudioClassifier();
   ~AudioClassifier();
   /*
   * @brief: check whether current audio contains music
   * @return: true-contain music; false-does not contain music
   * @param audio: audio signal(1 sec)
   * @param len: length of audio signal(number of samples)
   */
   bool musicCheck(unsigned char* audio, unsigned short len);
   /*
   * @brief: get the real time check result
   * @return: true-contain music; false-dose not contain music
   */
   bool getRealTimeCheckRet();
};

#endif // !AUDIO_CLASSIFIER_H
