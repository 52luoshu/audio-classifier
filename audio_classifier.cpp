#pragma once
#include "audio_classifier.h"
#include <string>

AudioClassifier::AudioClassifier() {
   mHamming = nullptr;
   mReferRet = nullptr;
   mStepCount = 0;
   mWindowLen = 0;
   mReferNum = 0;
   mReferCnt = 0;
   mReferRetInited = false;
   mRealRet = false;
}


AudioClassifier::~AudioClassifier() {
   if (mHamming) {
      delete[] mHamming;
   }
   if (mReferRet) {
      delete[] mReferRet;
   }
}

/*
* @brief: get the real time check result
* @return: true-contain music; false-dose not contain music
*/
bool AudioClassifier::getRealTimeCheckRet() {
   return mRealRet;
}

/*
* @brief: check whether current audio contains music
* @return: true-contain music; false-does not contain music
* @param audio: audio signal(1 sec)
* @param len: length of audio signal(number of samples)
*/
bool AudioClassifier::musicCheck(unsigned char* audio, unsigned short len) {
   algoParam params;
   float* emphasised = emphasis(audio, len, params.emphasisFactor);
   float** frames = frameSeg(emphasised, len, params.windowLen, params.stepLen);
   bool checkRet = false;
   bool lowEnergyCheckRet = musicCheckByLowEnergyFrameRate(frames, mWindowLen, mStepCount, params.lowEnergyFactors, params.lowEnergyCenters);
   if (lowEnergyCheckRet == true) {
      float** subBandCenters = (float**)array2dCreate(4, 3, sizeof(float));
      for (ushort i = 0; i < 4; i++) {
         for (ushort j = 0; j < 3; j++) {
            subBandCenters[i][j] = params.subBandCenters[i][j];
         }
      }
      bool subBandCheckRet = musicCheckByFft(emphasised, len, params.nfft, subBandCenters);
      array2dDelete((void**)subBandCenters, 4);

      if (subBandCheckRet == false) {
         checkRet = false;
      }
      else if (subBandCheckRet == true) {
         checkRet = true;
      }
   }
   else if (lowEnergyCheckRet == false) {
      checkRet = false;
   }
   mRealRet = checkRet;

   bool showRet = preRetRefer(checkRet, 30);
   if (emphasised) {
      delete[] emphasised;
   }
   array2dDelete((void**)frames, mStepCount);
   return showRet;
}

/*
* @brief: pre-emphasis
* @return: result after pre-emphasis
* @param datas: original audio signal
* @param len: length of original audio signal
* @param factor: factor for pre-emphasis
*/
float* AudioClassifier::emphasis(uchar* audio, ushort len, float factor) {
   float* empAudio = new float[len];
   empAudio[0] = (float)audio[0];
   for (ushort i = 0; i < len; i++) {
      empAudio[i] = (float)audio[i] - factor * audio[i - 1];
   }
   return empAudio;
}

/*
* @brief: add window, and divede audio to frames
* @return: audio frams
* @param audio: audio signal
* @param audioLen: length of window
* @param stepLen: lenth of window step
*/
float** AudioClassifier::frameSeg(float* audio, ushort audioLen, ushort windowLen, ushort stepLen) {
   if (mWindowLen != windowLen) {
      mWindowLen = windowLen;
      if (mHamming) {
         delete[] mHamming;
      }
      mHamming = hammingWindow(windowLen, 1.0);
   }
   mStepCount = (audioLen - windowLen) / stepLen + 1;

   float** frames = (float**)array2dCreate(mStepCount, windowLen, sizeof(float));
   for (uint i = 0; i < mStepCount; i++) {
      std::memcpy(frames[i], audio + (i*stepLen), windowLen * sizeof(float));
      hadmanMulti(frames[i], mHamming, windowLen);
   }
   return frames;
}

/*
* @brief: hamming window
* @return: pointer for hamming window
* @param len: window length
* @param amp: amplitude
*/
float* AudioClassifier::hammingWindow(ushort len, float amp) {
   float* hamming = new float[len];
   for (ushort i = 0; i < (len + 1) / 2; ++i) {
      hamming[i] = amp*(0.54 - 0.46*std::cos(2 * PI*i / (len - 1)));
      hamming[len - 1 - i] = hamming[i];
   }
   return hamming;
}


/*
* @brief: create 2d array
* @return: 2d array
* @param h: height of 2d array
* @param w: width of 2d array
* @param psize: size of data type
*/
void** AudioClassifier::array2dCreate(uint h, uint w, size_t psize) {
   void** array2d = new void*[h];
   for (uint i = 0; i < h; i++) {
      array2d[i] = (void*)new char[w*psize]();
   }
   return array2d;
}

/*
* @brief: delete ad array buffer
* @param buff: buffer for 2d array
* @param h: height of 2d array
*/
void AudioClassifier::array2dDelete(void** buff, uint h) {
   if (buff) {
      for (uint i = 0; i < h; i++) {
         if (buff[i]) {
            delete[]((uchar*)buff[i]);
         }
      }
      delete[] buff;
   }
}

/*
* @brief: hadman multiply
* @result: array1
* @param array1, array2: arrays for hadman multiply
* @param len: length of arrays
*/
void AudioClassifier::hadmanMulti(float* array1, float* array2, uint len) {
   for (uint i = 0; i < len; i++) {
      array1[i] *= array2[i];
   }
}

/*
* @brief: check audio type by low energy frame rate
* @return: retue-contain music; false-does not contain music
* @param frames: audio frames
* @param frmLen: length of frame
* @param frmNum: numbers of frames
* @param factors: factors for low energy frame rate
* @param centers: centers for low energy frame rate
*/
bool AudioClassifier::musicCheckByLowEnergyFrameRate(float** frames, ushort frmLen, uint frmNum, float factors[4], float centers[4][4]) {
   float rates[4] = { 0.0, 0.0, 0.0, 0.0 };
   for (ushort i = 0; i < 4; i++) {
      rates[i] = lowEnergyFrameRate(frames, frmLen, frmNum, factors[i]);
   }
   float dists[4][4] = {};
   float minDist[4] = { 1.0, 1.0, 1.0, 1.0 };
   for (ushort i = 0; i < 4; i++) {
      for (ushort j = 0; j < 4; j++) {
         dists[i][j] = distance(&rates[i], &centers[i][j], 1);
         minDist[i] = minDist[i] < dists[i][j] ? minDist[i] : dists[i][j];
      }
   }
   bool ret = false;
   if (minDist[3] != dists[3][0] && minDist[2] != dists[2][0] && minDist[1] != dists[1][0] && minDist[0] != dists[0][0]) {
      ret = true;
   }
   else {
      ret = false;
   }
   return ret;
}

/*
* @brief: calculate the average low energy frame rate of all frames
* @return: average low energy frame rate
* @param frames: audio frames
* @param frmLen: length of frame
* @param frmNum: numbers of frames
* @param factor: low energy frame rate factor
*/
float AudioClassifier::lowEnergyFrameRate(float** frames, ushort frmLen, uint frmNum, float factor) {
   double energy = 0.0, avgEnergy = 0.0;
   double* frameEnergy = new double[frmNum];
   std::memset(frameEnergy, 0, frmNum * sizeof(double));
   for (uint i = 0; i < frmNum; i++) {
      for (ushort j = 0; j < frmLen; j++) {
         frameEnergy[i] += std::pow(frames[i][j], 2);
      }
      energy += frameEnergy[i];
   }
   avgEnergy = factor * energy / frmNum;
   uint lowEnergyFrmCnt = 0;
   for (uint i = 0; i < frmNum; i++) {
      lowEnergyFrmCnt = frameEnergy[i] < avgEnergy ? (lowEnergyFrmCnt + 1) : lowEnergyFrmCnt;
   }
   if (frameEnergy) {
      delete[] frameEnergy;
   }
   float lowEnergyFrmRate = (float)lowEnergyFrmCnt / frmNum;
   return lowEnergyFrmRate;
}

/*
* @brief: calculate the distance of 2 arrays
* @return: distance of 2 arrays
* @param array1, array2: arrays for calculate
* @param len: length of arrays
*/
float AudioClassifier::distance(float* array1, float* array2, uint len) {
   float distance = 0.0;
   double tmp = 0.0;
   for (uint i = 0; i < len; i++) {
      tmp += std::pow((array1[i] - array2[i]), 2);
   }
   distance = std::sqrt(tmp);
   return distance;
}

/*
* @brief: check audio type by fft
* @return: true-contain music; false-does not contain music
* @param audio: audio signal
* @param len: length of audio signal
* @param nfft: nfft / 2 + 1 = sample rate
* @param centers: centers for sub-band energy ratio algotithm
*/
bool AudioClassifier::musicCheckByFft(float* audio, uint len, ushort nfft, float** centers) {
   float* fftRet = fft(audio, len, nfft);
   float* subBandRate = subBandEnergyRate(fftRet);
   ushort audioTypeNum = 4, subBandNum = 3;
   float* nearest = nearestCenter(subBandRate, centers, audioTypeNum, subBandNum);
   bool ret = false;
   if (nearest[0] == centers[0][0] && nearest[1] == centers[0][1] && nearest[2] == centers[0][2]) {
      ret = false;
   }
   else {
      ret = true;
   }
   if (fftRet) {
      delete[] fftRet;
   }
   if (subBandRate) {
      delete[] subBandRate;
   }
   return ret;
}

/*
* @brief: fast fourier transform
* @return: spectrum of signal
* @param frame: audio frame
* @param len: length of audio signal
* @param nfft: nfft / 2 + 1 = sample rate
*/
float* AudioClassifier::fft(float* frame, uint len, ushort nfft) {
   float* tmpDatas = new float[nfft];
   std::memset(tmpDatas, 0, nfft * sizeof(float));
   std::memcpy(tmpDatas, frame, len * sizeof(float));
   float* ret = new float[nfft / 2 + 1];
   uint bit = 1;
   while (nfft >> bit) { 
      ++bit; 
   }
   --bit;
   uint* rev = new uint[nfft];
   rev[0] = 0;
   for (int i = 0; i < nfft; i++) {
      rev[i] = (rev[i >> 1] >> 1) | (i & 1) << (bit - 1);
   }
   complex* dfft = new complex[nfft * 3];
   std::memset(dfft, 0, (nfft * 3) * sizeof(complex));
   for (uint i = 0; i < nfft; i++) {
      dfft[i].x = tmpDatas[rev[i]];
   }
   if (tmpDatas) {
      delete[] tmpDatas;
   }
   if (rev) {
      delete[] rev;
   }
   /* size of mid is half the length of the sequence to be merged */
   for (uint mid = 1; mid < nfft; mid *= 2) {
      /* the unit root, the coefficient 2 of PI has been reduced */
      complex temp(std::cos(PI / mid), std::sin(PI / mid));
      /* mid*2 is the length of the sequence to be merged,
         i is which bit is merged */
      for (int i = 0; i < nfft; i += mid * 2) {
         complex omega(1, 0);
         /* half scan */
         for (int j = 0; j < mid; j++, omega = omega * temp) {
            complex x = dfft[i + j], y = omega * dfft[i + j + mid];
            /* butterfly transform */
            dfft[i + j] = x + y;
            dfft[i + j + mid] = x - y;
         }
      }
   }
   for (int i = 0; i < nfft / 2 + 1; i++) {
      float x = dfft[i].x;
      float y = dfft[i].y;
      ret[i] = std::sqrt(x * x + y * y);
   }
   if (dfft) {
      delete[] dfft;
   }
   return ret;
}

/*
* @brief: calculate sub-band energy rate
* @return: sub-suband energy rate
* @param spectrum: spectrum (fft result)
*/
float* AudioClassifier::subBandEnergyRate(float* spectrum) {
   ushort subBandNum = 3;
   double* subBandEnergy = new double[subBandNum];
   std::memset(subBandEnergy, 0, subBandNum * sizeof(double));
   ushort zeroFreqPoint = 16000;
   ushort freqNode[7] = { zeroFreqPoint, zeroFreqPoint + 20, zeroFreqPoint + 64, zeroFreqPoint + 1200, zeroFreqPoint + 4000,zeroFreqPoint + 8000,zeroFreqPoint + 16000 };
   double energy = 0.0;
   for (ushort i = zeroFreqPoint; i < freqNode[6]; i++) {
      double sampleEnergy = std::pow(spectrum[i], 2);
      energy += sampleEnergy;
      if (i >= freqNode[1] && i < freqNode[2]) {         /* (20, 64) */
         subBandEnergy[0] += sampleEnergy;
      }
      else if (i >= freqNode[2] && i < freqNode[3]) {    /* (64, 1200) */
         subBandEnergy[1] += sampleEnergy;
      }
      else if (i >= freqNode[3] && i < freqNode[4]) {    /* (1200, 4000) */
         subBandEnergy[2] += sampleEnergy;
      }
   }
   float* subEnergyRate = new float[subBandNum];
   std::memset(subEnergyRate, 0, subBandNum * sizeof(float));
   for (ushort i = 0; i < subBandNum; i++) {
      subEnergyRate[i] = subBandEnergy[i] / energy;
   }
   if (subBandEnergy) {
      delete[] subBandEnergy;
   }
   return subEnergyRate;
}


/*
* @brief: calculate the nearest center
* @return: pointer of nearest center
* @param data: data to check
* @param centers: centers for check
* @param centerNum: numbers of center
* @param centerLen: length of center
*/
float* AudioClassifier::nearestCenter(float* data, float** centers, ushort centerNum, ushort centerLen) {
   float* distances = new float[centerNum];
   std::memset(distances, 0, centerNum * sizeof(float));
   for (ushort i = 0; i < centerNum; i++) {
      distances[i] = distance(data, centers[i], centerLen);
   }
   ushort nearestIdx = 0;
   float minDist = distances[0];
   for (ushort i = 0; i < centerNum; i++) {
      if (minDist >= distances[i]) {
         minDist = distances[i];
         nearestIdx = i;
      }
   }
   if (distances) {
      delete[] distances;
   }
   return centers[nearestIdx];
}

/*
* @brief: reference the previous result
* @return: result after reference the previous result
* @param currentRet: current result
* @param referNum: number to reference
*/
bool AudioClassifier::preRetRefer(bool currentRet, ushort referNum) {
   /* non-reference */
   if (referNum == 0) {
      return currentRet;
   }
   /* number to reference has reset, reinitialize the reference queue */
   if (mReferNum != referNum) {
      mReferNum = referNum;
      if (mReferRet) {
         delete[] mReferRet;
      }
      mReferRet = new bool[referNum];
      std::memset(mReferRet, 0, referNum * sizeof(bool));
      mReferRetInited = false;
      mReferCnt = 0;
   }
   bool ret = false;
   /* the reference queue has not completed a round of updates */
   if (mReferCnt < mReferNum) {
      ret = currentRet;
   }
   /* the reference queue has completed updates */
   else if (mReferCnt >= mReferNum) {
      /* loop */
      if (mReferCnt >= 60000) {
         mReferCnt = mReferNum;
      }
      float unitRate = 1.0 / (referNum + 1);
      float trueRate = 0.0, falseRate = 0.0;
      ushort repeat = 1;
      if (mReferRet[0] == true) { 
         trueRate += unitRate;
      }
      else if (mReferRet[0] == false) {
         falseRate += unitRate;
      }
      /* reference */
      for (ushort i = 1; i < mReferNum; i++) {
         if (mReferRet[i] == mReferRet[i - 1]) {
            repeat++;
         }
         else if (mReferRet[i] != mReferRet[i - 1]) {
            repeat = 1;
         }
         if (mReferRet[i] == true) {
            trueRate += repeat * unitRate;
         }
         else if (mReferRet[i] == false) {
            falseRate += repeat * unitRate;
         }
      }
      if (currentRet == true) {
         trueRate += unitRate;
      }
      else if (currentRet == false) {
         falseRate += unitRate;
      }
      ret = trueRate > falseRate ? true : false;
   }
   /* replace the reference queue */
   for (ushort i = 0; i < mReferNum - 1; i++) {
      mReferRet[i] = mReferRet[i + 1];
   }
   mReferRet[mReferNum - 1] = currentRet;
   mReferCnt++;
   return ret;
}