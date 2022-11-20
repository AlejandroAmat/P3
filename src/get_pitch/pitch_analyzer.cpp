/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"
using namespace std;
static int numero= 0;
/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {
    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l]
      r[l] =0;
      for(unsigned int n = l; n < x.size(); n++){
          r[l] += x[n]*x[n-l];
      }
      r[l] /= x.size();
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \TODO Implement the Hamming window
     
      if(frameLen%2==0){ //if framelen is even
        unsigned int half = frameLen/2;

        for(unsigned int n =0; n<half; n++){ //Hamming Window for first half
        window[n] = 0.5 * (1-cos(2*M_PI*(n+1)/frameLen+1));
        }

        unsigned int idx = half-1;
        for(unsigned int n=half; n<frameLen; n++){ //Symmentric window for the second half
          window[n]= window[idx];
          idx--;
        }


      }
      else{

          unsigned int half = (frameLen + 1 )/2;

        for(unsigned int n =0; n<half; n++){ //Hamming Window for first half
        window[n] = 0.5 * (1-cos(2*M_PI*(n+1)/frameLen+1));
        }

        int idx = half-2;
        for(unsigned int n=half; n<frameLen; n++){ //Symmentric window for the second half
          window[n]= window[idx];
          idx--;
        }



      }

    //  for(unsigned int n =0; n<frameLen; n++){
    //   window[n] = 0.53836 - 0.46164*(cos((2*M_PI*n)/(frameLen-1)));
    //  }
      
      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, float zcr) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    numero++;
    cout<<"("<<numero+1<<") pot:"<<pot<<", zcr:"<<zcr<<", r1norm:"<<r1norm<<", maxnorm:"<<rmaxnorm<< endl;
      //if((rmaxnorm >= 0.72 || r1norm >0.9) && pot>-15.7  ) return false;
      // if(pot<-10) return true;
    // if((r1norm>0.86 && ((pot>-36 ||  rmaxnorm>0.6) && zcr<120))) return false;
    bool unvoiced = true;
    if(r1norm>0.72668 && zcr <149 && pot>-45 &&  rmaxnorm >0.34) unvoiced=false;;

    return unvoiced;
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    // Signal Normalitzation
    // max_signal = *std::max_element(x.begin(), x.end());
    // for(auto it = x.begin(); it<x.end(); it++){
    //   *it /= max_signal;
    // }

    //Center Clipping
     float max_signal = *std::max_element(x.begin(), x.end());
    float Cl = 0.1133*max_signal;
    for(auto it = x.begin(); it<x.end(); it++){
        if(abs(*it)<Cl) *it=1e-10;
        else if(*it > Cl) *it+=1.317*Cl;
        else if(*it < -Cl) *it-=1.3117*Cl;
    }

   

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);

    //Compute correlation
    autocorrelation(x, r);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.

    // for (iR=iRMax = r.begin()+npitch_min; iR< r.begin() + npitch_max ;iR++ ){
    //   if(*iR > *iRMax) iRMax = iR;
    // }
     for (iR=iRMax = r.begin()+npitch_min; iR< r.begin() + npitch_max ;iR++ ){
      if(*iR > *iRMax) iRMax = iR;
    }

    unsigned int lag = iRMax - r.begin();

    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    //CRUCES POR 0

       double zcr =0;
       for(unsigned int i =1; i < x.size();i++){
      if((x[i]>0 && x[i-1]<0) || (x[i]<0 && x[i-1]>0) ) //falta considerar caso 0
   	    zcr++;
      }
      

    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0], zcr))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }
}

