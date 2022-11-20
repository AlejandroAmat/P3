/// @file

#include <iostream>
#include <fstream>
#include <string.h>
#include <errno.h>

#include "wavfile_mono.h"
#include "pitch_analyzer.h"

#include "docopt.h"

#define FRAME_LEN   0.030 /* 30 ms. */
#define FRAME_SHIFT 0.015 /* 15 ms. */

using namespace std;
using namespace upc;

static const char USAGE[] = R"(
get_pitch - Pitch Estimator 

Usage:
    get_pitch [options] <input-wav> <output-txt>
    get_pitch (-h | --help)
    get_pitch --version

Options:
    -p FLOAT, --umbralPot=FLOAT    umbral de potencia para estimacion pitch [default: -45]
    -z FLOAT, --umbralZCR=FLOAT    umbral de ZCR para estimacion pitch [default: 164]
    -n FLOAT, --umbralR1=FLOAT     umbral de R1 para estimacion pitch [default: 0.52]
    -u FLOAT, --umbralRmax==FLOAT  umbral de RMax para estimacion pitch  [default: 0.405]

Arguments:
    input-wav   Wave file with the audio signal
    output-txt  Output file: ASCII file with the result of the estimation:
                    - One line per frame with the estimated f0
                    - If considered unvoiced, f0 must be set to f0 = 0
)";






int main(int argc, const char *argv[]) {
	/// \TODO 
	///  Modify the program syntax and the call to **docopt()** in order to
	///  add options and arguments to the program.
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
        {argv + 1, argv + argc},	// array of arguments, without the program name
        true,    // show help if requested
        "2.0");  // version string

	std::string input_wav = args["<input-wav>"].asString();
	std::string output_txt = args["<output-txt>"].asString();
  float umbralPot=stof(args["--umbralPot"].asString());
  float umbralZCR=stof(args["--umbralZCR"].asString());
  float umbralR1=stof(args["--umbralR1"].asString());
  float umbralRmax=stof(args["--umbralRmax"].asString());




  // Read input sound file
  unsigned int rate;
  vector<float> x;
  if (readwav_mono(input_wav, rate, x) != 0) {
    cerr << "Error reading input file " << input_wav << " (" << strerror(errno) << ")\n";
    return -2;
  }

  int n_len = rate * FRAME_LEN;
  int n_shift = rate * FRAME_SHIFT;
  
  // Define analyzer
  PitchAnalyzer analyzer(n_len, rate, umbralPot, umbralZCR, umbralR1, umbralRmax, PitchAnalyzer::HAMMING, 50, 500);

  /// \TODO
  /// Preprocessing (Clipping) is done inside pitch_analyzer
  
  // Iterate for each frame and save values in f0 vector
  vector<float>::iterator iX;
  vector<float> f0;
  for (iX = x.begin(); iX + n_len < x.end(); iX = iX + n_shift) {
    float f = analyzer(iX, iX + n_len);
    f0.push_back(f);
  }

  /// \TODO
  /// Postprocess the estimation in order to supress errors. For instance, a median filter
  /// or time-warping may be used.



//Median filter of 3 coefficients
  vector<float> filter;
  vector<float> f0_filtered;

  f0_filtered.push_back(f0[0]); //We do not take into account f0(0) anf f0(size-1) as we apply a filter that gets f[n-1], f[n] and f[n+1];
  for(unsigned int l =1; l<f0.size()-1; l++){

    for(int k =-1; k<2; k++){
      filter.push_back(f0[l+k])  ;
      }

    sort(filter.begin(), filter.end());
    f0_filtered.push_back(filter[1]);
    filter.clear();

  }
  f0_filtered.push_back(f0[f0.size()-1]);

  // Write f0 contour into the output file
  ofstream os(output_txt);
  if (!os.good()) {
    cerr << "Error reading output file " << output_txt << " (" << strerror(errno) << ")\n";
    return -3;
  }

  os << 0 << '\n'; //pitch at t=0
  for (iX = f0_filtered.begin(); iX != f0_filtered.end(); ++iX) 
    os << *iX << '\n';
  os << 0 << '\n';//pitch at t=Dur

  return 0;
}
