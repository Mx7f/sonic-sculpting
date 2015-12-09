#include <G3D/G3DAll.h>
#ifndef util_h
#define util_h

typedef float Sample;

//https://en.wikipedia.org/wiki/Piano_key_frequencies
static double pianoKeyNumberToFrequency(int index) {
  alwaysAssertM(index > 0 && index <= 88, "Piano Key index needs to be in [1,88]");
  return pow(2.0, (double(index)-49.0)/12.0) * 440.0;
}

G3D_DECLARE_ENUM_CLASS(PianoKey,C,D_b,D,E_b,E,F,F_s,G,G_s,A,B_b,B);
static int pianoKeyNameAndOctaveToIndex(PianoKey p, int octave) {
  int index = (12*octave) + p - 8;
  return index;
}

static double getFrequencyFromKey(PianoKey p, int octave) {
  return pianoKeyNumberToFrequency(pianoKeyNameAndOctaveToIndex(p, octave));
}

#define METERS_PER_SAMPLE_WINDOW 0.1f

#endif
