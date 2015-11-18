#ifndef AudioSample_h
#define AudioSample_h
#include <G3D/G3DAll.h>
#include "util.h"
struct AudioSample {
    Array<Sample> buffer;
    /** In Hz */
    int sampleRate;
    int sampleCount() {
        return buffer.size();
    }
    static shared_ptr<AudioSample> createSine(int sampleRate, double frequency, int sampleCountDuration, float fadeOutProportion) {
        shared_ptr<AudioSample> s(new AudioSample());
        s->sampleRate = sampleRate;
        s->buffer.resize(sampleCountDuration);
        int fadeOutBeginSample = sampleCountDuration * fadeOutProportion;
        float volume = 0.1f; // Allow ten simulataneous samples without clipping
        for (int i = 0; i < sampleCountDuration; ++i) {
            float sinValue = sin(2.0f * pif() * float(i) * float(frequency) / sampleRate);
            float fade = clamp(1.0f - (float(i - fadeOutBeginSample) / (sampleCountDuration - fadeOutBeginSample)), 0.0f, 1.0f);
            s->buffer[i] = sinValue * volume * fade;
        }
        return s;
    }
	static shared_ptr<AudioSample> createEmpty(int sampleRate) {
		shared_ptr<AudioSample> s(new AudioSample());
		s->sampleRate = sampleRate;
		return s;
	}
};
#endif