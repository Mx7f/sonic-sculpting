#ifndef Synthesizer_h
#define Synthesizer_h
#include<G3D/G3DAll.h>
#include "AudioSample.h"
#include <mutex>
struct SoundInstance {
    shared_ptr<AudioSample> audioSample;
    int currentPosition;
    /** Returns true if finished */
    bool play(Array<float>& buffer) {
        // TODO: double check this
        for (int i = 0; i < buffer.size(); ++i) {
            if (currentPosition >= 0) {
                buffer[i] += audioSample->buffer[currentPosition];
            }
            ++currentPosition;
            if (currentPosition == audioSample->sampleCount()) {
                return true;
            }
        }
        return false;
    }
    SoundInstance() {}
    SoundInstance(const shared_ptr<AudioSample>& audioSample, int currentPosition) :
        audioSample(audioSample), currentPosition(currentPosition) {}
};

class Synthesizer {
public:
    
    static shared_ptr<Synthesizer> global;
private:
    std::mutex mutex;
    Array<SoundInstance> m_sounds;
    double sampleCount;
    double lastSampleCount;
public:
    Synthesizer() : sampleCount(0.0), lastSampleCount(0.0) {}
    void queueSound(shared_ptr<AudioSample> audioSample, int delay = 0);

    double currentSampleCount() const {
        return sampleCount;
    }

    double tick() {
        double result = sampleCount - lastSampleCount;
        lastSampleCount = sampleCount;
        return result;
    }

    void synthesize(Array<float>& samples);
};
#endif