#ifndef Synthesizer_h
#define Synthesizer_h

#include <mutex>
#include "SoundInstance.h"

class Synthesizer {
public:
    
    static shared_ptr<Synthesizer> global;
private:
    std::mutex mutex;
    Array<SoundInstance> m_sounds;
    int m_soundCount;
    double sampleCount;
    double lastSampleCount;
public:
    Synthesizer() : sampleCount(0.0), lastSampleCount(0.0), m_soundCount(0) {}
    void queueSound(shared_ptr<AudioSample> audioSample, int delay = 0);

    double currentSampleCount() const {
        return sampleCount;
    }

    int totalSoundCount() const {
        return m_soundCount;
    }

    double tick() {
        double result = sampleCount - lastSampleCount;
        lastSampleCount = sampleCount;
        return result;
    }

    void getSoundInstances(Array<SoundInstance>& soundInstances);

    void synthesize(Array<float>& samples);
};
#endif