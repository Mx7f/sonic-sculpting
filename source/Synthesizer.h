#ifndef Synthesizer_h
#define Synthesizer_h

#include <mutex>
#include "SoundInstance.h"

class Synthesizer {
public:
    
    static shared_ptr<Synthesizer> global;
private:
    // Needed to keep driver from crashing OS X
    // Actual problem: without these in place, the last reference to texture/framebuffer of a SoundInstance
    // will be deleted on the sound thread, not the GL thread, which causes the driver to flip out.
    // Proper fix would be to move those sounds to a delete queue, will do after presentation deadline.
    // Issue #44 https://github.com/Mx7f/sonic-sculpting/issues/44
    Array<shared_ptr<Framebuffer>>  m_totalHack;
    Array<shared_ptr<Texture>>      m_totalHack2;


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
    /** See documentation for m_totalHack as for why this is necessary in the short term. 
        https://github.com/Mx7f/sonic-sculpting/issues/44 
    */
    void flushHackQueue();

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
