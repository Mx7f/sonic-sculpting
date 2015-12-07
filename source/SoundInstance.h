#ifndef SoundInstance_h
#define SoundInstance_h
#include<G3D/G3DAll.h>
#include "AudioSample.h"

struct SoundInstance {
    shared_ptr<AudioSample> audioSample;
    int currentPosition;
    String name;
    /** Returns true if finished */
    bool play(Array<float>& buffer) {
        if (currentPosition == audioSample->sampleCount()) {
            return true;
        }
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
    void initialize();
    SoundInstance() {}
    SoundInstance(const shared_ptr<AudioSample>& audioSample, int currentPosition, const String& name) :
        audioSample(audioSample), currentPosition(currentPosition), name(name) {
        initialize();
    }
    shared_ptr<Texture> displayTexture() const {
        return m_displayTexture;
    }
    shared_ptr<Framebuffer> m_framebuffer;
    shared_ptr<Texture> m_rawAudioTexture;
    shared_ptr<Texture> m_displayTexture;

private:

};
#endif
