#include "Synthesizer.h"

shared_ptr<Synthesizer> Synthesizer::global = shared_ptr<Synthesizer>(new Synthesizer());

void Synthesizer::queueSound(shared_ptr<AudioSample> audioSample, int delay) {
    mutex.lock(); {
      
        m_sounds.append(SoundInstance(audioSample, -delay, format("%d", m_soundCount)));
	m_totalHack.append(m_sounds[m_sounds.size()-1].m_framebuffer);
	m_totalHack2.append(m_sounds[m_sounds.size()-1].m_rawAudioTexture, m_sounds[m_sounds.size()-1].m_displayTexture);
	++m_soundCount;
    } mutex.unlock();
}

void Synthesizer::flushHackQueue() {
    for (int i = m_totalHack.size() - 1; i >= 0; --i) {
        if (m_totalHack[i].unique()) {
            m_totalHack.remove(i);
        }
    }
    for (int i = m_totalHack2.size() - 1; i >= 0; --i) {
        if (m_totalHack2[i].unique()) {
            m_totalHack2.remove(i);
        }
    }
}

void Synthesizer::getSoundInstances(Array<SoundInstance>& soundInstances) {
    mutex.lock(); {
        for (auto s : m_sounds) {
            soundInstances.append(s);
        }
    } mutex.unlock();
}

void Synthesizer::synthesize(Array<Sample>& samples) {
    mutex.lock(); {
        int maxIndex = m_sounds.size() - 1;
        for (int i = maxIndex; i >= 0; --i) {
            if (m_sounds[i].play(samples)) {
                m_sounds.remove(i);
            }
        }
        sampleCount += double(samples.size());
    } mutex.unlock();
}

