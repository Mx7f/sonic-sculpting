/**
  \file App.h

 */
#ifndef App_h
#define App_h

#include <G3D/G3DAll.h>

#ifdef G3D_WINDOWS
 // Compile RtAudio with windows direct sound support
#   define __WINDOWS_DS__
#endif
#include "RtAudio.h"
#include "chuck_fft.h"
#include "EWMAFrequency.h"

/** Global array where we dump raw audio data from RtAudio */
Array<float> g_currentAudioBuffer;

class App : public GApp {
protected:
    RtAudio m_rtAudio;

    /** EWMA of RMS */
    float m_smoothedRootMeanSquare;

    /** 3 different rates of exponentially-weighted moving averages of frequencies */
    EWMAFrequency m_fastMovingAverage;
    EWMAFrequency m_slowMovingAverage;
    EWMAFrequency m_glacialMovingAverage;


    /** Settings for RtAudio. We never need to change the defaults */
    struct AudioSettings {
      int numChannels;
      int sampleRate;
      RtAudioFormat rtAudioFormat;
      
      AudioSettings() :
          numChannels(1),
          sampleRate(48000),
          rtAudioFormat(RTAUDIO_FLOAT32) {}
    } m_audioSettings;

    /** How many slices of time to save */
    int m_maxSavedTimeSlices;

    /** CPU storage of raw samples */
    Array<float> m_cpuRawAudioData;
    /** CPU storage of fft samples */
    Array<complex> m_cpuFrequencyAudioData;

    /** GPU storage of raw samples */
    shared_ptr<Texture> m_rawAudioTexture;
    /** GPU storage of fft samples */
    shared_ptr<Texture> m_frequencyAudioTexture;


    /** Set all our audio textures on \param Args */
    void setAudioShaderArgs(Args& args);

    /** Do all of the interaction with RtAudio that we need to to set up realtime audio capture */
    void initializeAudio();

    /** Called from onInit */
    void makeGUI();

    /** Called once a frame to get the latest audio data and compute statistics such as RMS */
    void updateAudioData();


public:

    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;
    virtual void onGraphics3D(RenderDevice* rd, Array< shared_ptr<Surface> >& surface) override;
    virtual void onGraphics2D(RenderDevice* rd, Array< shared_ptr<Surface2D> >& surface2D) override;

    virtual bool onEvent(const GEvent& e) override;
    virtual void onCleanup() override;

    virtual void onUserInput(UserInput* ui) override;

    virtual void onSimulation(RealTime rdt, SimTime sdt, SimTime idt) override;

};

#endif
